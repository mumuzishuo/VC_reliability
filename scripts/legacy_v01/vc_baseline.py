#!/usr/bin/env python3
"""First descriptive baseline for the VC reliability pilot (Python 3.9+).

Run in this order on an allocated CPU node:
  python vc_baseline.py predict
  python vc_baseline.py evaluate
Dependencies: numpy, scipy, h5py; matplotlib is optional (plots only).
Inputs: ROOT/processed/v01/{sources,evaluation_only}/*_responses.h5.
Outputs: ROOT/results/baseline_v01/ (never overwrite existing results).

predict reads ONLY K562, RPE1 and Jurkat outcomes. For each observed perturbation
p, it predicts the arithmetic mean response across these three environments,
weighting environments equally. This is a transparent context-transfer baseline,
not GEARS, a learned virtual cell, or a predictor for unseen perturbation genes.
It does not use HepG2 controls or outcomes. It freezes model predictions, source
features, tie-breaking keys, protocol and file digests BEFORE evaluate can read
HepG2 outcomes. Directory separation alone is not an access-security mechanism.

Primary risk rule H_p = mean_{s,g}[(delta_{s,p,g} - mean_s delta_{s,p,g})**2].
Comparators (directions fixed in advance): 1 / minimum source cell count;
and RMS predicted response. Higher scores always mean higher proposed risk.
These are rankings, NOT failure probabilities or proven confidence estimates.
The primary review flag marks the highest-risk ceil(20% * N) predictions.

evaluate computes per-perturbation MSE across ALL fixed shared expression genes,
response Pearson correlation, and a zero-response baseline. Risk-coverage orders
predictions from low to high risk; AURC is the arithmetic mean of the retained
mean MSE over all k/N, k=1..N. At fixed coverage c, retain ceil(c*N), report actual
coverage, and compare with expected random selection (overall mean MSE). Ties
use source-stage seeded keys independent of target errors. Oracle sorting uses
the target errors and is an evaluation-only lower reference, never a method.

This run is exploratory and descriptive: no confidence intervals, p-values,
learned risk model, observation-noise correction, DEG selection or Systema
perturbation-specific metric. Target response-strength tertiles and cell-count
strata are post-reveal diagnostics, not risk inputs. Report all fixed rules;
do not select/flip a rule after seeing HepG2 and call it independent validation.
Tasks share controls, genes and pathways; the N perturbations are not N
independent biological replicates. A simple mean beating zero, or a positive
risk association, is not enough to claim general biological reliability.
HepG2/Jurkat external checksums remain unverified in the existing data record.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import spearmanr

SOURCES = ("K562", "RPE1", "Jurkat")
RULES = ("source_heterogeneity", "low_source_support", "predicted_rms")
SEED = 20261101
COVERAGES = (0.5, 0.8, 0.9)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def strings(dataset):
    return np.asarray(dataset.asstr()[:])


def read_response(path, role, environment):
    with h5py.File(path, "r") as f:
        require(f.attrs["role"] == role, "Incorrect data role: " + str(path))
        require(f.attrs["environment"] == environment, "Incorrect environment")
        perturbations, gene_ids = strings(f["perturbation"]), strings(f["gene_id"])
        delta = np.asarray(f["delta"][:], dtype=np.float64)
        n = np.asarray(f["n_cells"][:], dtype=np.int64)
        low = np.asarray(f["low_control_fraction"][:], dtype=np.float64)
        info = {"manifest_sha256": str(f.attrs["manifest_sha256"]),
                "normalization": str(f.attrs["normalization"])}
    require(delta.shape == (len(perturbations), len(gene_ids)), "Response shape mismatch")
    require(len(perturbations) == len(set(perturbations)), "Duplicate perturbations")
    require(len(gene_ids) == len(set(gene_ids)), "Duplicate expression gene IDs")
    require(n.shape == low.shape == (len(perturbations),), "Cell metadata shape mismatch")
    require(bool(np.isfinite(delta).all()), "Nonfinite responses")
    require(bool((n > 0).all()), "Invalid perturbation cell counts")
    require(bool((np.isfinite(low) & (low >= 0) & (low <= 1)).all()), "Invalid low-control fraction")
    return perturbations, gene_ids, delta, n, low, info


def align(reference_p, reference_g, p, g, delta, n, low):
    require(set(reference_p) == set(p), "Perturbation sets differ")
    require(set(reference_g) == set(g), "Expression gene sets differ")
    p_map, g_map = {s: i for i, s in enumerate(p)}, {s: i for i, s in enumerate(g)}
    pi = np.asarray([p_map[s] for s in reference_p])
    gi = np.asarray([g_map[s] for s in reference_g])
    return delta[np.ix_(pi, gi)], n[pi], low[pi]


def commit(temp, final):
    # An existing final file is never silently replaced.
    os.link(temp, final)
    temp.unlink()


def write_table(path, rows):
    require(bool(rows), "Refusing an empty table: " + str(path))
    temp = path.with_name(path.name + ".partial")
    with temp.open("x", newline="", encoding="utf-8") as h:
        writer = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("NA" if isinstance(v, (float, np.floating)) and
                                  not np.isfinite(v) else v) for k, v in row.items()})
    commit(temp, path)


def write_predictions(path, p, g, prediction, features, tie_key, info):
    temp = path.with_name(path.name + ".partial")
    with h5py.File(temp, "x") as f:
        for name, values in [("perturbation", p), ("gene_id", g)]:
            f.create_dataset(name, data=np.asarray(values, dtype=object),
                             dtype=h5py.string_dtype("utf-8"))
        f.create_dataset("source_mean", data=prediction, compression="gzip", compression_opts=1)
        f.create_dataset("tie_key", data=tie_key)
        group = f.create_group("risk_features")
        for name, values in features.items():
            group.create_dataset(name, data=values)
        for name, value in info.items():
            f.attrs[name] = value
    commit(temp, path)


def risk_order(score, tie_key):
    require(bool(np.isfinite(score).all()), "Nonfinite risk score")
    return np.lexsort((tie_key, score))


def predict(root, out):
    require(not out.exists(), "Output directory already exists; do not repeat predict. "
            "Use evaluate or inspect the existing run.")
    inputs, arrays, counts, low_fractions = {}, [], [], []
    ref_p = ref_g = ref_info = None
    for env in SOURCES:
        path = root / "processed" / "v01" / "sources" / (env + "_responses.h5")
        p, g, delta, n, low, info = read_response(path, "sources", env)
        if ref_p is None:
            ref_p, ref_g, ref_info = p, g, info
        else:
            require(info == ref_info, "Sources use different preprocessing manifests")
            delta, n, low = align(ref_p, ref_g, p, g, delta, n, low)
        arrays.append(delta)
        counts.append(n)
        low_fractions.append(low)
        inputs[env] = {"path": str(path.resolve()), "sha256": digest(path)}
    stack = np.stack(arrays)
    prediction = stack.mean(axis=0)
    features = {
        "source_heterogeneity": np.mean((stack - prediction[None, :, :]) ** 2, axis=(0, 2)),
        "low_source_support": 1.0 / np.min(np.stack(counts), axis=0),
        "predicted_rms": np.sqrt(np.mean(prediction ** 2, axis=1)),
    }
    tie_key = np.random.default_rng(SEED).permutation(len(ref_p))
    order = risk_order(features[RULES[0]], tie_key)
    rank = np.empty(len(ref_p), dtype=int)
    rank[order] = np.arange(len(ref_p))
    flagged = rank >= len(ref_p) - math.ceil(0.2 * len(ref_p))
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for i, p in enumerate(ref_p):
        row = {"perturbation": p, "target_environment": "HepG2",
               **{r: float(features[r][i]) for r in RULES},
               "primary_risk_percentile": float(100 * (rank[i] + 1) / len(ref_p)),
               "primary_review_flag": int(flagged[i]), "tie_key": int(tie_key[i])}
        row.update({"n_source_" + env: int(counts[j][i]) for j, env in enumerate(SOURCES)})
        row.update({"low_control_fraction_source_" + env: float(low_fractions[j][i])
                    for j, env in enumerate(SOURCES)})
        rows.append(row)
    write_table(out / "risk_predictions.tsv", rows)
    write_predictions(out / "predictions.h5", ref_p, ref_g, prediction, features, tie_key, ref_info)
    protocol = {
        "protocol": "source-mean-baseline-v0.1", "target": "HepG2", "sources": list(SOURCES),
        "n_perturbations": len(ref_p), "n_expression_genes": len(ref_g),
        "primary_risk": RULES[0], "all_rules": list(RULES), "higher_score_means": "higher risk",
        "source_weighting": "equal environment weights; same perturbation in each source",
        "primary_loss": "mean squared error across all fixed expression genes",
        "AURC": "mean over k=1..N of mean MSE among the k lowest-risk predictions",
        "fixed_coverages": list(COVERAGES), "coverage_rounding": "ceil(c * N)",
        "random_reference": "exact expected retained MSE = all-task mean MSE",
        "tie_break": "seeded permutation frozen before reading target outcomes", "seed": SEED,
        "secondary_diagnostics": ["true response RMS tertiles (evaluation only)",
                                  "target cell count <50 vs >=50 (evaluation only)",
                                  "target low-control cells absent vs present (evaluation only)"],
        "source_inputs": inputs, "preprocessing": ref_info,
        "script_sha256": digest(Path(__file__)),
        "prediction_sha256": digest(out / "predictions.h5"),
        "risk_table_sha256": digest(out / "risk_predictions.tsv"),
        "limitations": ["developmental descriptive pilot, not independent validation",
                        "no uncertainty intervals or hypothesis tests",
                        "no correction for observation noise or shared-control dependence",
                        "no perturbation-specific Systema metric yet",
                        "risk scores are not error probabilities",
                        "target counts are eligibility/evaluation metadata, not risk features",
                        "HepG2/Jurkat external raw-file checksums remain unverified"]}
    with (out / "frozen_protocol.json").open("x", encoding="utf-8") as h:
        json.dump(protocol, h, indent=2)
        h.write("\n")
    print(f"PREDICT_OK: {len(ref_p)} x {len(ref_g)}; target outcomes NOT read", flush=True)
    print("Frozen predictions and risks:", out, flush=True)


def ratio_gain(retained_error, all_error):
    return float(1 - retained_error / all_error) if all_error > 0 else float("nan")


def spearman(score, error):
    if len(score) < 3 or np.all(score == score[0]) or np.all(error == error[0]):
        return float("nan")
    # Tuple access works with both scipy 1.9.1 on the HPC and newer versions.
    return float(spearmanr(score, error)[0])


def curve(error, order):
    return np.cumsum(error[order], dtype=np.float64) / np.arange(1, len(error) + 1)


def risk_summary(error, features, tie_key, stratum="all"):
    rows, curves = [], {}
    all_error = float(np.mean(error))
    oracle = float(np.mean(curve(error, risk_order(error, tie_key))))
    for name in RULES:
        score = features[name]
        rc = curve(error, risk_order(score, tie_key))
        curves[name] = rc
        aurc = float(np.mean(rc))
        row = {"stratum": stratum, "risk_rule": name, "n_tasks": len(error),
               "score_constant": int(np.all(score == score[0])),
               "spearman_risk_vs_mse": spearman(score, error),
               "all_mean_mse": all_error, "AURC": aurc,
               "expected_random_AURC": all_error, "oracle_AURC": oracle,
               "AURC_relative_improvement": ratio_gain(aurc, all_error)}
        for coverage in COVERAGES:
            k = math.ceil(coverage * len(error))
            key = str(int(100 * coverage))
            row.update({"retained_n_" + key: k, "actual_coverage_" + key: k / len(error),
                        "retained_mean_mse_" + key: float(rc[k - 1]),
                        "relative_improvement_" + key: ratio_gain(rc[k - 1], all_error)})
        rows.append(row)
    return rows, curves


def pearson_rows(a, b):
    ac, bc = a - a.mean(axis=1, keepdims=True), b - b.mean(axis=1, keepdims=True)
    numerator = np.einsum("ij,ij->i", ac, bc)
    denominator = np.sqrt(np.einsum("ij,ij->i", ac, ac) * np.einsum("ij,ij->i", bc, bc))
    return np.divide(numerator, denominator, out=np.full(len(a), np.nan), where=denominator > 0)


def plots(out, features, errors, risk_curves, zero_error, truth_rms):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("PLOTS_SKIPPED: matplotlib missing; numeric outputs are complete", flush=True)
        return
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    coverage = np.arange(1, len(errors) + 1) / len(errors)
    labels = ("Source heterogeneity", "Low source support", "Predicted RMS")
    for rule, label in zip(RULES, labels):
        axes[0].plot(coverage, risk_curves[rule], label=label)
    axes[0].axhline(errors.mean(), color="black", linestyle="--", label="Expected random")
    axes[0].set(xlabel="Coverage (fraction retained)", ylabel="Mean response MSE",
                title="Risk–coverage: source-mean baseline")
    axes[0].legend(fontsize=7)
    dots = axes[1].scatter(np.log1p(features[RULES[0]]), np.log1p(errors), c=truth_rms,
                           s=8, alpha=.55, cmap="viridis")
    axes[1].set(xlabel="log1p(source heterogeneity)", ylabel="log1p(response MSE)",
                title="Primary risk rule vs observed error")
    fig.colorbar(dots, ax=axes[1], label="Observed response RMS (evaluation only)")
    axes[2].bar(["Zero response", "Source mean"], [zero_error.mean(), errors.mean()],
                 color=["#929baa", "#3278b7"])
    axes[2].set(ylabel="All-task mean MSE", title="Prediction baseline comparison")
    fig.savefig(out / "baseline_overview.png", dpi=180)
    fig.savefig(out / "baseline_overview.pdf")
    plt.close(fig)


def evaluate(root, out):
    protocol_path = out / "frozen_protocol.json"
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    require(not (out / "evaluation_summary.json").exists(), "Evaluation already completed")
    require(protocol["script_sha256"] == digest(Path(__file__)), "Script changed after freezing")
    for filename, key in [("predictions.h5", "prediction_sha256"),
                          ("risk_predictions.tsv", "risk_table_sha256")]:
        require(digest(out / filename) == protocol[key], "Frozen artifact changed: " + filename)
    for env, item in protocol["source_inputs"].items():
        require(digest(item["path"]) == item["sha256"], "Source changed after prediction: " + env)
    with h5py.File(out / "predictions.h5", "r") as f:
        ref_p, ref_g = strings(f["perturbation"]), strings(f["gene_id"])
        prediction = f["source_mean"][:].astype(np.float64)
        tie_key = f["tie_key"][:]
        features = {r: f["risk_features"][r][:] for r in RULES}
    target = root / "processed" / "v01" / "evaluation_only" / "HepG2_responses.h5"
    p, g, truth, n, low, info = read_response(target, "evaluation_only", "HepG2")
    require(info == protocol["preprocessing"], "Target preprocessing differs from sources")
    truth, n, low = align(ref_p, ref_g, p, g, truth, n, low)
    error = np.mean((prediction - truth) ** 2, axis=1)
    zero_error = np.mean(truth ** 2, axis=1)
    truth_rms = np.sqrt(zero_error)
    correlations = pearson_rows(prediction, truth)
    baseline_rows = []
    for model, losses in [("zero_response", zero_error), ("source_mean", error)]:
        finite = correlations[np.isfinite(correlations)] if model == "source_mean" else []
        baseline_rows.append({"model": model, "n_tasks": len(error), "n_genes": len(ref_g),
                              "mean_mse": float(losses.mean()), "median_mse": float(np.median(losses)),
                              "mean_response_pearson": float(np.mean(finite)) if len(finite) else float("nan"),
                              "defined_pearson_tasks": len(finite),
                              "relative_mse_improvement_vs_zero": ratio_gain(float(losses.mean()), float(zero_error.mean()))})
    write_table(out / "baseline_metrics.tsv", baseline_rows)
    per_task = []
    for i, p in enumerate(ref_p):
        per_task.append({"perturbation": p, "target_environment": "HepG2",
                         **{r: float(features[r][i]) for r in RULES},
                         "source_mean_mse": float(error[i]), "zero_response_mse": float(zero_error[i]),
                         "source_mean_pearson": float(correlations[i]),
                         "zero_mse_minus_source_mean_mse": float(zero_error[i] - error[i]),
                         "observed_rms_evaluation_only": float(truth_rms[i]),
                         "target_n_cells_evaluation_only": int(n[i]),
                         "target_low_control_fraction_evaluation_only": float(low[i])})
    write_table(out / "prediction_level_results.tsv", per_task)
    risk_rows, risk_curves = risk_summary(error, features, tie_key)
    write_table(out / "risk_metrics.tsv", risk_rows)
    curve_rows = []
    oracle_curve = curve(error, risk_order(error, tie_key))
    for i in range(len(error)):
        curve_rows.append({"retained_n": i + 1, "coverage": (i + 1) / len(error),
                           **{r: float(risk_curves[r][i]) for r in RULES},
                           "expected_random": float(error.mean()),
                           "oracle_evaluation_only": float(oracle_curve[i])})
    write_table(out / "risk_coverage.tsv", curve_rows)
    tertile = np.digitize(truth_rms, np.quantile(truth_rms, [1 / 3, 2 / 3]), right=True)
    strata = [("observed_rms_" + name, tertile == i)
              for i, name in enumerate(["low", "middle", "high"])]
    strata += [("target_cells_lt50", n < 50), ("target_cells_ge50", n >= 50),
               ("target_low_control_absent", low == 0), ("target_low_control_present", low > 0)]
    diagnostics = []
    for name, mask in strata:
        if mask.sum() >= 3:
            rows, _ = risk_summary(error[mask], {r: features[r][mask] for r in RULES},
                                   tie_key[mask], name)
            diagnostics.extend(rows)
    if diagnostics:
        write_table(out / "stratified_diagnostics.tsv", diagnostics)
    plots(out, features, error, risk_curves, zero_error, truth_rms)
    summary = {"status": "EVALUATE_OK", "n_tasks": len(error), "n_genes": len(ref_g),
               "source_mean_mse": float(error.mean()), "zero_response_mse": float(zero_error.mean()),
               "primary_risk_spearman": risk_rows[0]["spearman_risk_vs_mse"],
               "primary_risk_improvement_80": risk_rows[0]["relative_improvement_80"],
               "source_mean_fraction_better_than_zero": float(np.mean(error < zero_error)),
               "target_file_sha256": digest(target), "frozen_protocol_sha256": digest(protocol_path),
               "warning": "Descriptive pilot only; no independent validation or observation-noise correction"}
    # Write JSON with null for any undefined scalar, never non-standard NaN tokens.
    summary = {k: (None if isinstance(v, float) and not np.isfinite(v) else v)
               for k, v in summary.items()}
    with (out / "evaluation_summary.json").open("x", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, allow_nan=False)
        f.write("\n")
    print("EVALUATE_OK:", json.dumps(summary, allow_nan=False), flush=True)
    print("Tables and optional plots:", out, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=["predict", "evaluate"])
    parser.add_argument("--root", type=Path, default=Path.home() / "VC_reliability")
    parser.add_argument("--output", type=Path, help="Default ROOT/results/baseline_v01")
    args = parser.parse_args()
    out = args.output or args.root / "results" / "baseline_v01"
    (predict if args.mode == "predict" else evaluate)(args.root, out)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        sys.exit(1)
