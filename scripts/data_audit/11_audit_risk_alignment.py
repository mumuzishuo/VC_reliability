#!/usr/bin/env python3
"""Retrospective audit: absolute MSE ranking versus failure relative to zero.

This reads a previously frozen source-only risk table and previously evaluated
HepG2 task errors. It must not be used to select or refit a prospective risk rule.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr


RULES = ("source_heterogeneity", "low_source_support", "predicted_rms")
COVERAGE = 0.8


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def table(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    require(bool(rows), f"Empty table: {path}")
    return rows


def classification(score, failure):
    """Average-tie AUROC and noninterpolated average precision; high = risk."""
    npos = int(failure.sum())
    nneg = len(failure) - npos
    if npos == 0 or nneg == 0:
        return None, None
    ranks = rankdata(score, method="average")
    auc = (ranks[failure].sum() - npos * (npos + 1) / 2) / (npos * nneg)
    ix = np.argsort(-score, kind="stable")
    sorted_scores, sorted_y = score[ix], failure[ix]
    end = np.r_[np.flatnonzero(sorted_scores[1:] != sorted_scores[:-1]), len(ix) - 1]
    tp = np.cumsum(sorted_y)[end]
    ap = np.sum(np.diff(np.r_[0, tp]) * tp / (end + 1)) / npos
    return float(auc), float(ap)


def check_freeze(protocol_path, baseline_code, prediction_h5, risk_table):
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    require(protocol["target"] == "HepG2", "Expected the original HepG2 pilot")
    require(protocol["primary_risk"] == RULES[0], "Original primary risk changed")
    for filename, expected in (
        (baseline_code, protocol["script_sha256"]),
        (prediction_h5, protocol["prediction_sha256"]),
        (risk_table, protocol["risk_table_sha256"]),
    ):
        require(sha256(filename) == expected, f"Frozen file hash mismatch: {filename.name}")
    source_inputs = protocol["source_inputs"]
    require(set(source_inputs) == {"K562", "RPE1", "Jurkat"}, "Unexpected source files")
    for env, item in source_inputs.items():
        require(sha256(Path(item["path"])) == item["sha256"],
                f"Source response changed after freezing: {env}")
    return protocol, len(source_inputs)


def audit(per_task, risks, protocol):
    targets = [row["perturbation"] for row in per_task]
    require(len(targets) == len(set(targets)) == protocol["n_perturbations"],
            "Invalid or changed task cohort")
    require({r["target_environment"] for r in per_task} == {"HepG2"},
            "Task errors are not from HepG2")
    risk_by_gene = {row["perturbation"]: row for row in risks}
    require(len(risk_by_gene) == len(risks) and set(risk_by_gene) == set(targets),
            "Risk and evaluated task labels differ")
    frozen = [risk_by_gene[gene] for gene in targets]
    require({r["target_environment"] for r in frozen} == {"HepG2"},
            "Risk rows have the wrong target environment")
    tie = np.array([int(row["tie_key"]) for row in frozen], dtype=np.int64)
    n = len(targets)
    require(set(tie.tolist()) == set(range(n)), "Frozen tie keys are not a permutation")
    e = np.array([float(row["source_mean_mse"]) for row in per_task])
    z = np.array([float(row["zero_response_mse"]) for row in per_task])
    require(np.isfinite(e).all() and np.isfinite(z).all() and
            (e >= 0).all() and (z >= 0).all(), "Invalid task errors")
    d = e - z
    failure = d > 0
    k = math.ceil(COVERAGE * n)
    primary = np.array([float(row[RULES[0]]) for row in frozen])
    primary_order = np.lexsort((tie, primary))
    flagged = np.zeros(n, dtype=np.int64)
    # The original review flag rounds the excluded 20% UP; the coverage
    # evaluation rounds the retained 80% UP. For N=1264 these differ by one.
    flagged[primary_order[n - math.ceil(0.2 * n):]] = 1
    require(np.array_equal(flagged, [int(row["primary_review_flag"]) for row in frozen]),
            "Frozen primary review flags differ from the source-only rank")
    output = []
    for name in RULES:
        score = np.array([float(row[name]) for row in frozen])
        ref = np.array([float(row[name]) for row in per_task])
        require(np.isfinite(score).all() and np.allclose(score, ref, rtol=0, atol=1e-12),
                f"Frozen source risk score changed: {name}")
        ix = np.lexsort((tie, score))[:k]
        policy = z.copy()
        policy[ix] = e[ix]
        auc, ap = classification(score, failure)
        output.append({
            "rule": name,
            "score_direction": "higher means higher original risk",
            "n_tasks": n,
            "n_retained": k,
            "actual_coverage": k / n,
            "all_mean_mse": float(e.mean()),
            "all_zero_mean_mse": float(z.mean()),
            "all_mean_D": float(d.mean()),
            "all_median_D": float(np.median(d)),
            "failure_prevalence_D_gt0": float(failure.mean()),
            "spearman_score_vs_mse": float(spearmanr(score, e)[0]),
            "spearman_score_vs_D": float(spearmanr(score, d)[0]),
            "auroc_D_gt0": auc,
            "auprc_AP_D_gt0": ap,
            "retained_mean_mse": float(e[ix].mean()),
            "retained_zero_mean_mse": float(z[ix].mean()),
            "retained_mean_D": float(d[ix].mean()),
            "retained_failure_rate": float(failure[ix].mean()),
            "retained_vs_all_mse_improvement": float(1 - e[ix].mean() / e.mean()),
            "policy_mse_retained_prediction_else_zero": float(policy.mean()),
            "policy_minus_always_prediction": float(policy.mean() - e.mean()),
            "policy_minus_always_zero": float(policy.mean() - z.mean()),
            "random_same_coverage_expected_policy_mse": float(z.mean() + k / n * d.mean()),
        })
    return output, {"n_tasks": n, "n_primary_review_flags": int(flagged.sum()),
                    "n_excluded_at_80_percent_coverage": n - k,
                    "n_failures": int(failure.sum()),
                    "n_better_than_zero": int((d < 0).sum()), "n_ties": int((d == 0).sum())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-task", required=True, type=Path)
    parser.add_argument("--risk-table", required=True, type=Path)
    parser.add_argument("--frozen-protocol", required=True, type=Path)
    parser.add_argument("--baseline-code", required=True, type=Path)
    parser.add_argument("--prediction-h5", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output_dir.exists(), "Output directory exists; use a new version")
    protocol, source_count = check_freeze(
        args.frozen_protocol, args.baseline_code, args.prediction_h5, args.risk_table
    )
    records, cohort = audit(table(args.per_task), table(args.risk_table), protocol)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    with (args.output_dir / "risk_alignment_summary.tsv").open(
        "x", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    report = {
        "scope": "retrospective HepG2 pilot; do not use for prospective rule selection",
        "failure_definition": "MSE(source_mean) - MSE(zero) > 0",
        "frozen_protocol_sha256": sha256(args.frozen_protocol),
        "n_frozen_source_files_sha256_verified": source_count,
        "source_only_scores_match_frozen_risk_table": True,
        "prediction_and_baseline_script_sha256_verified": True,
        "cohort": cohort,
        "risk_rules": list(RULES),
        "coverage": COVERAGE,
        "limitation": "hash consistency does not establish independent preregistration or external validity",
    }
    (args.output_dir / "risk_alignment_summary.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print("FREEZE_HASHES_OK", json.dumps(cohort), flush=True)
    for row in records:
        print("RISK_ALIGNMENT", row["rule"], json.dumps({
            key: row[key] for key in ("spearman_score_vs_mse", "spearman_score_vs_D",
                                  "auroc_D_gt0", "retained_failure_rate",
                                  "policy_mse_retained_prediction_else_zero")
        }), flush=True)


if __name__ == "__main__":
    main()
