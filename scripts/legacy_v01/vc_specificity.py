#!/usr/bin/env python3
"""Retrospective specificity audit of the three existing internal folds.

Primary retrieval uses cosine of control-relative response vectors.
Secondary retrieval subtracts each matrix's mean across all 1264 perturbations.
The target centroid is evaluation-only. All models use the same candidate/query
set within each target and representation. Responses of zero norm in truth or
any experimental half are excluded and documented. Zero-norm queries receive
chance credit via a uniform tie. Ties within 1e-12 cosine are averaged over
uniform random tie breaking, including MRR (not reciprocal of mean rank).
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import vc_baseline as bb
import vc_internal as vi
import vc_preprocess as pp

TIE_TOL = 1e-12
NORM_TOL = 1e-12
REPRESENTATIONS = ("control_relative_cosine", "perturbation_centered_cosine")
MODELS = ("zero", "paired_raw", "common_raw", "paired_calibrated", "common_calibrated")
require = bb.require


def scale_from_sources(a, b):
    power = np.mean(a * a) + np.mean(b * b)
    return float(np.clip(2 * np.mean(a * b) / power, 0, 1)) if power > 0 else 0.0


def represent(x, name):
    if name == REPRESENTATIONS[0]:
        return x
    # Avoid numerical cancellation creating a direction for constant queries.
    if np.array_equal(x, np.broadcast_to(x[0], x.shape)):
        return np.zeros_like(x)
    return x - x.mean(axis=0, keepdims=True)


def cosine_scores(query, candidates):
    qn = np.linalg.norm(query, axis=1)
    cn = np.linalg.norm(candidates, axis=1)
    require(np.isfinite(query).all() and np.isfinite(candidates).all(), "Nonfinite retrieval input")
    require((cn > NORM_TOL).all(), "Zero-norm candidates must be excluded by the common mask")
    q = np.divide(query, qn[:, None], out=np.zeros_like(query), where=qn[:, None] > NORM_TOL)
    return np.clip(q @ (candidates / cn[:, None]).T, -1, 1), qn <= NORM_TOL


def retrieval(query, candidates):
    require(query.shape == candidates.shape and len(query) > 1, "Retrieval needs paired labels and N>1")
    scores, undefined = cosine_scores(query, candidates)
    correct = np.diag(scores)
    better = np.sum(scores > correct[:, None] + TIE_TOL, axis=1)
    tied = np.sum(np.abs(scores - correct[:, None]) <= TIE_TOL, axis=1)
    require((tied > 0).all(), "Correct candidate missing from tie group")
    n = len(query)
    harmonic = np.r_[0., np.cumsum(1 / np.arange(1, n + 1))]
    top1 = np.clip(1 - better, 0, tied) / tied
    top10 = np.clip(min(10, n) - better, 0, tied) / tied
    mrr = (harmonic[better + tied] - harmonic[better]) / tied
    percentile = (better + (tied - 1) / 2) / (n - 1)
    rows = dict(top1=top1, top10=top10, mrr=mrr, rank_percentile=percentile,
                n_better=better, n_tied=tied, zero_norm_query=undefined)
    summary = {key: float(rows[key].mean()) for key in ("top1", "top10", "mrr", "rank_percentile")}
    summary.update(n_candidates=n, n_queries=n, n_zero_norm_queries=int(undefined.sum()),
                   chance_top1=1/n, chance_top10=min(10, n)/n,
                   chance_mrr=float(harmonic[-1]/n), chance_rank_percentile=0.5)
    return summary, rows


def freeze(root, out, manifest):
    old = root / "results/internal_v01"
    original = vi.verify_protocol(root, old)
    require(not out.exists(), f"Output exists: {out}. Select a NEW --output.")
    for env in vi.ENVS:
        require(vi.stage_done(old, "split_" + env), f"{env}: missing verified split files")
    out.mkdir(parents=True)
    hashes = {}
    for target in vi.ENVS:
        source = [env for env in vi.ENVS if env != target]
        a, b = [vi.load_source(root, env, manifest)[0] for env in source]
        paired = (a + b) / 2
        common = paired.mean(axis=0)
        alpha = scale_from_sources(a, b)
        beta = scale_from_sources(a.mean(axis=0), b.mean(axis=0))
        # Each coefficient is the minimizer of the two directional source
        # transfer losses for its own predictor family, constrained to [0,1].
        path = out / f"{target}_frozen.h5"
        pp.write_h5(path, dict(perturbation=manifest["perturbations"],
                              gene_id=manifest["gene_ids"], paired=paired, common=common),
                    dict(target=target, sources=",".join(source),
                         alpha_paired=alpha, alpha_common=beta))
        hashes[path.name] = bb.digest(path)
        print(f"FREEZE_OK {target}: paired alpha={alpha:.6f}; common alpha={beta:.6f}", flush=True)
    vi.json_write(out / "protocol.json", {
        "analysis": "retrospective internal specificity audit, not independent validation",
        "script_sha256": bb.digest(Path(__file__)),
        "previous_protocol_sha256": bb.digest(old / "frozen_protocol.json"),
        "split_repeats": original["split_repeats"],
        "frozen_prediction_sha256": hashes,
        "norm_threshold": NORM_TOL, "cosine_tie_tolerance": TIE_TOL,
        "primary_metric": "paired_calibrated MSE minus common_calibrated MSE on ALL inherited tasks",
        "secondary_metrics": "raw cosine retrieval; perturbation-centered cosine retrieval; half-split reference",
        "common_response": "equal-weight mean over source perturbations, then equal-weight mean over source environments",
        "calibration": "symmetric single-source transfer within sources, constrained scalar [0,1], applied to two-source mean",
        "evaluation_cohort": "manifest cohort selected using all four environments; no target response used to set coefficients",
        "candidate_mask": "same queries/candidates across models and splits within representation; nonzero full truth and every experimental half",
        "target_centering": "evaluation-only target centroid across full inherited perturbation panel; never a risk input",
        "tie_credit": "expected top-k and reciprocal rank under uniform tie breaking",
        "halves": "existing perturbation-stratified cells and GEM-stratified controls; two directions per repeat are dependent",
        "inference": "descriptive per-environment comparisons, no independent-screen confidence intervals",
    })


def evaluate(root, out, manifest):
    internal = root / "results/internal_v01"
    protocol = json.loads((out / "protocol.json").read_text())
    require(protocol["script_sha256"] == bb.digest(Path(__file__)), "Audit code changed")
    require(protocol["previous_protocol_sha256"] == bb.digest(internal / "frozen_protocol.json"),
            "Previous protocol changed")
    for name, digest in protocol["frozen_prediction_sha256"].items():
        require(bb.digest(out / name) == digest, "Frozen prediction changed")
    model_rows, task_rows, summary_rows, rank_rows, exclusion_rows, batch_rows = [], [], [], [], [], []
    for target in vi.ENVS:
        with h5py.File(out / f"{target}_frozen.h5", "r") as f:
            require(np.array_equal(bb.strings(f["perturbation"]), manifest["perturbations"])
                    and np.array_equal(bb.strings(f["gene_id"]), manifest["gene_ids"]),
                    "Frozen label order mismatch")
            pred = f["paired"][:]
            common = np.broadcast_to(f["common"][:], pred.shape)
            alpha, beta = f.attrs["alpha_paired"], f.attrs["alpha_common"]
        # Prediction files and both source-only coefficients already exist.
        truth, _ = vi.load_source(root, target, manifest)
        predictions = dict(zero=np.zeros_like(pred), paired_raw=pred, common_raw=common,
                           paired_calibrated=alpha*pred, common_calibrated=beta*common)
        errors = {name: np.mean((p - truth)**2, axis=1) for name, p in predictions.items()}
        for name in MODELS:
            err, zero = errors[name], errors["zero"]
            baseline = errors["common_calibrated"]
            model_rows.append(dict(target=target, model=name, n_tasks=len(truth),
                                   alpha_paired=alpha, alpha_common=beta, mse=float(err.mean()),
                                   zero_mse=float(zero.mean()), gain_vs_zero=float(1-err.mean()/zero.mean()) if zero.mean() else np.nan,
                                   fraction_better_zero=float(np.mean(err < zero)),
                                   mse_minus_common_calibrated=float((err-baseline).mean()),
                                   fraction_better_common_calibrated=float(np.mean(err < baseline))))
        for i, p in enumerate(manifest["perturbations"]):
            task_rows.append(dict(target=target, perturbation=p,
                                  **{name+"_mse": errors[name][i] for name in MODELS},
                                  paired_minus_common_calibrated=errors["paired_calibrated"][i]-errors["common_calibrated"][i]))
        with h5py.File(internal / f"{target}_halves.h5", "r") as f:
            for key, values in (("perturbation", manifest["perturbations"]), ("gene_id", manifest["gene_ids"])):
                require(np.array_equal(bb.strings(f[key]), values), "Half-split order mismatch")
            require(f.attrs["protocol_sha256"] == protocol["previous_protocol_sha256"],
                    "Half-split protocol mismatch")
            ds = f["delta"]
            require(ds.shape[:2] == (protocol["split_repeats"], 2), "Unexpected half dimensions")
            # Batch composition audit runs on all tasks, independently of retrieval exclusions.
            nb = f["n_cells_by_batch"][:].astype(np.float64)
            weights = nb / nb.sum(axis=3, keepdims=True)
            tv = 0.5 * np.abs(weights[:, 0] - weights[:, 1]).sum(axis=2)
            for r in range(len(tv)):
                batch_rows.append(dict(target=target, repeat=r+1,
                                       median_half_GEM_weight_TV=float(np.median(tv[r])),
                                       p95_half_GEM_weight_TV=float(np.quantile(tv[r], .95)),
                                       max_half_GEM_weight_TV=float(tv[r].max())))
            for representation in REPRESENTATIONS:
                actual = represent(truth, representation)
                keep = np.linalg.norm(actual, axis=1) > NORM_TOL
                require(np.isfinite(actual).all(), "Nonfinite full target")
                # A fixed eligible panel is shared by every model and half comparison.
                for r in range(ds.shape[0]):
                    for side in range(2):
                        half = represent(ds[r, side].astype(np.float64), representation)
                        require(np.isfinite(half).all(), "Nonfinite experimental half")
                        keep &= np.linalg.norm(half, axis=1) > NORM_TOL
                require(keep.sum() > 1, "Fewer than two eligible retrieval candidates")
                ids = np.asarray(manifest["perturbations"])[keep]
                for p, included in zip(manifest["perturbations"], keep):
                    exclusion_rows.append(dict(target=target, representation=representation,
                                               perturbation=p, included=int(included)))
                for name in MODELS[1:]:
                    q = represent(predictions[name], representation)[keep]
                    summary, rows = retrieval(q, actual[keep])
                    prefix = dict(target=target, representation=representation, model=name, repeat=0, direction="full")
                    summary_rows.append(dict(**prefix, **summary))
                    for i, p in enumerate(ids):
                        rank_rows.append(dict(**prefix, perturbation=p, **{k: v[i] for k,v in rows.items()}))
                for r in range(ds.shape[0]):
                    a = represent(ds[r, 0].astype(np.float64), representation)[keep]
                    b = represent(ds[r, 1].astype(np.float64), representation)[keep]
                    for direction, q, c in (("half0_to_half1", a, b), ("half1_to_half0", b, a)):
                        summary, _ = retrieval(q, c)
                        summary_rows.append(dict(target=target, representation=representation,
                                                 model="experimental_half_reference", repeat=r+1,
                                                 direction=direction, **summary))
        paired_mse = errors["paired_calibrated"].mean()
        common_mse = errors["common_calibrated"].mean()
        print(f"MSE_OK {target}: calibrated paired={paired_mse:.8f}; "
              f"calibrated common={common_mse:.8f}; paired-minus-common={paired_mse-common_mse:+.8f}", flush=True)
        for representation in REPRESENTATIONS:
            a = next(row for row in summary_rows if row["target"] == target
                     and row["representation"] == representation and row["model"] == "paired_calibrated")
            print(f"RETRIEVAL_OK {target} {representation}: N={a['n_candidates']}; "
                  f"Top1={a['top1']:.4%}; Top10={a['top10']:.4%}; "
                  f"MRR={a['mrr']:.6f}; chanceMRR={a['chance_mrr']:.6f}", flush=True)
    for name, rows in (("model_comparison", model_rows), ("task_errors", task_rows),
                       ("retrieval_summary", summary_rows), ("retrieval_tasks", rank_rows),
                       ("retrieval_cohort", exclusion_rows), ("half_batch_composition", batch_rows)):
        bb.write_table(out / (name + ".tsv"), rows)
    vi.json_write(out / "complete.json", {"protocol_sha256": bb.digest(out / "protocol.json"),
                                         "outputs": {p.name: bb.digest(p) for p in out.glob("*.tsv")}})
    print(f"ALL_OK {out}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.home() / "VC_reliability")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = args.output or args.root / "results/specificity_v01"
    manifest = json.loads((args.root / "processed/v01/manifest.json").read_text())
    freeze(args.root, out, manifest)
    evaluate(args.root, out, manifest)


if __name__ == "__main__":
    main()
