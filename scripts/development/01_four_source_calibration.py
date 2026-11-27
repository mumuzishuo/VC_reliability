#!/usr/bin/env python3
"""Build four-source predictions and assess calibration in four internal folds."""

import argparse
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "legacy_v01"))
import vc_baseline as baseline
import vc_diagnostics as diagnostics
import vc_internal as internal
import vc_preprocess as preprocess


SOURCE_ROLES = {"K562": "sources", "RPE1": "sources", "Jurkat": "sources",
                "HepG2": "evaluation_only"}
SEED = 20261127
COVERAGE = 0.8
NORMALIZATION = "ln(1 + 10000 * X / obs.UMI_count)"
require = baseline.require


def read_sources(root, checksum_path, provenance_path):
    processed = root / "processed/v01"
    manifest_path = processed / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_hash = baseline.digest(manifest_path)
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    require(provenance["manifest_file_sha256"] == manifest_hash,
            "The preprocessing manifest differs from the response provenance audit")
    require(manifest["script_sha256"] == baseline.digest(preprocess.__file__),
            "Archived preprocessing code does not match the preprocessing manifest")
    require(manifest["normalization"] == NORMALIZATION
            and manifest["min_cells"] == 30 and manifest["min_controls"] == 1,
            "Unexpected v01 preprocessing settings")
    tasks, genes = manifest["perturbations"], manifest["gene_ids"]
    require(bool(tasks) and bool(genes) and len(set(tasks)) == len(tasks)
            and len(set(genes)) == len(genes), "Empty or duplicated manifest axes")

    expected = {}
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        digest, filename = line.split(maxsplit=1)
        require(filename not in expected, "Duplicate checksum entry: " + filename)
        expected[filename] = digest

    deltas = np.empty((len(SOURCE_ROLES), len(tasks), len(genes)), dtype=np.float64)
    counts = np.empty((len(SOURCE_ROLES), len(tasks)), dtype=np.int64)
    inputs = {"processed/v01/manifest.json": manifest_hash}
    for i, (environment, role) in enumerate(SOURCE_ROLES.items()):
        path = processed / role / f"{environment}_responses.h5"
        key = path.relative_to(root).as_posix()
        actual = baseline.digest(path)
        require(expected.get(key) == actual, "Response file differs from audit: " + key)
        p, g, response, n, low, info = baseline.read_response(path, role, environment)
        require(info == {"manifest_sha256": manifest_hash, "normalization": NORMALIZATION},
                "Response metadata mismatch: " + key)
        deltas[i], counts[i], _ = baseline.align(tasks, genes, p, g, response, n, low)
        require((counts[i] >= manifest["min_cells"]).all(), "Low cell support: " + key)
        inputs[key] = actual
    return tasks, genes, deltas, counts, inputs


def fit_scales(deltas):
    """Fit a separate constrained scalar for each predictor across source folds."""
    n_sources, n_tasks, _ = deltas.shape
    require(n_sources >= 2, "Calibration needs at least two source contexts")
    total = deltas.sum(axis=0)
    centroids = deltas.mean(axis=1)
    centroid_total = centroids.sum(axis=0)
    paired_cross = paired_power = common_cross = common_power = 0.0
    for truth, centroid in zip(deltas, centroids):
        paired = (total - truth) / (n_sources - 1)
        common = (centroid_total - centroid) / (n_sources - 1)
        paired_cross += np.einsum("ij,ij->", paired, truth)
        paired_power += np.einsum("ij,ij->", paired, paired)
        common_cross += n_tasks * np.dot(common, centroid)
        common_power += n_tasks * np.dot(common, common)

    result = {}
    for name, cross, power in (("paired", paired_cross, paired_power),
                               ("common", common_cross, common_power)):
        unconstrained = float(cross / power) if power > 0 else None
        scale = float(np.clip(unconstrained, 0, 1)) if unconstrained is not None else 0.0
        result[name] = {"scale": scale, "cross_product": float(cross),
                        "predictor_sum_squares": float(power)}
    return result


def predict(deltas, counts):
    require(len(deltas) >= 3, "The source transfer risk score needs at least three contexts")
    paired = deltas.mean(axis=0)
    common = paired.mean(axis=0)
    scales = fit_scales(deltas)
    excess = np.zeros(len(paired))
    for held_out, truth in enumerate(deltas):
        sources = deltas[np.arange(len(deltas)) != held_out]
        inner_scale = fit_scales(sources)["paired"]["scale"]
        inner_prediction = inner_scale * sources.mean(axis=0)
        excess += diagnostics.ms(inner_prediction - truth) - diagnostics.ms(truth)
    scores = {
        "source_heterogeneity": deltas.var(axis=0).mean(axis=1),
        "low_source_support": 1.0 / counts.min(axis=0),
        "predicted_rms": np.sqrt(diagnostics.ms(paired)),
        "source_held_out_excess_mse": excess / len(deltas),
    }
    return paired, common, scales, scores


def risk_rows(environment, scores, tie_key, error, zero):
    excess = error - zero
    failure = excess > 0
    k = math.ceil(COVERAGE * len(error))
    coverage = k / len(error)
    for name, score in scores.items():
        retained = baseline.risk_order(score, tie_key)[:k]
        retained_mse = float(error[retained].mean())
        policy_mse = float(zero.mean() + excess[retained].sum() / len(error))
        auroc, ap = diagnostics.classifier_metrics(score, failure)
        yield dict(
            held_out=environment, score=name, n_tasks=len(error), n_retained=k,
            actual_coverage=coverage, failure_prevalence=float(failure.mean()),
            auroc_D_gt_0=auroc, auprc_average_precision=ap,
            spearman_score_vs_mse=baseline.spearman(score, error),
            spearman_score_vs_D=baseline.spearman(score, excess),
            retained_mse=retained_mse, retained_zero_mse=float(zero[retained].mean()),
            retained_failure_rate=float(failure[retained].mean()),
            retained_mse_gain_vs_all=baseline.ratio_gain(retained_mse, float(error.mean())),
            policy_mse=policy_mse, policy_minus_always_predict=policy_mse - float(error.mean()),
            policy_minus_always_zero=policy_mse - float(zero.mean()),
            random_same_coverage_policy_mse=float(zero.mean() + coverage * excess.mean()),
        )


def internal_folds(tasks, deltas, counts):
    models, losses, risks = [], [], []
    for fold, environment in enumerate(SOURCE_ROLES):
        keep = np.arange(len(deltas)) != fold
        paired, common, scales, scores = predict(deltas[keep], counts[keep])
        alpha, beta = scales["paired"]["scale"], scales["common"]["scale"]
        truth = deltas[fold]
        zero = diagnostics.ms(truth)
        errors = {
            "zero_response": zero,
            "paired_raw": diagnostics.ms(paired - truth),
            "common_raw": diagnostics.ms(common - truth),
            "paired_calibrated": diagnostics.ms(alpha * paired - truth),
            "common_calibrated": diagnostics.ms(beta * common - truth),
        }
        for model, error in errors.items():
            models.append(dict(held_out=environment, model=model, alpha_paired=alpha,
                               alpha_common=beta, **diagnostics.summarize(error, zero),
                               median_mse=float(np.median(error)),
                               gain_vs_zero=baseline.ratio_gain(float(error.mean()), float(zero.mean()))))
            for i, task in enumerate(tasks):
                losses.append(dict(held_out=environment, perturbation=task, model=model,
                                   mse=float(error[i]), zero_mse=float(zero[i]),
                                   D=float(error[i] - zero[i]),
                                   improvement_vs_zero=baseline.ratio_gain(error[i], zero[i])))
        tie_key = np.random.default_rng(SEED + fold + 1).permutation(len(tasks))
        risks.extend(risk_rows(environment, scores, tie_key, errors["paired_calibrated"], zero))
    return models, losses, risks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--response-checksums", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir
    require(not out.exists(), "Use a new output directory: " + str(out))
    tasks, genes, deltas, counts, inputs = read_sources(
        args.legacy_root, args.response_checksums, args.provenance,
    )
    paired, common, scales, scores = predict(deltas, counts)
    models, losses, risks = internal_folds(tasks, deltas, counts)

    out.mkdir(parents=True, exist_ok=False)
    preprocess.write_h5(
        out / "predictions.h5",
        {"perturbation": tasks, "gene_id": genes, "paired_raw": paired,
         "common_raw": common, "source_n_cells": counts,
         "tie_key": np.random.default_rng(SEED).permutation(len(tasks)), **scores},
        {"status": "INTERNAL_DEVELOPMENT_ONLY", "sources": ",".join(SOURCE_ROLES),
         "alpha_paired": scales["paired"]["scale"], "alpha_common": scales["common"]["scale"],
         "normalization": NORMALIZATION},
    )
    baseline.write_table(out / "internal_model_metrics.tsv", models)
    baseline.write_table(out / "internal_task_errors.tsv", losses)
    baseline.write_table(out / "internal_risk_metrics.tsv", risks)
    internal.json_write(out / "parameters.json", {
        "source_file_roles": SOURCE_ROLES,
        "input_sha256_relative_to_legacy_root": inputs,
        "input_audit_sha256": {str(path): baseline.digest(path) for path in
                               (args.response_checksums, args.provenance)},
        "code_sha256": {Path(path).name: baseline.digest(path) for path in
                        (__file__, baseline.__file__, diagnostics.__file__, internal.__file__,
                         preprocess.__file__)},
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__,
                     "scipy": scipy.__version__, "h5py": h5py.__version__},
        "seed": SEED, "nominal_coverage": COVERAGE, "coverage_rounding": "ceil(0.8 * n_tasks)",
        "scale_bounds": [0, 1], "zero_power_scale": 0,
        "risk_event": "paired_calibrated MSE minus zero response MSE > 0",
        "risk_policy": "retain lowest scores and predict zero on the other tasks",
        "risk_scores": list(scores), "higher_scores_mean_higher_proposed_risk": True,
    })
    summary = {
        "status": "INTERNAL_DEVELOPMENT_ONLY", "n_tasks": len(tasks), "n_genes": len(genes),
        "four_source_calibration": scales,
        "internal_calibrated_paired": [
            {**row, "gain_vs_zero": row["gain_vs_zero"] if math.isfinite(row["gain_vs_zero"]) else None}
            for row in models if row["model"] == "paired_calibrated"
        ],
    }
    internal.json_write(out / "summary.json", summary)
    print(json.dumps(summary, indent=2, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
