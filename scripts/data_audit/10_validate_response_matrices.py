#!/usr/bin/env python3
"""Check the arithmetic of stored GEM-matched perturbation responses.

This reads precomputed summaries only. It never trains or selects a model and
does not read target responses for use as predictors.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def parse_input(value):
    env, separator, path = value.partition("=")
    if not separator or not env or not path:
        raise argparse.ArgumentTypeError("Expected ENVIRONMENT=PATH")
    return env, Path(path)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def inspect(path, env, feature_ids, task_names, min_controls, chunk_size, tolerance):
    rows = []
    with h5py.File(path, "r") as handle:
        if handle["gene_id"].asstr()[...].tolist() != feature_ids:
            raise ValueError(f"{env}: expression gene IDs or order changed")
        if handle["perturbation"].asstr()[...].tolist() != task_names:
            raise ValueError(f"{env}: perturbation labels or order changed")
        gems = handle["gem_group"].asstr()[...].tolist()
        if len(gems) != len(set(gems)) or not gems:
            raise ValueError(f"{env}: invalid GEM labels")
        n = handle["n_cells"][...]
        by_gem = handle["n_cells_by_batch"][...]
        control_n = handle["control_cell_count_by_batch"][...]
        control_means = handle["control_mean_by_batch"][...]
        if (
            n.shape != (len(task_names),)
            or by_gem.shape != (len(task_names), len(gems))
            or control_n.shape != (len(gems),)
            or control_means.shape != (len(gems), len(feature_ids))
            or np.any(n <= 0)
            or np.any(by_gem < 0)
            or np.any(control_n <= 0)
            or not np.array_equal(by_gem.sum(axis=1), n)
            or not np.isfinite(control_means).all()
        ):
            raise ValueError(f"{env}: invalid stored cell counts or control means")
        for key in ("perturbed_mean_log1p", "matched_control_mean_log1p", "delta"):
            if handle[key].shape != (len(task_names), len(feature_ids)):
                raise ValueError(f"{env}: unexpected shape for {key}")
        if handle["low_control_fraction"].shape != (len(task_names),):
            raise ValueError(f"{env}: unexpected low_control_fraction shape")

        for start in range(0, len(task_names), chunk_size):
            end = min(start + chunk_size, len(task_names))
            weights = by_gem[start:end].astype(np.float64)
            matched_expected = weights @ control_means.astype(np.float64)
            matched_expected /= n[start:end, None]
            matched = handle["matched_control_mean_log1p"][start:end]
            perturbed = handle["perturbed_mean_log1p"][start:end]
            delta = handle["delta"][start:end]
            low = handle["low_control_fraction"][start:end]
            if not all(np.isfinite(array).all() for array in (matched, perturbed, delta, low)):
                raise ValueError(f"{env}: nonfinite stored response at rows {start}:{end}")

            matched_error = np.max(np.abs(matched_expected - matched), axis=1)
            delta_expected = perturbed.astype(np.float64) - matched.astype(np.float64)
            delta_error = np.max(np.abs(delta_expected - delta), axis=1)
            low_expected = weights[:, control_n < min_controls].sum(axis=1) / n[start:end]
            low_error = np.abs(low_expected - low)
            for offset, task in enumerate(task_names[start:end]):
                rows.append((
                    env, task, int(n[start + offset]),
                    float(matched_error[offset]), float(delta_error[offset]),
                    float(low_error[offset]), float(low_expected[offset]),
                ))

    max_matched = max(row[3] for row in rows)
    max_delta = max(row[4] for row in rows)
    max_low = max(row[5] for row in rows)
    summary = {
        "file_name": path.name,
        "n_tasks": len(rows),
        "n_features": len(feature_ids),
        "n_gems": len(gems),
        "min_controls_threshold": min_controls,
        "n_tasks_with_any_gem_below_threshold": sum(row[6] > 0 for row in rows),
        "max_abs_matched_control_difference": max_matched,
        "max_abs_delta_difference": max_delta,
        "max_abs_low_control_fraction_difference": max_low,
        "pass": (max_matched <= tolerance and max_delta <= tolerance and max_low <= 1e-6),
    }
    return summary, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=parse_input, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--tasks", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--tolerance", type=float, default=5e-6)
    args = parser.parse_args()
    if args.chunk_size < 1 or not 0 < args.tolerance < 1e-3:
        parser.error("Invalid chunk size or arithmetic tolerance")
    if len({env for env, _ in args.input}) != len(args.input):
        parser.error("Duplicate environment label")

    features = [row["gene_id"] for row in read_tsv(args.features)]
    tasks = [row["gene"] for row in read_tsv(args.tasks) if row["in_main_panel"] == "1"]
    if not features or not tasks or len(set(features)) != len(features) or len(set(tasks)) != len(tasks):
        parser.error("Empty or duplicated reference gene or task axis")
    manifest_bytes = args.manifest.read_bytes()
    manifest = json.loads(manifest_bytes)
    min_controls = int(manifest["min_controls"])
    if min_controls < 1:
        parser.error("manifest min_controls must be positive")
    provenance_report = json.loads(args.provenance.read_text(encoding="utf-8"))
    if provenance_report.get("manifest_file_sha256") != hashlib.sha256(manifest_bytes).hexdigest():
        parser.error("The manifest has changed since the provenance audit")
    provenance = provenance_report["environments"]
    if set(provenance) != {env for env, _ in args.input}:
        parser.error("Provenance audit and input environment labels differ")
    for env, _ in args.input:
        if not provenance[env]["checks"] or not all(provenance[env]["checks"].values()):
            parser.error(f"Provenance count/axis checks did not all pass for {env}")

    summaries = {}
    all_rows = []
    for env, path in args.input:
        summaries[env], rows = inspect(
            path, env, features, tasks, min_controls, args.chunk_size, args.tolerance
        )
        all_rows.extend(rows)
        print(f"MATRIX_CHECK {env} {json.dumps(summaries[env])}", flush=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "protocol": "precomputed response arithmetic audit; no model fitting",
        "arithmetic_tolerance": args.tolerance,
        "low_control_fraction_tolerance": 1e-6,
        "control_matching_rule": "sum_GEM (perturbation_cells_GEM / perturbation_cells_total) * control_mean_GEM",
        "delta_rule": "perturbed_mean_log1p - matched_control_mean_log1p",
        "environments": summaries,
    }
    (args.output_dir / "response_matrix_validation.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    with (args.output_dir / "response_matrix_validation_by_task.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "environment", "perturbation", "n_cells", "matched_control_max_abs_diff",
            "delta_max_abs_diff", "low_control_fraction_abs_diff", "low_control_fraction",
        ))
        writer.writerows(all_rows)
    if not all(record["pass"] for record in summaries.values()):
        raise SystemExit("Response matrix arithmetic differs from the declared rules; inspect QC output")


if __name__ == "__main__":
    main()
