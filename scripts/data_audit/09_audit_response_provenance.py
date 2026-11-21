#!/usr/bin/env python3
"""Check existing response metadata and GEM counts against the audited inputs."""

import argparse
import hashlib
import json
import math
from importlib import import_module
from pathlib import Path

import h5py
import numpy as np


axis_audit = import_module("08_audit_existing_responses")
RESPONSE_ATTRIBUTES = (
    "environment",
    "role",
    "protocol",
    "normalization",
    "manifest_sha256",
    "stored_count_fraction_min",
    "stored_count_fraction_max",
    "warning",
)


def json_value(value):
    if isinstance(value, (np.ndarray, np.generic)):
        value = value.tolist()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, list):
        return [json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def load_count_references(gem_dir, environments, tasks):
    controls = {env: {} for env in environments}
    for row in axis_audit.read_tsv(gem_dir / "gem_control_summary_all_environments.tsv"):
        env, gem = row["environment"], row["gem_group"]
        if env not in controls:
            continue
        count = int(row["control_cells"])
        if gem in controls[env] or count < 0:
            raise ValueError(f"Invalid control count for {env}, GEM {gem}")
        controls[env][gem] = count

    task_index = {row["gene"]: i for i, row in enumerate(tasks)}
    gem_index = {
        env: {gem: i for i, gem in enumerate(groups)}
        for env, groups in controls.items()
    }
    perturbed = {
        env: np.zeros((len(tasks), len(controls[env])), dtype=np.int64)
        for env in environments
    }
    path = gem_dir / "perturbation_cell_counts_by_gem_all_environments.tsv"
    for row in axis_audit.read_tsv(path):
        env, gene = row["environment"], row["gene"]
        if env not in controls or gene not in task_index:
            continue
        gem = row["gem_group"]
        i, j = task_index[gene], gem_index[env][gem]
        count = int(row["perturbation_cells"])
        if count <= 0 or perturbed[env][i, j] != 0:
            raise ValueError(f"Invalid perturbation count for {env}, {gene}, GEM {gem}")
        perturbed[env][i, j] = count

    return controls, gem_index, perturbed


def inspect_response(path, file_sha, feature_ids, tasks, gem_controls, gem_index, counts, env):
    with h5py.File(path, "r") as h5:
        attrs = {
            key: json_value(h5.attrs[key])
            for key in RESPONSE_ATTRIBUTES
            if key in h5.attrs
        }
        gems = h5["gem_group"].asstr()[...].tolist()
        task_names = [row["gene"] for row in tasks]
        genes_match = h5["gene_id"].asstr()[...].tolist() == feature_ids
        tasks_match = h5["perturbation"].asstr()[...].tolist() == task_names
        gems_match = len(gems) == len(gem_controls) and set(gems) == set(gem_controls)
        control_n = h5["control_cell_count_by_batch"][...]
        perturbation_n = h5["n_cells_by_batch"][...]
        total_n = h5["n_cells"][...]
        historical_low = h5["low_control_fraction"][...]

    expected_total = np.array([int(row[env]) for row in tasks])
    checks = {
        "gene_order_matches": genes_match,
        "task_order_matches": tasks_match,
        "gem_ids_match": gems_match,
        "controls_match": None,
        "perturbations_by_gem_match": None,
        "perturbation_totals_match": (
            bool(np.array_equal(total_n, expected_total)) if tasks_match else None
        ),
        "manifest_hash_matches_file_bytes": attrs.get("manifest_sha256") == file_sha,
    }
    low_diagnostics = {}
    if gems_match:
        expected_controls = np.array([gem_controls[gem] for gem in gems])
        checks["controls_match"] = bool(np.array_equal(control_n, expected_controls))
        if tasks_match:
            expected_counts = counts[:, [gem_index[gem] for gem in gems]]
            checks["perturbations_by_gem_match"] = bool(
                np.array_equal(perturbation_n, expected_counts)
            )
            if (
                historical_low.shape == expected_total.shape
                and np.all(np.isfinite(historical_low))
                and np.all(expected_total > 0)
            ):
                for threshold in (10, 20):
                    calculated = expected_counts[:, expected_controls < threshold].sum(axis=1)
                    calculated = calculated / expected_total
                    low_diagnostics[f"max_abs_difference_controls_lt{threshold}"] = float(
                        np.max(np.abs(calculated - historical_low))
                    )

    return {
        "file": path.name,
        "attributes": attrs,
        "checks": checks,
        "low_control_fraction_comparisons": low_diagnostics,
    }


def summarise_raw_audit(path):
    rows = axis_audit.read_tsv(path)
    columns = list(rows[0]) if rows else []
    differences = {}
    for name in columns:
        if "diff" not in name.lower() and "error" not in name.lower():
            continue
        values = []
        for row in rows:
            try:
                value = float(row[name])
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                values.append(abs(value))
        if values:
            differences[name] = {"finite_rows": len(values), "max_abs": max(values)}
    return {"rows": len(rows), "columns": columns, "difference_columns": differences}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, action="append", type=axis_audit.parse_input)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--tasks", required=True, type=Path)
    parser.add_argument("--gem-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--raw-audit", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    environments = [name for name, _ in args.input]
    if len(set(environments)) != len(environments):
        parser.error("Environment labels must be unique")
    features = [row["gene_id"] for row in axis_audit.read_tsv(args.features)]
    tasks = [
        row for row in axis_audit.read_tsv(args.tasks) if row["in_main_panel"] == "1"
    ]
    if not features or not tasks:
        parser.error("The reference feature and main task tables must be nonempty")

    manifest_bytes = args.manifest.read_bytes()
    file_sha = hashlib.sha256(manifest_bytes).hexdigest()
    manifest = json.loads(manifest_bytes)
    controls, gem_index, counts = load_count_references(args.gem_dir, environments, tasks)
    report = {
        "manifest_file_sha256": file_sha,
        "manifest_top_level_keys": sorted(manifest) if isinstance(manifest, dict) else None,
        "raw_audit": summarise_raw_audit(args.raw_audit),
        "environments": {},
    }
    for env, path in args.input:
        report["environments"][env] = inspect_response(
            path, file_sha, features, tasks, controls[env], gem_index[env], counts[env], env
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    for env, result in report["environments"].items():
        print(f"CHECK {env} {json.dumps(result['checks'])}")
        print(f"METADATA {env} {json.dumps(result['attributes'])}")
        if result["low_control_fraction_comparisons"]:
            print(f"LOW_CONTROL {env} {json.dumps(result['low_control_fraction_comparisons'])}")
    print(f"RAW_AUDIT {json.dumps(report['raw_audit'])}")


if __name__ == "__main__":
    main()
