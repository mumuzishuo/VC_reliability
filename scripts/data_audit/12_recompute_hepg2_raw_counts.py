#!/usr/bin/env python3
"""Audit selected HepG2 response rows directly against raw single-cell counts.

The prepare phase selects labels using only cell and GEM/control metadata. It
does not open raw X or read any target response. The audit phase independently
normalizes float64 counts, aggregates perturbations and GEM controls, and only
then reads the stored response arrays for comparison. This is an arithmetic
audit, not a prospective assessment of prediction reliability.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import h5py
import numpy as np

ENVIRONMENT = "HepG2"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def strings(dataset):
    return np.asarray(dataset.asstr()[:], dtype=str)


def write_json(path, value):
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def load_state(args):
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    require(manifest["target"] == ENVIRONMENT, "Unexpected target environment")
    require(manifest["script_sha256"] == digest(args.preprocess_code),
            "Preprocessing code differs from manifest")
    item = manifest["inputs"][ENVIRONMENT]
    raw_path = Path(item["path"])
    stat = raw_path.stat()
    require(stat.st_size == item["size_bytes"] and stat.st_mtime_ns == item["mtime_ns"],
            "Raw input size or modification time differs from manifest")
    with h5py.File(args.response, "r") as h:
        require(h.attrs["role"] == "evaluation_only" and
                h.attrs["environment"] == ENVIRONMENT and
                h.attrs["manifest_sha256"] == digest(args.manifest),
                "Target response provenance differs from manifest")
        labels, features, gems = (strings(h[key]) for key in
                                  ("perturbation", "gene_id", "gem_group"))
        require(np.array_equal(labels, manifest["perturbations"]),
                "Perturbation order differs from manifest")
        require(np.array_equal(features, manifest["gene_ids"]),
                "Feature order differs from manifest")
        controls = h["control_cell_count_by_batch"][:]
        by_gem = h["n_cells_by_batch"][:]
        n_cells = h["n_cells"][:]
        require(by_gem.shape == (len(labels), len(gems)) and
                np.array_equal(by_gem.sum(axis=1), n_cells),
                "Saved task cell/GEM counts are inconsistent")
        require(controls.shape == (len(gems),) and
                (controls >= manifest["min_controls"]).all(),
                "Saved control counts are inconsistent")
        require(h["delta"].shape == (len(labels), len(features)),
                "Target response dimensions differ")
    return manifest, raw_path, labels, features, gems, controls, by_gem, n_cells


def choose_labels(labels, gems, controls, by_gem, n_cells, cutoff):
    chosen = {}

    def add(index, reason):
        chosen.setdefault(str(labels[index]), []).append(reason)

    order = sorted(range(len(labels)), key=lambda i: (int(n_cells[i]), str(labels[i])))
    for q in range(5):
        add(order[((len(order) - 1) * q) // 4], "cell_support_quintile_" + str(q))
    frac = by_gem.max(axis=1) / n_cells
    add(sorted(range(len(labels)), key=lambda i: (-frac[i], str(labels[i])))[0],
        "max_single_gem_fraction")
    low_gems = [j for j in range(len(gems)) if controls[j] < cutoff]
    for j in low_gems:
        eligible = [i for i in range(len(labels)) if by_gem[i, j] > 0]
        require(eligible, "No selected-panel perturbations in low-control GEM " + gems[j])
        ranked = sorted(eligible, key=lambda i: (-by_gem[i, j] / n_cells[i],
                                                  -by_gem[i, j], str(labels[i])))
        for rank, index in enumerate(ranked[:2], start=1):
            add(index, "low_control_gem_{}_rank_{}".format(gems[j], rank))
    return chosen, low_gems


def prepare(args, state):
    require(not args.output_dir.exists(), "Output directory already exists; use a new version")
    manifest, raw_path, labels, features, gems, controls, by_gem, n_cells = state
    chosen, low_gems = choose_labels(labels, gems, controls, by_gem, n_cells,
                                    args.low_control_cutoff)
    args.output_dir.mkdir(parents=True)
    selection = args.output_dir / "selection.tsv"
    indexes = {p: i for i, p in enumerate(labels)}
    with selection.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("perturbation", "n_cells", "gem_groups",
                         "low_control_cell_fraction", "selection_reasons"))
        for name in sorted(chosen):
            i = indexes[name]
            writer.writerow((name, int(n_cells[i]), int((by_gem[i] > 0).sum()),
                             float(by_gem[i, low_gems].sum() / n_cells[i]),
                             "|".join(sorted(chosen[name]))))
    write_json(args.output_dir / "selection.json", {
        "purpose": "arithmetic QC; never use selection or target results to tune risk rules",
        "rule": "5 deterministic cell-support quantiles; maximum single-GEM fraction; "
                "top two cell-fractions for each GEM with fewer than cutoff controls",
        "low_control_cutoff": args.low_control_cutoff,
        "low_control_gems": {str(gems[j]): int(controls[j]) for j in low_gems},
        "n_selected": len(chosen),
        "selection_tsv_sha256": digest(selection),
        "manifest_sha256": digest(args.manifest),
        "response_sha256": digest(args.response),
        "preprocess_code_sha256": digest(args.preprocess_code),
        "raw_file_name": raw_path.name,
        "raw_size_bytes": raw_path.stat().st_size,
        "raw_mtime_ns": raw_path.stat().st_mtime_ns,
        "n_features": len(features),
        "n_tasks": len(labels),
    })
    print("SELECTION_FROZEN", json.dumps({"n_selected": len(chosen),
                                            "low_control_gems": [gems[j] for j in low_gems]}),
          flush=True)


def metadata_reader(path):
    spec = importlib.util.spec_from_file_location("vc_preprocess_audit", str(path))
    require(spec is not None and spec.loader is not None, "Cannot import preprocessing metadata reader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.metadata


def audit(args, state):
    manifest, raw_path, labels, features, gems, controls, by_gem, n_cells = state
    selection = args.output_dir / "selection.tsv"
    meta = json.loads((args.output_dir / "selection.json").read_text(encoding="utf-8"))
    require(meta["selection_tsv_sha256"] == digest(selection) and
            meta["manifest_sha256"] == digest(args.manifest) and
            meta["response_sha256"] == digest(args.response) and
            meta["preprocess_code_sha256"] == digest(args.preprocess_code) and
            meta["low_control_cutoff"] == args.low_control_cutoff,
            "Frozen selection, response, code or manifest has changed")
    require(raw_path.name == meta["raw_file_name"] and
            raw_path.stat().st_size == meta["raw_size_bytes"] and
            raw_path.stat().st_mtime_ns == meta["raw_mtime_ns"],
            "Raw input differs from frozen selection")
    expected, _ = choose_labels(labels, gems, controls, by_gem, n_cells,
                                args.low_control_cutoff)
    with selection.open(encoding="utf-8", newline="") as handle:
        selected = list(csv.DictReader(handle, delimiter="\t"))
    require([r["perturbation"] for r in selected] == sorted(expected) and
            all(r["selection_reasons"] == "|".join(sorted(expected[r["perturbation"]]))
                for r in selected), "Saved selection differs from deterministic metadata rule")
    require(not (args.output_dir / "per_task.tsv").exists() and
            not (args.output_dir / "audit_summary.json").exists(),
            "Audit results already exist; refusing overwrite")

    # The legacy reader decodes AnnData categories, while normalization and
    # GEM-weighted aggregation below are independent of the legacy estimator.
    cell_labels, cell_gems, raw_ids, umi = metadata_reader(args.preprocess_code)(raw_path)
    require(np.array_equal(np.sort(np.intersect1d(raw_ids, features)), np.sort(features)),
            "Raw expression feature IDs do not match target response")
    positions = {g: i for i, g in enumerate(raw_ids)}
    cols = np.asarray([positions[g] for g in features], dtype=np.int64)
    require(len(cell_labels) == manifest["inputs"][ENVIRONMENT]["n_cells"],
            "Raw cell count differs from manifest")
    target_index = {str(p): i for i, p in enumerate(labels)}
    chosen_names = [r["perturbation"] for r in selected]
    chosen_index = [target_index[p] for p in chosen_names]
    needed = by_gem[chosen_index].sum(axis=0) > 0
    gem_index = {str(g): j for j, g in enumerate(gems)}
    control = cell_labels == manifest["control_label"]
    valid_gem = np.isin(cell_gems, gems)
    raw_control_names, raw_control_counts = np.unique(cell_gems[control], return_counts=True)
    raw_controls = dict(zip(raw_control_names, raw_control_counts))
    require(set(gems) == {str(g) for g, count in raw_controls.items()
                          if count >= manifest["min_controls"]},
            "Raw eligible GEMs differ from target response")
    require(all(raw_controls.get(str(g), 0) == int(n) for g, n in zip(gems, controls)),
            "Raw GEM control counts differ from target response")
    selected_mask = np.isin(cell_labels, chosen_names) & valid_gem
    control_mask = control & np.isin(cell_gems, gems[needed])
    rows_needed = selected_mask | control_mask
    p_lookup = {name: i for i, name in enumerate(chosen_names)}
    p_sums = np.zeros((len(chosen_names), len(features)), dtype=np.float64)
    c_sums = np.zeros((len(gems), len(features)), dtype=np.float64)
    observed_p = np.zeros((len(chosen_names), len(gems)), dtype=np.int64)
    observed_c = np.zeros(len(gems), dtype=np.int64)
    ratio_min, ratio_max = np.inf, -np.inf
    n_raw_rows = 0

    # Target response values have not been read at this point.
    with h5py.File(raw_path, "r") as h:
        x = h["X"]
        require(isinstance(x, h5py.Dataset) and
                x.shape == (len(cell_labels), len(raw_ids)), "Expected dense raw X")
        for start in range(0, len(cell_labels), args.chunk_size):
            end = min(start + args.chunk_size, len(cell_labels))
            keep = rows_needed[start:end]
            if not keep.any():
                continue
            raw = x[start:end, :][keep]
            totals = umi[start:end][keep]
            ratio = raw.sum(axis=1, dtype=np.float64) / totals
            require(np.isfinite(raw).all() and (raw >= 0).all() and
                    np.equal(raw, np.floor(raw)).all() and
                    np.isfinite(ratio).all() and (ratio <= 1.0001).all(),
                    "Invalid raw count or UMI denominator")
            ratio_min = min(ratio_min, float(ratio.min()))
            ratio_max = max(ratio_max, float(ratio.max()))
            y = np.log1p(10000.0 * raw[:, cols].astype(np.float64) / totals[:, None])
            b = cell_gems[start:end][keep]
            p = cell_labels[start:end][keep]
            for batch in np.unique(b):
                j = gem_index[str(batch)]
                in_batch = b == batch
                is_control = in_batch & (p == manifest["control_label"])
                if is_control.any():
                    c_sums[j] += y[is_control].sum(axis=0, dtype=np.float64)
                    observed_c[j] += int(is_control.sum())
                for name in np.unique(p[in_batch & ~is_control]):
                    i = p_lookup[str(name)]
                    in_task = in_batch & (p == name)
                    p_sums[i] += y[in_task].sum(axis=0, dtype=np.float64)
                    observed_p[i, j] += int(in_task.sum())
            n_raw_rows += len(raw)
            if end == len(cell_labels) or (start // args.chunk_size) % 100 == 0:
                print("HepG2 raw cells scanned: {}/{}".format(end, len(cell_labels)), flush=True)

    require(np.array_equal(observed_p, by_gem[chosen_index]),
            "Raw task GEM counts differ from target response")
    require(np.array_equal(observed_c[needed], controls[needed]),
            "Raw GEM controls used differ from target response")
    c_means = np.zeros_like(c_sums)
    c_means[needed] = c_sums[needed] / observed_c[needed, None]
    p_means = p_sums / n_cells[chosen_index, None]
    matched = (observed_p / n_cells[chosen_index, None]) @ c_means
    delta = p_means - matched

    rows = []
    with h5py.File(args.response, "r") as h:
        control_difference = max(float(np.max(np.abs(c_means[j] -
                          h["control_mean_by_batch"][j, :]))) for j in np.flatnonzero(needed))
        for i, target_row in enumerate(chosen_index):
            diffs = {}
            for key, computed, field in (
                ("perturbed_mean", p_means[i], "perturbed_mean_log1p"),
                ("matched_control", matched[i], "matched_control_mean_log1p"),
                ("delta", delta[i], "delta"),
            ):
                reference = h[field][target_row, :]
                diffs[key + "_max_abs_diff"] = float(np.max(np.abs(computed - reference)))
            low = int(observed_p[i, controls < args.low_control_cutoff].sum())
            low_fraction_difference = abs(low / n_cells[target_row] -
                                          float(h["low_control_fraction"][target_row]))
            rows.append({"perturbation": chosen_names[i],
                         "n_cells": int(n_cells[target_row]),
                         "gem_groups": int(np.count_nonzero(observed_p[i])),
                         "cells_in_low_control_gems": low,
                         "low_control_fraction_max_abs_diff": low_fraction_difference,
                         **diffs})

    summary = {"environment": ENVIRONMENT, "status": "PASS", "n_tasks": len(rows),
               "n_features": len(features), "n_gems_used": int(needed.sum()),
               "n_raw_rows_checked": n_raw_rows, "max_control_mean_abs_diff": control_difference,
               "max_perturbed_mean_abs_diff": max(r["perturbed_mean_max_abs_diff"] for r in rows),
               "max_matched_control_abs_diff": max(r["matched_control_max_abs_diff"] for r in rows),
               "max_delta_abs_diff": max(r["delta_max_abs_diff"] for r in rows),
               "max_low_control_fraction_abs_diff": max(r["low_control_fraction_max_abs_diff"] for r in rows),
               "raw_count_over_umi_min": ratio_min, "raw_count_over_umi_max": ratio_max,
               "absolute_tolerance": args.absolute_tolerance,
               "selection_json_sha256": digest(args.output_dir / "selection.json"),
               "selection_tsv_sha256": digest(selection),
               "numpy_version": np.__version__, "h5py_version": h5py.__version__,
               "python_version": sys.version.split()[0]}
    require(np.isfinite(list(summary[key] for key in ("max_control_mean_abs_diff",
                      "max_perturbed_mean_abs_diff", "max_matched_control_abs_diff",
                      "max_delta_abs_diff", "max_low_control_fraction_abs_diff"))).all(),
            "Nonfinite audit difference")
    summary["status"] = "PASS" if all(summary[key] <= args.absolute_tolerance for key in (
        "max_control_mean_abs_diff", "max_perturbed_mean_abs_diff",
        "max_matched_control_abs_diff", "max_delta_abs_diff",
        "max_low_control_fraction_abs_diff")) else "FAIL"
    with (args.output_dir / "per_task.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output_dir / "audit_summary.json", summary)
    print("HEPG2_RAW_AUDIT", json.dumps(summary, sort_keys=True), flush=True)
    require(summary["status"] == "PASS", "Raw and saved target responses differ beyond tolerance")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "audit"))
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--response", required=True, type=Path)
    parser.add_argument("--preprocess-code", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--low-control-cutoff", type=int, default=30)
    parser.add_argument("--absolute-tolerance", type=float, default=2e-6)
    parser.add_argument("--chunk-size", type=int, default=512)
    args = parser.parse_args()
    require(args.low_control_cutoff > 0 and args.absolute_tolerance > 0 and
            args.chunk_size > 0, "Cutoff, tolerance and chunk size must be positive")
    state = load_state(args)
    if args.mode == "prepare":
        prepare(args, state)
    else:
        audit(args, state)


if __name__ == "__main__":
    main()
