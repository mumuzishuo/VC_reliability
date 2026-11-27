#!/usr/bin/env python3
"""Summarize X-Atlas task support across GEM batches and guide pairs."""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

MIN_CELLS = 30
MIN_CONTROLS = 30


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_cohort(directory, expected_name):
    names = ("summary.json", "gem_control_counts.tsv", "task_support.tsv",
             "task_gem_counts.tsv")
    paths = {name: directory / name for name in names}
    summary = json.loads(paths["summary.json"].read_text(encoding="utf-8"))
    if summary["status"] != "GUIDE_METADATA_ONLY" or summary["cell_line"] != expected_name:
        raise ValueError("Unexpected cohort summary in " + str(directory))

    controls = {}
    for row in read_tsv(paths["gem_control_counts.tsv"]):
        gem = row["gem"]
        if gem in controls:
            raise ValueError("Duplicate GEM: " + gem)
        controls[gem] = int(row["non_targeting_cells"])
    if len(controls) != summary["gems_in_guide_csv"]:
        raise ValueError("GEM count differs from summary")
    if sum(controls.values()) != summary["non_targeting_cells"]:
        raise ValueError("Control count differs from summary")

    tasks = {}
    for row in read_tsv(paths["task_support.tsv"]):
        gene_id = row["gene_id"]
        if gene_id in tasks:
            raise ValueError("Duplicate task: " + gene_id)
        tasks[gene_id] = row
    if len(tasks) != summary["main_panel_tasks"]:
        raise ValueError("Task count differs from summary")

    by_task = defaultdict(dict)
    for row in read_tsv(paths["task_gem_counts.tsv"]):
        gene_id, gem = row["gene_id"], row["gem"]
        if gene_id not in tasks or gem not in controls or gem in by_task[gene_id]:
            raise ValueError("Unknown or repeated task/GEM combination")
        if int(row["non_targeting_cells_in_gem"]) != controls[gem]:
            raise ValueError("Task/GEM control count differs from GEM table")
        n = int(row["guide_validated_cells"])
        if n <= 0:
            raise ValueError("Task/GEM count must be positive")
        by_task[gene_id][gem] = n

    rows = []
    eligible = set()
    multi_gem = set()
    multi_pair = set()
    for gene_id, task in sorted(tasks.items()):
        counts = by_task[gene_id]
        n = int(task["guide_validated_cells"])
        pairs = int(task["observed_guide_pairs"])
        low = sum(value for gem, value in counts.items() if controls[gem] < MIN_CONTROLS)
        if n != sum(counts.values()) or len(counts) != int(task["observed_gems"]):
            raise ValueError("Task cell/GEM count differs from task table: " + gene_id)
        if n - low != int(task["cells_in_gems_with_ge30_controls"]):
            raise ValueError("Well-controlled task count differs: " + gene_id)
        if n >= MIN_CELLS:
            eligible.add(gene_id)
            if len(counts) >= 2:
                multi_gem.add(gene_id)
            if pairs >= 2:
                multi_pair.add(gene_id)
        largest = max(counts.values()) / n if n else None
        effective = n * n / sum(value * value for value in counts.values()) if n else None
        rows.append({
            "cell_line": expected_name,
            "gene_id": gene_id,
            "panel_symbol": task["panel_symbol"],
            "cells": n,
            "observed_guide_pairs": pairs,
            "observed_gems": len(counts),
            "largest_gem_cell_fraction": largest,
            "effective_gem_count": effective,
            "cells_in_gems_lt30_controls": low,
            "minimum_controls_in_used_gem": min((controls[g] for g in counts), default=None),
        })
    if len(eligible) != summary["main_panel_targets_ge30_cells"]:
        raise ValueError("Eligible task count differs from summary")
    return rows, eligible, multi_gem, multi_pair, {str(p): digest(p) for p in paths.values()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hct116-dir", type=Path, required=True)
    parser.add_argument("--hek293t-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    h = read_cohort(args.hct116_dir, "HCT116")
    k = read_cohort(args.hek293t_dir, "HEK293T")
    h_panel = {r["gene_id"]: r["panel_symbol"] for r in h[0]}
    k_panel = {r["gene_id"]: r["panel_symbol"] for r in k[0]}
    if h_panel != k_panel:
        raise ValueError("Cohorts have different task panels")
    results = {"status": "GUIDE_METADATA_ONLY", "min_cells": MIN_CELLS,
               "min_controls_per_gem": MIN_CONTROLS, "cohorts": {}, "both": {}}
    for name, data in (("HCT116", h), ("HEK293T", k)):
        results["cohorts"][name] = {
            "tasks_ge30_cells": len(data[1]),
            "tasks_ge30_cells_in_at_least_two_gems": len(data[2]),
            "tasks_ge30_cells_with_at_least_two_guide_pairs": len(data[3]),
        }
    results["both"] = {
        "tasks_ge30_cells": len(h[1] & k[1]),
        "tasks_ge30_cells_in_at_least_two_gems": len(h[2] & k[2]),
        "tasks_ge30_cells_with_at_least_two_guide_pairs": len(h[3] & k[3]),
    }

    args.output_dir.mkdir(parents=True, exist_ok=False)
    columns = list(h[0][0])
    with (args.output_dir / "task_sample_support.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(h[0] + k[0])
    (args.output_dir / "summary.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    parameters = {"inputs_sha256": {**h[4], **k[4]},
                  "interpretation": "sample is a GEM capture batch, not a biological replicate"}
    (args.output_dir / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
