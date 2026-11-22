#!/usr/bin/env python3
"""Compare guide-supported tasks in HCT116 and HEK293T by exact Ensembl ID."""

import argparse
import csv
import json
from importlib import import_module
from pathlib import Path

audit = import_module("05_xatlas_observed_support")
SUPPORT_COLUMNS = (
    "guide_validated_cells",
    "cells_in_gems_with_controls",
    "cells_in_gems_with_ge30_controls",
)


def read_support(run, cell_line):
    summary = json.loads((run / "summary.json").read_text(encoding="utf-8"))
    if summary["cell_line"] != cell_line or summary["status"] != "GUIDE_METADATA_ONLY":
        raise ValueError(f"Unexpected audit cohort in {run}")
    tasks = {}
    with (run / "task_support.tsv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            gene_id = row["gene_id"]
            if gene_id in tasks:
                raise ValueError(f"Duplicate target ID {gene_id} in {run}")
            tasks[gene_id] = row
    if len(tasks) != summary["main_panel_tasks"]:
        raise ValueError(f"Task table and summary disagree in {run}")
    return tasks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hct116_dir", type=Path)
    parser.add_argument("hek293t_dir", type=Path)
    args = parser.parse_args()

    hct116 = read_support(args.hct116_dir, "HCT116")
    hek293t = read_support(args.hek293t_dir, "HEK293T")
    hct_symbols = {gene_id: row["panel_symbol"] for gene_id, row in hct116.items()}
    hek_symbols = {gene_id: row["panel_symbol"] for gene_id, row in hek293t.items()}
    if hct_symbols != hek_symbols:
        raise ValueError("HCT116 and HEK293T must use the same main panel")

    minimum = audit.MIN_TASK_CELLS
    summary = {
        "status": "GUIDE_METADATA_ONLY",
        "main_panel_tasks": len(hct116),
        "diagnostic_min_cells_per_task": minimum,
        "diagnostic_min_controls_per_gem": audit.CONTROL_DIAGNOSTIC,
        "tasks_meeting_cell_cutoff": {},
    }
    for column in SUPPORT_COLUMNS:
        hct_ids = {g for g, row in hct116.items() if int(row[column]) >= minimum}
        hek_ids = {g for g, row in hek293t.items() if int(row[column]) >= minimum}
        summary["tasks_meeting_cell_cutoff"][column] = {
            "HCT116": len(hct_ids), "HEK293T": len(hek_ids), "both": len(hct_ids & hek_ids),
        }

    out = args.hek293t_dir
    fields = ["gene_id", "panel_symbol"]
    rows = []
    for cell_line in ("HCT116", "HEK293T"):
        fields.extend(f"{cell_line}_{column}" for column in SUPPORT_COLUMNS)
    for gene_id, symbol in sorted(hct_symbols.items()):
        counts = [run[gene_id][column] for run in (hct116, hek293t) for column in SUPPORT_COLUMNS]
        rows.append([gene_id, symbol, *counts])
    audit.design.write_table(out / "paired_task_support.tsv", fields, rows)
    audit.save_json(out / "paired_support_summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
