#!/usr/bin/env python3
"""Summarize observed cells, guide pairs and controls from X-Atlas guide metadata."""

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from importlib import import_module
from pathlib import Path

# Reuse the panel reader and TSV writer from the guide-design audit.
design = import_module("02_xatlas_guide_overlap")
HEADER = (
    "cell_barcode", "num_features", "feature_call", "num_umis",
    "gene_target", "pass_guide_filter",
)
CONTROL = "Non-Targeting"
MIN_TASK_CELLS = 30
CONTROL_DIAGNOSTIC = 30


def pair_key(text):
    return tuple(sorted(text.strip().split("|")))


def read_library(path):
    library = {}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            pair = pair_key(row["id_ab"])
            gene_id = row["target_gene_id"].strip()
            symbol = row["target_gene_name"].strip()
            if len(pair) != 2 or not pair[0] or pair[0] == pair[1] or not symbol:
                raise ValueError(f"Invalid guide pair for {symbol!r}")
            if gene_id != CONTROL and not design.ENSEMBL.fullmatch(gene_id):
                raise ValueError(f"Invalid target ID {gene_id!r}")
            if pair in library and library[pair] != (gene_id, symbol):
                raise ValueError("Conflicting library annotations for the same guide pair")
            library[pair] = (gene_id, symbol)
    return library


def scan_calls(path, library, panel, cell_line):
    seen = set()
    excluded = Counter()
    cells = Counter()
    controls = Counter()
    control_pairs = defaultdict(set)
    targeting = Counter()
    observed_ids = set()
    task_gems = defaultdict(Counter)
    task_pairs = defaultdict(Counter)
    sample_prefix = cell_line + "_"

    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, strict=True)
        if tuple(next(reader, ())) != HEADER:
            raise ValueError("Guide call header differs from the verified six-column schema")
        for barcode, num_features, feature_call, _, gene_target, pass_filter in reader:
            if barcode in seen:
                raise ValueError(f"Repeated cell barcode at CSV line {reader.line_num}")
            seen.add(barcode)
            prefix, _, gem = barcode.partition("-")
            if not prefix or not gem.startswith(sample_prefix) or gem == sample_prefix:
                excluded["missing_or_unexpected_gem"] += 1
                continue
            cells[gem] += 1

            flag = pass_filter.strip().lower()
            if flag not in ("true", "1"):
                reason = "failed_guide_filter" if flag in ("false", "0") else "unknown_filter_value"
                excluded[reason] += 1
                continue
            if num_features.strip() != "2":
                excluded["pass_flag_without_two_guides"] += 1
                continue
            pair = pair_key(feature_call)
            if pair not in library:
                excluded["guide_pair_missing_in_library"] += 1
                continue
            gene_id, symbol = library[pair]
            label = gene_target.strip()
            if gene_id == CONTROL:
                label, symbol = label.casefold(), symbol.casefold()
            if label != symbol:
                excluded["gene_label_differs_from_library"] += 1
                continue

            if gene_id == CONTROL:
                controls[gem] += 1
                control_pairs[gem].add(pair)
            else:
                targeting[gem] += 1
                observed_ids.add(gene_id)
                if gene_id in panel:
                    task_gems[gene_id][gem] += 1
                    task_pairs[gene_id][pair] += 1

    if not cells:
        raise ValueError(f"No {cell_line} GEM suffixes found")
    total = len(seen)
    if total != sum(excluded.values()) + sum(controls.values()) + sum(targeting.values()):
        raise RuntimeError("Cell counts do not balance")
    return (total, excluded, cells, controls, control_pairs, targeting,
            observed_ids, task_gems, task_pairs)


def write_results(out, panel, library, counts, cell_line):
    (total, excluded, cells, controls, control_pairs, targeting,
     observed_ids, task_gems, task_pairs) = counts
    designed = Counter(gene_id for gene_id, _ in library.values() if gene_id in panel)
    panel_cells_by_gem = Counter()
    panel_counts = Counter()
    task_rows = []
    for gene_id, symbol in sorted(panel.items()):
        by_gem = task_gems[gene_id]
        panel_cells_by_gem.update(by_gem)
        n_cells = sum(by_gem.values())
        n_controlled = sum(n for gem, n in by_gem.items() if controls[gem] > 0)
        n_well_supported = sum(n for gem, n in by_gem.items() if controls[gem] >= CONTROL_DIAGNOSTIC)
        panel_counts["main_panel_targets_observed"] += n_cells > 0
        panel_counts["main_panel_targets_ge30_cells"] += n_cells >= MIN_TASK_CELLS
        panel_counts["main_panel_targets_ge30_cells_in_gems_with_controls"] += (
            n_controlled >= MIN_TASK_CELLS
        )
        panel_counts["main_panel_targets_ge30_cells_in_gems_with_ge30_controls"] += (
            n_well_supported >= MIN_TASK_CELLS
        )
        task_rows.append((gene_id, symbol, designed[gene_id], len(task_pairs[gene_id]),
                          n_cells, len(by_gem), n_controlled, n_well_supported))

    design.write_table(
        out / "gem_control_counts.tsv",
        ("gem", "guide_csv_cells", "non_targeting_cells", "observed_control_pairs",
         "targeting_cells", "main_panel_cells"),
        ((gem, cells[gem], controls[gem], len(control_pairs[gem]), targeting[gem],
          panel_cells_by_gem[gem]) for gem in sorted(cells)),
    )
    design.write_table(
        out / "task_support.tsv",
        ("gene_id", "panel_symbol", "designed_guide_pairs", "observed_guide_pairs",
         "guide_validated_cells", "observed_gems", "cells_in_gems_with_controls",
         "cells_in_gems_with_ge30_controls"),
        task_rows,
    )
    design.write_table(
        out / "task_gem_counts.tsv",
        ("gene_id", "gem", "guide_validated_cells", "non_targeting_cells_in_gem"),
        ((gene_id, gem, n, controls[gem]) for gene_id in sorted(task_gems)
         for gem, n in sorted(task_gems[gene_id].items())),
    )
    design.write_table(
        out / "task_guide_counts.tsv",
        ("gene_id", "guide_pair", "guide_validated_cells"),
        ((gene_id, "|".join(pair), n) for gene_id in sorted(task_pairs)
         for pair, n in sorted(task_pairs[gene_id].items())),
    )

    return {
        "status": "GUIDE_METADATA_ONLY",
        "cell_line": cell_line,
        "unique_cells_in_guide_csv": total,
        "excluded_cells_by_first_reason": dict(sorted(excluded.items())),
        "gems_in_guide_csv": len(cells),
        "non_targeting_cells": sum(controls.values()),
        "min_non_targeting_cells_per_gem": min(controls[gem] for gem in cells),
        "gems_with_zero_non_targeting_cells": sum(controls[gem] == 0 for gem in cells),
        "gems_with_lt30_non_targeting_cells": sum(controls[gem] < CONTROL_DIAGNOSTIC for gem in cells),
        "targeting_cells": sum(targeting.values()),
        "observed_target_gene_ids": len(observed_ids),
        "main_panel_tasks": len(panel),
        "main_panel_targets_in_design": len(designed),
        "main_panel_cells": sum(panel_cells_by_gem.values()),
        **panel_counts,
    }


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-line", required=True)
    parser.add_argument("--guide-calls", type=Path, required=True)
    parser.add_argument("--guide-library", type=Path, required=True)
    parser.add_argument("--task-panel", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    panel, _ = design.panel_ids(args.task_panel, 1264)
    library = read_library(args.guide_library)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    print(f"Scanning {args.guide_calls.name}", flush=True)
    counts = scan_calls(args.guide_calls, library, panel, args.cell_line)
    summary = write_results(args.output_dir, panel, library, counts, args.cell_line)
    save_json(args.output_dir / "summary.json", summary)
    save_json(args.output_dir / "parameters.json", {
        "cell_line": args.cell_line,
        "guide_calls": str(args.guide_calls),
        "guide_library": str(args.guide_library),
        "task_panel": str(args.task_panel),
        "python_version": sys.version.split()[0],
        "cell_key": "Unique full barcode. GEM is the sample suffix after the first hyphen.",
        "guide_rule": "True/1 filter, two guides, exact unordered library pair, matching target label.",
        "target_id": "Use the guide library Ensembl ID. The cell gene label serves as a check.",
        "control_rule": "Library target_gene_id Non-Targeting with a case-insensitive label check.",
        "diagnostic_min_task_cells": MIN_TASK_CELLS,
        "diagnostic_min_controls_per_gem": CONTROL_DIAGNOSTIC,
    })
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
