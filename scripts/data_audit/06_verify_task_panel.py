#!/usr/bin/env python3
"""Reproduce the four-environment perturbation panel and audit GEM support."""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


ENVIRONMENTS = ("HepG2", "Jurkat", "K562_essential", "RPE1")
CELL_THRESHOLDS = (10, 20, 30, 50, 100)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def read_gene_counts(mapping_dir):
    targets = {}
    controls = Counter()
    for env in ENVIRONMENTS:
        targets[env] = {}
        for row in read_tsv(mapping_dir / f"gene_summary_{env}.tsv"):
            gene, gene_id = row["gene"], row["gene_id"]
            cells = int(row["cell_count"])
            if row["environment"] != env or cells < 1:
                raise ValueError(f"Unexpected environment or cell count in {env}")
            if row["is_control"] == "True":
                controls[env] += cells
                continue
            if row["is_control"] != "False" or not gene_id or gene in targets[env]:
                raise ValueError(f"Invalid or ambiguous target {gene!r} in {env}")
            targets[env][gene] = (gene_id, cells)
    return targets, controls


def build_panel(targets, min_cells):
    shared_genes = set.intersection(*(set(targets[env]) for env in ENVIRONMENTS))
    panel = []
    for gene in sorted(shared_genes):
        gene_ids, cell_counts = zip(*(targets[env][gene] for env in ENVIRONMENTS))
        if len(set(gene_ids)) != 1:
            raise ValueError(f"Inconsistent gene IDs for shared target {gene!r}")
        smallest = min(cell_counts)
        panel.append(
            (gene, gene_ids[0], *cell_counts, smallest, int(smallest >= min_cells))
        )
    if len({row[1] for row in panel}) != len(panel):
        raise ValueError("Multiple shared targets use the same gene_id")
    return panel


def audit_gem_support(gem_dir, targets, expected_controls, panel):
    gem_counts = {}
    control_totals = Counter()
    for row in read_tsv(gem_dir / "gem_control_summary_all_environments.tsv"):
        env, gem = row["environment"], row["gem_group"]
        controls = int(row["control_cells"])
        perturbed = int(row["non_control_cells"])
        total = int(row["total_cells"])
        if (
            env not in targets or (env, gem) in gem_counts or total < 1
            or min(controls, perturbed) < 0 or controls + perturbed != total
        ):
            raise ValueError(f"Invalid GEM summary row for {env}/{gem}")
        gem_counts[env, gem] = (controls, perturbed)
        control_totals[env] += controls
    if control_totals != expected_controls:
        raise ValueError("GEM control counts differ from the mapping totals")

    main_targets = [(row[0], row[1]) for row in panel if row[-1]]
    main_genes = {gene for gene, _ in main_targets}
    expected = {
        (env, gene, gene_id): cells
        for env, genes in targets.items()
        for gene, (gene_id, cells) in genes.items()
    }
    observed = Counter()
    observed_by_gem = Counter()
    min_controls = {}
    groups_used = Counter()
    cells_below = {threshold: Counter() for threshold in (1, 10, 20)}
    seen = set()
    for row in read_tsv(gem_dir / "perturbation_cell_counts_by_gem_all_environments.tsv"):
        env, gem = row["environment"], row["gem_group"]
        gene, gene_id = row["gene"], row["gene_id"]
        key = (env, gene, gene_id)
        cells = int(row["perturbation_cells"])
        if (env, gem) not in gem_counts or key not in expected:
            raise ValueError(f"Unmapped target or GEM for {env}/{gene}/{gem}")
        controls = gem_counts[env, gem][0]
        if cells < 1 or int(row["control_cells_same_gem"]) != controls or (*key, gem) in seen:
            raise ValueError(f"Invalid or duplicate perturbation row for {env}/{gene}/{gem}")
        seen.add((*key, gem))
        observed[key] += cells
        observed_by_gem[env, gem] += cells
        if gene in main_genes:
            min_controls[key] = min(min_controls.get(key, controls), controls)
            groups_used[key] += 1
            for threshold, counter in cells_below.items():
                if controls < threshold:
                    counter[key] += cells

    if observed != expected:
        changed = sum(
            observed[key] != expected.get(key, 0) for key in observed.keys() | expected.keys()
        )
        raise ValueError(f"{changed} target counts differ between GEM and mapping tables")
    if any(observed_by_gem[key] != perturbed for key, (_, perturbed) in gem_counts.items()):
        raise ValueError("Perturbation counts do not sum to the GEM totals")

    low_groups = []
    support = []
    summary = {}
    for env in ENVIRONMENTS:
        counts = [n for (label, _), (n, _) in gem_counts.items() if label == env]
        low_groups.extend(
            (env, gem, controls, perturbed)
            for (label, gem), (controls, perturbed) in sorted(gem_counts.items())
            if label == env and controls < 20
        )
        affected_lt10 = 0
        affected_lt20 = 0
        max_fraction_lt20 = 0.0
        for gene, gene_id in main_targets:
            key = (env, gene, gene_id)
            total = expected[key]
            n_lt10 = cells_below[10][key]
            n_lt20 = cells_below[20][key]
            affected_lt10 += n_lt10 > 0
            affected_lt20 += n_lt20 > 0
            fraction_lt20 = n_lt20 / total
            max_fraction_lt20 = max(max_fraction_lt20, fraction_lt20)
            support.append((
                env, gene, gene_id, total, groups_used[key], min_controls[key],
                cells_below[1][key], n_lt10, n_lt20, n_lt10 / total, fraction_lt20,
            ))
        summary[env] = {
            "gem_groups": len(counts),
            "zero_control_gem_groups": counts.count(0),
            "gem_groups_lt10_controls": sum(n < 10 for n in counts),
            "gem_groups_lt20_controls": sum(n < 20 for n in counts),
            "main_tasks_in_gems_lt10_controls": affected_lt10,
            "main_tasks_in_gems_lt20_controls": affected_lt20,
            "max_perturbation_cell_fraction_in_gems_lt20_controls": max_fraction_lt20,
        }
    return low_groups, support, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping-dir", required=True, type=Path)
    parser.add_argument("--gem-dir", type=Path, help="Include the Step 6 GEM support audit")
    parser.add_argument("--min-cells", default=30, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.min_cells < 1:
        parser.error("--min-cells must be positive")

    targets, controls = read_gene_counts(args.mapping_dir)
    panel = build_panel(targets, args.min_cells)
    summary = {
        "shared_non_control_targets": len(panel),
        "min_cells_in_every_environment": args.min_cells,
        "main_panel_targets": sum(row[-1] for row in panel),
        "targets_by_min_cells_in_every_environment": {
            n: sum(row[-2] >= n for row in panel) for n in CELL_THRESHOLDS
        },
    }
    if args.gem_dir:
        low_groups, support, summary["gem_control_diagnostics"] = audit_gem_support(
            args.gem_dir, targets, controls, panel
        )

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(
        out / "shared_perturbation_tasks.tsv",
        ("gene", "gene_id", *ENVIRONMENTS, "min_cells", "in_main_panel"),
        panel,
    )
    if args.gem_dir:
        write_tsv(
            out / "low_control_gems.tsv",
            ("environment", "gem_group", "control_cells", "non_control_cells"),
            low_groups,
        )
        write_tsv(
            out / "main_task_control_support.tsv",
            (
                "environment", "gene", "gene_id", "perturbation_cells",
                "gems_used", "min_controls_in_used_gems", "cells_in_zero_control_gems",
                "cells_in_gems_lt10_controls", "cells_in_gems_lt20_controls",
                "fraction_in_gems_lt10_controls", "fraction_in_gems_lt20_controls",
            ),
            support,
        )

    report = json.dumps(summary, indent=2, allow_nan=False) + "\n"
    (out / "task_panel_summary.json").write_text(report, encoding="utf-8")
    print(report, end="")


if __name__ == "__main__":
    main()
