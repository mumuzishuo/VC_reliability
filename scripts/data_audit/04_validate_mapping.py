#!/usr/bin/env python3
"""Validate exported guide-to-gene and cross-environment mapping tables."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_control(gene: str) -> bool:
    normalised = re.sub(r"[^a-z0-9]", "", gene.lower())
    return normalised in {"nontargeting", "negativecontrol", "vehicle"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    mapping_dir = args.mapping_dir
    gene_rows = read_tsv(mapping_dir / "gene_summary_all_environments.tsv")
    guide_rows = read_tsv(mapping_dir / "guide_to_gene_all_environments.tsv")
    cross_rows = read_tsv(mapping_dir / "cross_environment_gene_id_check.tsv")

    genes_by_environment: defaultdict[str, set[str]] = defaultdict(set)
    non_control_genes_by_environment: defaultdict[str, set[str]] = defaultdict(set)
    guides_by_environment: defaultdict[str, set[str]] = defaultdict(set)
    cells_by_environment: Counter[str] = Counter()
    control_cells_by_environment: Counter[str] = Counter()

    for row in gene_rows:
        environment = row["environment"]
        gene = row["gene"]
        genes_by_environment[environment].add(gene)
        cells_by_environment[environment] += int(row["cell_count"])
        if is_control(gene):
            control_cells_by_environment[environment] += int(row["cell_count"])
        else:
            non_control_genes_by_environment[environment].add(gene)

    for row in guide_rows:
        guides_by_environment[row["environment"]].add(row["sgID_AB"])

    coverage_distribution = Counter()
    consistent_distribution = Counter()
    multiple_gene_id_genes = []

    for row in cross_rows:
        gene = row["gene"]
        if is_control(gene):
            continue

        n_environments = int(row["n_environments"])
        coverage_distribution[n_environments] += 1

        if row["status"] == "consistent":
            consistent_distribution[n_environments] += 1
        else:
            multiple_gene_id_genes.append(gene)

    shared_all_environments = consistent_distribution[4]

    summary = {
        "environment_gene_category_counts": {
            environment: len(genes)
            for environment, genes in sorted(genes_by_environment.items())
        },
        "environment_non_control_gene_counts": {
            environment: len(genes)
            for environment, genes in sorted(
                non_control_genes_by_environment.items()
            )
        },
        "environment_guide_counts": {
            environment: len(guides)
            for environment, guides in sorted(guides_by_environment.items())
        },
        "gene_summary_cell_counts": dict(sorted(cells_by_environment.items())),
        "control_cell_counts": dict(
            sorted(control_cells_by_environment.items())
        ),
        "non_control_gene_coverage_distribution": dict(
            sorted(coverage_distribution.items())
        ),
        "consistent_gene_id_coverage_distribution": dict(
            sorted(consistent_distribution.items())
        ),
        "shared_non_control_genes_all_four_environments": shared_all_environments,
        "genes_with_multiple_gene_ids": sorted(set(multiple_gene_id_genes)),
        "n_genes_with_multiple_gene_ids": len(set(multiple_gene_id_genes)),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
