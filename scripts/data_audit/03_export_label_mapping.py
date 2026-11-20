#!/usr/bin/env python3
"""Export decoded guide-to-gene and cross-environment mapping tables."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import h5py
import numpy as np


REQUIRED_COLUMNS = ["sgID_AB", "gene", "gene_id", "gem_group"]


def decode_scalar(value: Any) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        value = value.item()
    return str(value)


def read_column(node: h5py.Dataset | h5py.Group) -> list[str] | None:
    if isinstance(node, h5py.Dataset):
        raw_values = np.asarray(node[...])

        if "categories" in node.attrs:
            reference = node.attrs["categories"]
            if isinstance(reference, np.ndarray):
                reference = reference.item()
            categories_node = node.file[reference]
            categories = np.asarray(categories_node[...])

            values = []
            for code in raw_values:
                code_value = int(code)
                if code_value < 0 or code_value >= len(categories):
                    values.append("<NA>")
                else:
                    values.append(decode_scalar(categories[code_value]))
            return values

        return [decode_scalar(value) for value in raw_values]

    if "codes" not in node or "categories" not in node:
        return None

    codes = np.asarray(node["codes"][...])
    categories = np.asarray(node["categories"][...])
    values = []
    for code in codes:
        code_value = int(code)
        if code_value < 0 or code_value >= len(categories):
            values.append("<NA>")
        else:
            values.append(decode_scalar(categories[code_value]))
    return values


def is_control(gene: str) -> bool:
    normalised = re.sub(r"[^a-z0-9]", "", gene.lower())
    return normalised in {"nontargeting", "negativecontrol", "vehicle"}


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def profile_dataset(
    label: str,
    path: Path,
    output_dir: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        values = {}
        for column in REQUIRED_COLUMNS:
            if column not in obs:
                raise KeyError(f"{column} is missing from {path}")
            decoded = read_column(obs[column])
            if decoded is None:
                raise ValueError(f"{column} could not be decoded in {path}")
            values[column] = decoded

    n_cells = len(values["gene"])
    if any(len(column_values) != n_cells for column_values in values.values()):
        raise ValueError(f"Metadata columns have different lengths in {path}")

    guide_counts: Counter[tuple[str, str, str]] = Counter()
    guide_gems: defaultdict[tuple[str, str, str], set[str]] = defaultdict(set)
    gene_counts: Counter[tuple[str, str]] = Counter()
    gene_gems: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    gene_guides: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    gene_ids_by_gene: defaultdict[str, set[str]] = defaultdict(set)

    for guide, gene, gene_id, gem in zip(
        values["sgID_AB"],
        values["gene"],
        values["gene_id"],
        values["gem_group"],
    ):
        guide_key = (guide, gene, gene_id)
        gene_key = (gene, gene_id)

        guide_counts[guide_key] += 1
        guide_gems[guide_key].add(gem)

        gene_counts[gene_key] += 1
        gene_gems[gene_key].add(gem)
        gene_guides[gene_key].add(guide)
        gene_ids_by_gene[gene].add(gene_id)

    guide_rows = []
    for (guide, gene, gene_id), cell_count in sorted(
        guide_counts.items(), key=lambda item: (item[0][1], item[0][0])
    ):
        guide_rows.append(
            [
                label,
                guide,
                gene,
                gene_id,
                str(cell_count),
                str(len(guide_gems[(guide, gene, gene_id)])),
                str(is_control(gene)),
            ]
        )

    gene_rows = []
    for (gene, gene_id), cell_count in sorted(gene_counts.items()):
        gene_rows.append(
            [
                label,
                gene,
                gene_id,
                str(cell_count),
                str(len(gene_guides[(gene, gene_id)])),
                str(len(gene_gems[(gene, gene_id)])),
                str(is_control(gene)),
            ]
        )

    write_tsv(
        output_dir / f"guide_to_gene_{label}.tsv",
        [
            "environment",
            "sgID_AB",
            "gene",
            "gene_id",
            "cell_count",
            "gem_group_count",
            "is_control",
        ],
        guide_rows,
    )
    write_tsv(
        output_dir / f"gene_summary_{label}.tsv",
        [
            "environment",
            "gene",
            "gene_id",
            "cell_count",
            "guide_count",
            "gem_group_count",
            "is_control",
        ],
        gene_rows,
    )

    guide_records = [
        {
            "environment": row[0],
            "sgID_AB": row[1],
            "gene": row[2],
            "gene_id": row[3],
            "cell_count": row[4],
            "gem_group_count": row[5],
            "is_control": row[6],
        }
        for row in guide_rows
    ]
    gene_records = [
        {
            "environment": row[0],
            "gene": row[1],
            "gene_id": row[2],
            "cell_count": row[3],
            "guide_count": row[4],
            "gem_group_count": row[5],
            "is_control": row[6],
        }
        for row in gene_rows
    ]
    return guide_records, gene_records


def write_cross_environment_check(
    gene_records: list[dict[str, str]],
    output_dir: Path,
) -> None:
    environments_by_gene: defaultdict[str, set[str]] = defaultdict(set)
    gene_ids_by_gene: defaultdict[str, set[str]] = defaultdict(set)

    for record in gene_records:
        gene = record["gene"]
        environments_by_gene[gene].add(record["environment"])
        gene_ids_by_gene[gene].add(record["gene_id"])

    rows = []
    for gene in sorted(environments_by_gene):
        environments = sorted(environments_by_gene[gene])
        gene_ids = sorted(gene_ids_by_gene[gene])
        status = "consistent" if len(gene_ids) == 1 else "multiple_gene_ids"
        rows.append(
            [
                gene,
                "|".join(environments),
                str(len(environments)),
                "|".join(gene_ids),
                str(len(gene_ids)),
                status,
            ]
        )

    write_tsv(
        output_dir / "cross_environment_gene_id_check.tsv",
        [
            "gene",
            "environments",
            "n_environments",
            "gene_ids",
            "n_gene_ids",
            "status",
        ],
        rows,
    )


def parse_input(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Input must use LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("Input must use LABEL=PATH")
    return label, Path(raw_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        type=parse_input,
        help="Dataset label and path in the form LABEL=PATH",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    all_guides = []
    all_genes = []

    for label, path in args.input:
        if not path.is_file():
            raise FileNotFoundError(path)
        guides, genes = profile_dataset(label, path, args.output_dir)
        all_guides.extend(guides)
        all_genes.extend(genes)

    write_tsv(
        args.output_dir / "guide_to_gene_all_environments.tsv",
        [
            "environment",
            "sgID_AB",
            "gene",
            "gene_id",
            "cell_count",
            "gem_group_count",
            "is_control",
        ],
        [
            [
                record["environment"],
                record["sgID_AB"],
                record["gene"],
                record["gene_id"],
                record["cell_count"],
                record["gem_group_count"],
                record["is_control"],
            ]
            for record in all_guides
        ],
    )
    write_tsv(
        args.output_dir / "gene_summary_all_environments.tsv",
        [
            "environment",
            "gene",
            "gene_id",
            "cell_count",
            "guide_count",
            "gem_group_count",
            "is_control",
        ],
        [
            [
                record["environment"],
                record["gene"],
                record["gene_id"],
                record["cell_count"],
                record["guide_count"],
                record["gem_group_count"],
                record["is_control"],
            ]
            for record in all_genes
        ],
    )
    write_cross_environment_check(all_genes, args.output_dir)


if __name__ == "__main__":
    main()
