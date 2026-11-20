#!/usr/bin/env python3
"""Profile perturbation and control cell counts by GEM group."""

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
) -> tuple[list[list[str]], list[list[str]]]:
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

    gem_total: Counter[str] = Counter()
    gem_control: Counter[str] = Counter()
    gem_non_control_genes: defaultdict[str, set[str]] = defaultdict(set)
    gem_non_control_guides: defaultdict[str, set[str]] = defaultdict(set)

    gene_gem_counts: Counter[tuple[str, str, str]] = Counter()
    gene_total_counts: Counter[tuple[str, str]] = Counter()
    gene_gems: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    gene_guides: defaultdict[tuple[str, str], set[str]] = defaultdict(set)

    for guide, gene, gene_id, gem in zip(
        values["sgID_AB"],
        values["gene"],
        values["gene_id"],
        values["gem_group"],
    ):
        control = is_control(gene)
        gem_total[gem] += 1
        if control:
            gem_control[gem] += 1
        else:
            gem_non_control_genes[gem].add(gene)
            gem_non_control_guides[gem].add(guide)
            gene_gem_counts[(gene, gene_id, gem)] += 1
            gene_total_counts[(gene, gene_id)] += 1
            gene_gems[(gene, gene_id)].add(gem)
            gene_guides[(gene, gene_id)].add(guide)

    gem_rows = []
    for gem in sorted(gem_total, key=lambda value: int(value)):
        total = gem_total[gem]
        control = gem_control[gem]
        non_control = total - control
        gem_rows.append(
            [
                label,
                gem,
                str(total),
                str(control),
                str(non_control),
                f"{control / total:.8f}" if total else "0.00000000",
                str(len(gem_non_control_genes[gem])),
                str(len(gem_non_control_guides[gem])),
            ]
        )

    perturbation_rows = []
    for (gene, gene_id, gem), cell_count in sorted(gene_gem_counts.items()):
        control_count = gem_control[gem]
        perturbation_rows.append(
            [
                label,
                gene,
                gene_id,
                gem,
                str(cell_count),
                str(control_count),
                f"{control_count / (cell_count + control_count):.8f}"
                if cell_count + control_count
                else "0.00000000",
            ]
        )

    gene_rows = []
    for (gene, gene_id), cell_count in sorted(gene_total_counts.items()):
        gene_rows.append(
            [
                label,
                gene,
                gene_id,
                str(cell_count),
                str(len(gene_guides[(gene, gene_id)])),
                str(len(gene_gems[(gene, gene_id)])),
            ]
        )

    write_tsv(
        output_dir / f"gem_control_summary_{label}.tsv",
        [
            "environment",
            "gem_group",
            "total_cells",
            "control_cells",
            "non_control_cells",
            "control_fraction",
            "non_control_gene_count",
            "non_control_guide_count",
        ],
        gem_rows,
    )
    write_tsv(
        output_dir / f"perturbation_cell_counts_by_gem_{label}.tsv",
        [
            "environment",
            "gene",
            "gene_id",
            "gem_group",
            "perturbation_cells",
            "control_cells_same_gem",
            "control_fraction_in_pair",
        ],
        perturbation_rows,
    )
    write_tsv(
        output_dir / f"perturbation_cell_counts_{label}.tsv",
        [
            "environment",
            "gene",
            "gene_id",
            "total_perturbation_cells",
            "guide_count",
            "gem_group_count",
        ],
        gene_rows,
    )

    return gem_rows, perturbation_rows


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

    all_gem_rows = []
    all_perturbation_rows = []

    for label, path in args.input:
        if not path.is_file():
            raise FileNotFoundError(path)
        gem_rows, perturbation_rows = profile_dataset(
            label, path, args.output_dir
        )
        all_gem_rows.extend(gem_rows)
        all_perturbation_rows.extend(perturbation_rows)

    write_tsv(
        args.output_dir / "gem_control_summary_all_environments.tsv",
        [
            "environment",
            "gem_group",
            "total_cells",
            "control_cells",
            "non_control_cells",
            "control_fraction",
            "non_control_gene_count",
            "non_control_guide_count",
        ],
        all_gem_rows,
    )
    write_tsv(
        args.output_dir / "perturbation_cell_counts_by_gem_all_environments.tsv",
        [
            "environment",
            "gene",
            "gene_id",
            "gem_group",
            "perturbation_cells",
            "control_cells_same_gem",
            "control_fraction_in_pair",
        ],
        all_perturbation_rows,
    )


if __name__ == "__main__":
    main()
