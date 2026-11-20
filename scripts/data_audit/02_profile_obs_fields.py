#!/usr/bin/env python3
"""Profile perturbation and cell-level metadata fields from AnnData HDF5 files."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import h5py
import numpy as np


DEFAULT_COLUMNS = [
    "sgID_AB",
    "gene",
    "gene_id",
    "gem_group",
    "UMI_count",
    "mitopercent",
    "batch",
    "sample",
    "donor",
]


def normalise_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def decode_scalar(value: Any) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        value = value.item()
    return str(value)


def read_column(node: h5py.Dataset | h5py.Group) -> list[str] | None:
    """Decode ordinary, pandas-categorical and AnnData categorical columns."""
    if isinstance(node, h5py.Dataset):
        raw_values = np.asarray(node[...])

        # AnnData stores categorical observation values as integer codes.
        # The categories dataset is referenced from the dataset attributes.
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
def find_column(obs: h5py.Group, requested: str) -> str | None:
    requested_name = normalise_name(requested)
    for actual_name in obs.keys():
        if normalise_name(actual_name) == requested_name:
            return actual_name
    return None


def is_missing(value: str) -> bool:
    return value in {"", "<NA>", "NA", "NaN", "nan", "None", "none"}


def numeric_summary(values: list[str]) -> dict[str, float] | None:
    numeric_values = []
    for value in values:
        if is_missing(value):
            continue
        try:
            numeric_values.append(float(value))
        except ValueError:
            return None

    if not numeric_values:
        return None

    array = np.asarray(numeric_values, dtype=float)
    return {
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
        "mean": float(np.mean(array)),
        "median": float(np.median(array)),
    }


def summarise_values(values: list[str], top_n: int) -> dict[str, Any]:
    counts = Counter(values)
    return {
        "n_values": len(values),
        "n_missing": sum(is_missing(value) for value in values),
        "n_unique": len(counts),
        "top_values": [
            {"value": value, "count": int(count)}
            for value, count in counts.most_common(top_n)
        ],
        "numeric_summary": numeric_summary(values),
    }


def mapping_summary(
    labels: list[str] | None,
    genes: list[str] | None,
    gene_ids: list[str] | None,
    top_n: int,
) -> dict[str, Any] | None:
    if labels is None or genes is None or len(labels) != len(genes):
        return None

    label_gene_pairs = Counter(zip(labels, genes))
    genes_by_label: dict[str, set[str]] = defaultdict(set)
    labels_by_gene: dict[str, set[str]] = defaultdict(set)

    for label, gene in label_gene_pairs:
        if is_missing(label) or is_missing(gene):
            continue
        genes_by_label[label].add(gene)
        labels_by_gene[gene].add(label)

    result = {
        "n_unique_label_gene_pairs": len(label_gene_pairs),
        "n_labels_with_multiple_gene_values": sum(
            len(values) > 1 for values in genes_by_label.values()
        ),
        "n_genes_with_multiple_labels": sum(
            len(values) > 1 for values in labels_by_gene.values()
        ),
        "top_label_gene_pairs": [
            {
                "label": label,
                "gene": gene,
                "count": int(count),
            }
            for (label, gene), count in label_gene_pairs.most_common(top_n)
        ],
    }

    if gene_ids is not None and len(labels) == len(gene_ids):
        label_gene_id_pairs = Counter(zip(labels, gene_ids))
        result["n_unique_label_gene_id_pairs"] = len(label_gene_id_pairs)

    return result


def inspect_file(label: str, path: Path, requested_columns: list[str], top_n: int) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        obs = handle.get("obs")
        if not isinstance(obs, h5py.Group):
            raise ValueError(f"Missing obs group in {path}")

        resolved = {
            requested: find_column(obs, requested)
            for requested in requested_columns
        }
        values = {}
        summaries = {}

        for requested, actual in resolved.items():
            if actual is None:
                continue
            column_values = read_column(obs[actual])
            if column_values is None:
                continue
            values[requested] = column_values
            summaries[requested] = {
                "actual_column": actual,
                **summarise_values(column_values, top_n),
            }

        return {
            "label": label,
            "file_name": path.name,
            "file_size_bytes": path.stat().st_size,
            "obs_columns": sorted(obs.keys()),
            "resolved_columns": resolved,
            "column_summaries": summaries,
            "sgID_AB_mapping": mapping_summary(
                values.get("sgID_AB"),
                values.get("gene"),
                values.get("gene_id"),
                top_n,
            ),
        }


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
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--columns",
        default=",".join(DEFAULT_COLUMNS),
        help="Comma-separated obs columns to profile",
    )
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    requested_columns = [
        column.strip()
        for column in args.columns.split(",")
        if column.strip()
    ]

    summaries = []
    for label, path in args.input:
        if not path.is_file():
            raise FileNotFoundError(path)
        summaries.append(
            inspect_file(label, path, requested_columns, args.top_n)
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output_handle:
        json.dump(summaries, output_handle, indent=2, ensure_ascii=False)
        output_handle.write("\n")


if __name__ == "__main__":
    main()
