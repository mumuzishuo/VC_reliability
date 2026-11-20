#!/usr/bin/env python3
"""Inspect AnnData HDF5 structure without loading the expression matrix."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import h5py
import numpy as np


COLUMN_GROUPS = {
    "perturbation": [
        "sgidab",
        "perturbation",
        "perturbedgene",
        "targetgene",
        "gene",
        "genename",
    ],
    "gem_group": ["gemgroup", "gem", "gemgroupid"],
    "batch": ["batch", "batchid"],
    "sample": ["sample", "sampleid", "library"],
    "donor": ["donor", "donorid"],
    "umi": ["umicount", "ncounts", "totalcounts", "librarysize"],
    "mitochondrial_fraction": ["mitopercent", "mitochondrialfraction", "pctmt"],
}


def normalise_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def decode_scalar(value: Any) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        value = value.item()
    return str(value)


def read_column(node: h5py.Dataset | h5py.Group) -> list[str] | None:
    if isinstance(node, h5py.Dataset):
        return [decode_scalar(value) for value in node[...]]

    if "codes" not in node or "categories" not in node:
        return None

    codes = np.asarray(node["codes"][...])
    categories = np.asarray(node["categories"][...])
    values = []
    for code in codes:
        code_value = int(code)
        if code_value < 0:
            values.append("<NA>")
        else:
            values.append(decode_scalar(categories[code_value]))
    return values


def get_shape(node: h5py.Dataset | h5py.Group | None) -> list[int] | None:
    if node is None:
        return None

    if isinstance(node, h5py.Dataset):
        return [int(value) for value in node.shape]

    if "shape" in node.attrs:
        return [int(value) for value in node.attrs["shape"]]

    if "shape" in node:
        shape_node = node["shape"]
        if isinstance(shape_node, h5py.Dataset):
            return [int(value) for value in shape_node[...]]

    return None


def candidate_columns(group: h5py.Group | None) -> dict[str, str]:
    if group is None:
        return {}

    available = {normalise_name(name): name for name in group.keys()}
    matches = {}
    for category, candidates in COLUMN_GROUPS.items():
        for candidate in candidates:
            if candidate in available:
                matches[category] = available[candidate]
                break
    return matches


def summarise_column(
    group: h5py.Group,
    column_name: str,
    top_n: int,
) -> dict[str, Any]:
    values = read_column(group[column_name])
    if values is None:
        return {
            "column": column_name,
            "readable_values": False,
        }

    counts = Counter(values)
    top_values = [
        {"value": value, "count": int(count)}
        for value, count in counts.most_common(top_n)
    ]
    return {
        "column": column_name,
        "readable_values": True,
        "n_values": len(values),
        "n_unique": len(counts),
        "top_values": top_values,
    }


def inspect_file(label: str, path: Path, top_n: int) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        obs = handle.get("obs")
        var = handle.get("var")
        x_node = handle.get("X")
        matches = candidate_columns(obs if isinstance(obs, h5py.Group) else None)

        column_summaries = {}
        if isinstance(obs, h5py.Group):
            for category, column_name in matches.items():
                column_summaries[category] = summarise_column(
                    obs, column_name, top_n
                )

        return {
            "label": label,
            "file_name": path.name,
            "file_size_bytes": path.stat().st_size,
            "top_level_keys": sorted(handle.keys()),
            "x_shape": get_shape(x_node),
            "x_keys": sorted(x_node.keys()) if isinstance(x_node, h5py.Group) else [],
            "obs_columns": sorted(obs.keys()) if isinstance(obs, h5py.Group) else [],
            "var_columns": sorted(var.keys()) if isinstance(var, h5py.Group) else [],
            "candidate_columns": matches,
            "candidate_column_summaries": column_summaries,
            "has_raw_group": "raw" in handle,
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
    parser.add_argument("--top-n", type=int, default=10)
    args = parser.parse_args()

    summaries = []
    for label, path in args.input:
        if not path.is_file():
            raise FileNotFoundError(path)
        summaries.append(inspect_file(label, path, args.top_n))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output_handle:
        json.dump(summaries, output_handle, indent=2, ensure_ascii=False)
        output_handle.write("\n")


if __name__ == "__main__":
    main()
