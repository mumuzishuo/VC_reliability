#!/usr/bin/env python3
"""Inspect existing response HDF5 files and compare their axes to a task panel."""

import argparse
import csv
import json
from pathlib import Path

import h5py


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_references(feature_path, task_path):
    features = read_tsv(feature_path)
    tasks = read_tsv(task_path)
    main = [row for row in tasks if row["in_main_panel"] == "1"]
    refs = {
        "expression_gene_id": [row["gene_id"] for row in features],
        "main_target_symbol": [row["gene"] for row in main],
        "main_target_gene_id": [row["gene_id"] for row in main],
        "shared_target_symbol": [row["gene"] for row in tasks],
        "shared_target_gene_id": [row["gene_id"] for row in tasks],
    }
    for name, values in refs.items():
        if not values or any(not value for value in values) or len(set(values)) != len(values):
            raise ValueError(f"The reference axis {name} is empty or has invalid IDs")
    return refs


def describe_axis(labels, refs, sets):
    unique = set(labels)
    duplicates = len(labels) - len(unique)
    overlaps = {name: len(unique & expected) for name, expected in sets.items()}
    best = max(
        refs,
        key=lambda name: (
            overlaps[name] / len(unique | sets[name]),
            labels == refs[name],
        ),
    )
    if not overlaps[best]:
        return {
            "length": len(labels),
            "duplicate_labels": duplicates,
            "reference": None,
            "status": "no_overlap_with_references",
        }

    missing = len(sets[best] - unique)
    extra = len(unique - sets[best])
    if duplicates:
        status = "duplicate_labels"
    elif labels == refs[best]:
        status = "same_ids_same_order"
    elif missing == extra == 0:
        status = "same_ids_different_order"
    elif missing == 0:
        status = "contains_reference_with_extras"
    else:
        status = "missing_reference_ids"

    record = {
        "length": len(labels),
        "duplicate_labels": duplicates,
        "reference": best,
        "status": status,
        "shared": overlaps[best],
        "missing": missing,
        "extra": extra,
    }
    if best.endswith("_symbol"):
        record["missing_main_tasks"] = len(sets["main_target_symbol"] - unique)
    if best.endswith("_gene_id") and best != "expression_gene_id":
        record["missing_main_tasks"] = len(sets["main_target_gene_id"] - unique)
    return record


def inspect_file(path, refs, sets):
    label_limit = 2 * max(len(ref) for ref in refs.values())
    datasets = {}
    string_axes = {}
    with h5py.File(path, "r") as handle:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                datasets[name] = {
                    "shape": list(obj.shape) if obj.shape is not None else None,
                    "dtype": str(obj.dtype),
                    "attributes": sorted(obj.attrs.keys()),
                }
                if (
                    obj.ndim == 1
                    and 0 < obj.size <= label_limit
                    and h5py.check_string_dtype(obj.dtype) is not None
                ):
                    labels = obj.asstr(errors="replace")[...].tolist()
                    string_axes[name] = describe_axis(labels, refs, sets)

        handle.visititems(visit)
        root_attributes = sorted(handle.attrs.keys())

    return {
        "file": path.name,
        "file_bytes": path.stat().st_size,
        "root_attributes": root_attributes,
        "datasets": datasets,
        "string_axes": string_axes,
    }


def parse_input(value):
    label, separator, path = value.partition("=")
    if not separator or not label or not path:
        raise argparse.ArgumentTypeError("Input must use ENVIRONMENT=PATH")
    return label, Path(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, type=parse_input)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--tasks", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if len({label for label, _ in args.input}) != len(args.input):
        parser.error("Environment labels must be unique")

    refs = read_references(args.features, args.tasks)
    sets = {name: set(values) for name, values in refs.items()}
    report = {
        "reference_lengths": {name: len(values) for name, values in refs.items()},
        "environments": {},
    }
    for label, path in args.input:
        report["environments"][label] = inspect_file(path, refs, sets)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    for label, info in report["environments"].items():
        print(f"\n{label}  {info['file']}")
        print("  Root attributes", ", ".join(info["root_attributes"]) or "(none)")
        for name, dataset in info["datasets"].items():
            print(f"  {name}  shape={dataset['shape']}  dtype={dataset['dtype']}")
        for name, axis in info["string_axes"].items():
            print(
                f"  AXIS {name}  {axis['reference']}  {axis['status']}  "
                f"n={axis['length']}  missing={axis.get('missing', '-')}  "
                f"extra={axis.get('extra', '-')}  "
                f"missing_main_tasks={axis.get('missing_main_tasks', '-')}"
            )


if __name__ == "__main__":
    main()
