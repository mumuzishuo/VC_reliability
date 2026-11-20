#!/usr/bin/env python3
"""Audit expression feature IDs and map shared IDs to original matrix columns."""

import argparse
import csv
import json
import re
from collections import Counter
from importlib import import_module
from pathlib import Path

import h5py


# The earlier audits handle dense and sparse matrix metadata and both
# categorical encodings used in the project's H5AD files.
structure = import_module("01_inspect_h5ad_h5py")
metadata = import_module("02_profile_obs_fields")

MISSING_IDS = {"", "<NA>", "NA", "NaN", "nan", "None", "none"}
ENSEMBL_VERSION = re.compile(r"^(ENSG[0-9]+)\.[0-9]+$")


def index_gene_ids(gene_ids):
    counts = Counter(
        gene_id for gene_id in gene_ids if gene_id.strip() and gene_id not in MISSING_IDS
    )
    exact_indices = {
        gene_id: i for i, gene_id in enumerate(gene_ids) if counts[gene_id] == 1
    }
    versionless_counts = Counter(
        ENSEMBL_VERSION.sub(r"\1", gene_id)
        for gene_id in gene_ids
        if gene_id in counts
    )
    versionless_ids = {
        gene_id for gene_id, count in versionless_counts.items() if count == 1
    }
    details = {
        "missing_gene_id_rows": len(gene_ids) - sum(counts.values()),
        "duplicated_exact_gene_ids": sum(count > 1 for count in counts.values()),
        "versioned_ensembl_gene_ids": sum(
            ENSEMBL_VERSION.fullmatch(gene_id) is not None for gene_id in counts
        ),
        "duplicated_gene_ids_after_version_removal": sum(
            count > 1 for count in versionless_counts.values()
        ),
    }
    return exact_indices, versionless_ids, details


def inspect_file(path):
    with h5py.File(path, "r") as handle:
        x = handle["X"]
        shape = structure.get_shape(x)
        if shape is None and "h5sparse_shape" in x.attrs:
            shape = [int(n) for n in x.attrs["h5sparse_shape"]]
        gene_ids = metadata.read_column(handle["var/gene_id"])
        if gene_ids is None or shape is None or len(shape) != 2 or len(gene_ids) != shape[1]:
            raise ValueError(f"var/gene_id does not align with X columns in {path.name}")

        if isinstance(x, h5py.Dataset):
            layout = "dense"
            data = x
        else:
            layout = metadata.decode_scalar(
                x.attrs.get("encoding-type", x.attrs.get("h5sparse_format", "group"))
            )
            data = x.get("data")
        info = {
            "file_name": path.name,
            "x_shape": shape,
            "x_layout": layout,
            "x_value_dtype": str(data.dtype) if isinstance(data, h5py.Dataset) else None,
            "x_value_chunks": data.chunks if isinstance(data, h5py.Dataset) else None,
        }

    exact_indices, versionless_ids, id_info = index_gene_ids(gene_ids)
    info.update(id_info)
    return exact_indices, versionless_ids, info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, action="append", type=metadata.parse_input)
    parser.add_argument("--expected-count", type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    labels = [label for label, _ in args.input]
    if len(labels) < 2 or len(labels) != len(set(labels)):
        parser.error("Provide at least two distinct environment labels")

    indices_by_environment = {}
    versionless_sets = []
    environment_summary = {}
    for label, path in args.input:
        exact_indices, versionless_ids, info = inspect_file(path)
        indices_by_environment[label] = exact_indices
        versionless_sets.append(versionless_ids)
        environment_summary[label] = info

    common = sorted(set.intersection(*(set(indices) for indices in indices_by_environment.values())))
    summary = {
        "feature_rule": "exact var/gene_id with a unique matrix column in each input",
        "common_expression_genes_exact": len(common),
        "common_expression_genes_without_versions_diagnostic": len(
            set.intersection(*versionless_sets)
        ),
        "expected_common_expression_genes": args.expected_count,
        "matches_expected_count": (
            len(common) == args.expected_count if args.expected_count is not None else None
        ),
        "environments": environment_summary,
    }

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    with (out / "common_expression_features.tsv").open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *(f"{label}_column_0based" for label in labels)])
        for gene_id in common:
            writer.writerow([gene_id, *(indices_by_environment[label][gene_id] for label in labels)])

    report = json.dumps(summary, indent=2) + "\n"
    (out / "expression_feature_summary.json").write_text(report, encoding="utf-8")
    print(report, end="")


if __name__ == "__main__":
    main()
