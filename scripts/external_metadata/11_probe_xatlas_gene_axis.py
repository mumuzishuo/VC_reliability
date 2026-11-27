#!/usr/bin/env python3
"""Match the candidate expression IDs to the Figshare v3 H5AD var axis."""

import argparse
from collections import defaultdict
import csv
from importlib import import_module
import json
from pathlib import Path
import re
import subprocess
import sys

import h5py

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "legacy_v01"))
import vc_preprocess as legacy

ranges = import_module("08_probe_xatlas_expression_index")
CATALOG = Path("metadata/v01_external_catalog/xatlas_figshare_article.json")
ENSEMBL = re.compile(r"ENSG[0-9]+")
ID_FIELDS = ("gene_id", "gene_ids", "ensembl_id", "ensembl_gene_id", "feature_id", "feature_ids")


def inspect_var(var, gene_ids):
    index_name = var.attrs["_index"]
    if isinstance(index_name, bytes):
        index_name = index_name.decode("utf-8")
    fields = [name for name in dict.fromkeys((index_name, *ID_FIELDS)) if name in var]
    expected = set(gene_ids)
    matches = {}
    reports = {}
    n_var = None
    for name in fields:
        print(f"Reading var identifier field {name}", flush=True)
        labels = legacy.column(var, name)
        if n_var is None:
            n_var = len(labels)
        if len(labels) != n_var:
            raise ValueError(f"Unexpected length for var field {name}")
        positions = defaultdict(list)
        for column, label in enumerate(labels):
            if label in expected:
                positions[label].append(column)
        matches[name] = positions
        reports[name] = {
            "exact_ensembl_id_rows": sum(bool(ENSEMBL.fullmatch(label)) for label in labels),
            "candidate_ids_once": sum(len(columns) == 1 for columns in positions.values()),
            "candidate_ids_duplicated": sum(len(columns) > 1 for columns in positions.values()),
            "candidate_ids_missing": len(gene_ids) - len(positions),
        }

    if index_name not in matches or not n_var:
        raise ValueError("Missing or empty var index")
    # Prefer a field with more unique exact matches. Field order breaks ties.
    name = max(fields, key=lambda field: reports[field]["candidate_ids_once"])
    positions = matches[name]
    rows = []
    for gene_id in gene_ids:
        cols = positions.get(gene_id, [])
        if not cols:
            status = "missing"
        elif len(cols) > 1:
            status = "duplicate"
        elif any(mapping[gene_id] != cols for mapping in matches.values() if gene_id in mapping):
            status = "conflict"
        else:
            status = "exact"
        rows.append((gene_id, cols[0] if status == "exact" else "", status))
    n_exact = sum(status == "exact" for _, _, status in rows)
    report = {
        "status": "EXACT_VAR_AXIS_MATCH" if n_exact == len(rows) else "VAR_AXIS_REVIEW_NEEDED",
        "var_index": index_name,
        "n_var_rows": n_var,
        "var_fields_present": sorted(var.keys()),
        "var_identifier_field_used": name,
        "field_diagnostics": reports,
        "n_candidate_genes": len(gene_ids),
        "n_exact_unique": n_exact,
        "n_unresolved": len(rows) - n_exact,
    }
    return report, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-line", required=True, choices=("HCT116", "HEK293T"))
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("Select a new output directory")
    with args.features.open(encoding="utf-8", newline="") as stream:
        genes = [row["gene_id"] for row in csv.DictReader(stream, delimiter="\t")]
    if len(genes) != 6638 or len(set(genes)) != len(genes) or not all(map(ENSEMBL.fullmatch, genes)):
        raise ValueError("Expected the 6,638 unique Ensembl IDs from Step 30")

    record = json.loads(CATALOG.read_text(encoding="utf-8"))
    if record["doi"] != "10.25452/figshare.plus.29190726.v3":
        raise ValueError("Expected the pinned Figshare v3 record")
    filename = f"{args.cell_line}_filtered_dual_guide_cells.h5ad"
    asset = next(item for item in record["files"] if item["name"] == filename)
    url = f"https://api.figshare.com/v2/file/download/{asset['id']}"
    print(f"Probing var metadata for {args.cell_line}, Figshare file {asset['id']}", flush=True)
    with ranges.RangeReader(url, asset["size"]) as source:
        with h5py.File(source, "r") as handle:
            report, rows = inspect_var(handle["var"], genes)
        report.update(range_requests=source.requests, remote_bytes_read=source.bytes_read)
    report.update(
        cell_line=args.cell_line, figshare_doi=record["doi"], file_id=asset["id"],
        file_bytes=asset["size"], published_md5=asset["supplied_md5"], source_url=url,
        feature_axis_sha256=legacy.sha256(args.features),
        python_version=sys.version.split()[0], h5py_version=h5py.__version__,
        curl_version=subprocess.check_output(["curl", "--version"], text=True).splitlines()[0],
    )
    args.output_dir.mkdir(parents=True, exist_ok=False)
    with (args.output_dir / "feature_h5ad_mapping.tsv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_id", "h5ad_column_0based", "status"))
        writer.writerows(rows)
    summary = json.dumps(report, indent=2)
    (args.output_dir / "summary.json").write_text(summary + "\n", encoding="utf-8")
    print(summary, flush=True)


if __name__ == "__main__":
    main()
