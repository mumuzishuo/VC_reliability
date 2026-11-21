#!/usr/bin/env python3
"""Compare the audited 6,642 expression IDs to X-Atlas's feature dictionary.

This script reads only gene_metadata.parquet, not any cell expression or outcome.
"""

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


ENSEMBL_VERSION = re.compile(r"^(ENSG[0-9]+)\.[0-9]+$")
FIELDS = ("gene_id", "status", "n_external_rows", "external_gene_name",
          "external_gene_token_id", "versionless_external_ids")
OUTPUTS = ("feature_overlap.tsv", "summary.json", "parameters.json")


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def read_features(path, expected):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if "gene_id" not in (reader.fieldnames or []):
            raise ValueError("Feature table needs a gene_id column")
        ids = [row["gene_id"].strip() for row in reader]
    if len(ids) != expected or any(not value for value in ids) or len(set(ids)) != expected:
        raise ValueError("Expected %d distinct nonempty feature IDs" % expected)
    return ids


def versionless(value):
    return ENSEMBL_VERSION.sub(r"\1", value)


def build_reports(feature_ids, records):
    by_id = defaultdict(list)
    token_ids = defaultdict(set)
    missing_id_rows = 0
    missing_token_rows = 0
    for record in records:
        gene_id = str(record.get("ensembl_id") or "").strip()
        name = str(record.get("gene_name") or "").strip()
        token = record.get("gene_token_id")
        token = "" if token is None else str(token).strip()
        if not gene_id:
            missing_id_rows += 1
            continue
        if not token:
            missing_token_rows += 1
        else:
            token_ids[token].add(gene_id)
        by_id[gene_id].append((name, token))

    by_base = defaultdict(set)
    for gene_id in by_id:
        by_base[versionless(gene_id)].add(gene_id)

    rows = []
    status_counts = Counter()
    versioned_diagnostics = 0
    for gene_id in feature_ids:
        entries = by_id.get(gene_id, [])
        alternatives = sorted(by_base.get(versionless(gene_id), set()) - {gene_id})
        if len(entries) > 1:
            status = "ambiguous_external_id"
        elif len(entries) == 1:
            name, token = entries[0]
            if not token:
                status = "missing_external_token"
            elif len(token_ids[token]) > 1:
                status = "token_maps_to_multiple_ids"
            else:
                status = "exact_unique_token"
        elif alternatives:
            status = "no_exact_id_versionless_candidate"
            versioned_diagnostics += 1
        else:
            status = "no_external_id"
        status_counts[status] += 1
        rows.append((gene_id, status, len(entries),
                     entries[0][0] if len(entries) == 1 else "",
                     entries[0][1] if len(entries) == 1 else "",
                     "|".join(alternatives) if not entries else ""))

    summary = {
        "status": "FEATURE_DICTIONARY_ONLY",
        "matching_rule": "exact Ensembl ID; unique metadata row and unique token required for usable mapping; no automatic version stripping",
        "panel_expression_feature_count": len(feature_ids),
        "external_metadata_rows": len(records),
        "external_nonempty_unique_ids": len(by_id),
        "external_missing_id_rows": missing_id_rows,
        "external_missing_token_rows_with_id": missing_token_rows,
        "external_duplicated_exact_ids": sum(len(entries) > 1 for entries in by_id.values()),
        "external_versioned_ensembl_ids": sum(ENSEMBL_VERSION.fullmatch(gene_id) is not None for gene_id in by_id),
        "external_token_ids_used_by_multiple_gene_ids": sum(len(ids) > 1 for ids in token_ids.values()),
        "panel_exact_id_present": sum(len(by_id.get(gene_id, [])) > 0 for gene_id in feature_ids),
        "panel_exact_unique_token_mapping": status_counts["exact_unique_token"],
        "panel_unmatched_with_versionless_candidate": versioned_diagnostics,
        "panel_status_counts": dict(sorted(status_counts.items())),
        "unresolved": ["expression detection in each cell line", "observed guide and perturbation cell counts",
                       "controls per GEM/sample", "biological replicates", "external prediction performance"],
    }
    return rows, summary


def read_external(path):
    try:
        import pyarrow
        import pyarrow.parquet as pq
    except ImportError as error:
        raise RuntimeError("pyarrow is required to read the verified Parquet file") from error
    parquet = pq.ParquetFile(path)
    required = {"ensembl_id", "gene_name", "gene_token_id"}
    available = set(parquet.schema.names)
    if not required.issubset(available):
        raise ValueError("Missing expected Parquet fields: %r" % sorted(required - available))
    records = parquet.read(columns=sorted(required)).to_pylist()
    return records, pyarrow.__version__, parquet.metadata.num_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-metadata", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--expected-size", type=int, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--expected-features", type=int, default=6642)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    observed_hash = digest(args.gene_metadata)
    if observed_hash != args.expected_sha256.lower() or args.gene_metadata.stat().st_size != args.expected_size:
        raise ValueError("External metadata checksum or size differs from audited input")
    features = read_features(args.features, args.expected_features)
    records, arrow_version, parquet_rows = read_external(args.gene_metadata)
    if len(records) != parquet_rows:
        raise ValueError("Parquet row count differs from decoded records")
    rows, summary = build_reports(features, records)
    parameters = {
        "source": "https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/blob/main/metadata/gene_metadata.parquet",
        "gene_metadata_file": args.gene_metadata.name,
        "gene_metadata_sha256": observed_hash,
        "gene_metadata_size_bytes": args.expected_size,
        "features_file": args.features.name,
        "features_sha256": digest(args.features),
        "expected_feature_count": args.expected_features,
        "python_version": sys.version.split()[0],
        "pyarrow_version": arrow_version,
    }
    out = args.output_dir
    if out.exists():
        raise FileExistsError("Versioned output directory already exists: %s" % out)
    out.mkdir(parents=True, exist_ok=False)
    with (out / "feature_overlap.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        writer.writerows(rows)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (out / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")
    with (out / "SHA256SUMS").open("w", encoding="utf-8") as stream:
        for name in OUTPUTS:
            stream.write("%s  %s\n" % (digest(out / name), name))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
