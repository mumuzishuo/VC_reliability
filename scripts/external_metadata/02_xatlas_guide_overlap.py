#!/usr/bin/env python3
"""Audit X-Atlas guide *design* and compare target IDs with the 1,264-task panel.

Does not read any external single-cell expression or observed guide calls.
"""

import argparse
import csv
import hashlib
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED = {"target_gene_name", "target_gene_id", "id_a", "id_b", "id_ab"}
ENSEMBL = re.compile(r"^ENSG[0-9]{11}(?:\.[0-9]+)?$")
OUTPUTS = ("summary.json", "overlap_by_gene_id.tsv", "symbol_only_diagnostic.tsv", "parameters.json")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def md5_digest(path):
    h = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def panel_ids(path, expected_tasks):
    by_id = {}
    by_symbol = {}
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"gene", "gene_id", "in_main_panel"}.issubset(reader.fieldnames or []):
            raise ValueError("Task panel missing gene, gene_id or in_main_panel")
        for row in reader:
            if row["in_main_panel"] not in ("0", "1"):
                raise ValueError("Invalid task panel flag")
            if row["in_main_panel"] == "0":
                continue
            symbol, gene_id = row["gene"].strip(), row["gene_id"].strip()
            if not symbol or not gene_id or symbol in by_symbol or gene_id in by_id:
                raise ValueError("Ambiguous or missing main-panel target")
            by_id[gene_id] = symbol
            by_symbol[symbol] = gene_id
    if len(by_id) != expected_tasks:
        raise ValueError("Expected %d main-panel tasks; found %d" % (expected_tasks, len(by_id)))
    return by_id, by_symbol


def audit_library(path, expected_rows):
    ids_to_symbols = defaultdict(set)
    symbols_to_ids = defaultdict(set)
    pairs_by_id = defaultdict(set)
    all_ids = set()
    all_symbols = set()
    pair_targets = {}
    status = Counter()
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError("Guide library missing expected columns: %r" % (reader.fieldnames,))
        for row in reader:
            status["rows"] += 1
            symbol = (row["target_gene_name"] or "").strip()
            gene_id = (row["target_gene_id"] or "").strip()
            a, b, pair = ((row[field] or "").strip() for field in ("id_a", "id_b", "id_ab"))
            if not symbol:
                status["missing_symbol_rows"] += 1
            else:
                all_symbols.add(symbol)
            if not gene_id:
                status["missing_gene_id_rows"] += 1
            else:
                all_ids.add(gene_id)
            if not a or not b or not pair:
                status["missing_guide_pair_rows"] += 1
            if gene_id and not ENSEMBL.fullmatch(gene_id):
                status["non_ensembl_id_rows"] += 1
            if "." in gene_id and ENSEMBL.fullmatch(gene_id):
                status["versioned_ensembl_id_rows"] += 1
            if symbol and gene_id:
                ids_to_symbols[gene_id].add(symbol)
                symbols_to_ids[symbol].add(gene_id)
            if gene_id and pair:
                pairs_by_id[gene_id].add(pair)
                if pair in pair_targets:
                    status["duplicate_guide_pair_rows"] += 1
                    if pair_targets[pair] != gene_id:
                        raise ValueError("Guide pair maps to multiple target gene IDs: %r" % pair)
                else:
                    pair_targets[pair] = gene_id
    if status["rows"] != expected_rows:
        raise ValueError("Expected %d library rows; found %d" % (expected_rows, status["rows"]))
    return ids_to_symbols, symbols_to_ids, pairs_by_id, all_ids, all_symbols, status


def write_table(path, fields, data):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(data)


def report(out, params, panel_by_id, panel_by_symbol, ids_to_symbols, symbols_to_ids,
           pairs_by_id, all_ids, all_symbols, status):
    if any((out / name).exists() for name in (*OUTPUTS, "SHA256SUMS")):
        raise FileExistsError("Versioned output already exists; refusing overwrite")
    matched_ids = set(panel_by_id) & ids_to_symbols.keys()
    matched_symbols = set(panel_by_symbol) & symbols_to_ids.keys()
    symbol_only = sorted(
        s for s in matched_symbols if panel_by_symbol[s] not in symbols_to_ids[s]
    )
    id_rows = []
    inconsistent = []
    for gene_id in sorted(matched_ids):
        panel_symbol = panel_by_id[gene_id]
        source_symbols = sorted(ids_to_symbols[gene_id])
        symbol_agrees = panel_symbol in source_symbols
        if not symbol_agrees:
            inconsistent.append(gene_id)
        id_rows.append((panel_symbol, gene_id, "|".join(source_symbols),
                        len(pairs_by_id[gene_id]), int(symbol_agrees)))
    only_rows = [
        (s, panel_by_symbol[s], "|".join(sorted(symbols_to_ids[s])))
        for s in symbol_only
    ]
    counts = [len(pairs_by_id[gene_id]) for gene_id in matched_ids]
    summary = {
        "status": "GUIDE_DESIGN_ONLY",
        "matching_rule": "exact Ensembl gene_id; exact symbol is a separate diagnostic; no alias or version normalization",
        "main_panel_tasks": len(panel_by_id),
        "guide_library_rows": status["rows"],
        "unique_nonempty_target_gene_ids": len(all_ids),
        "unique_ensembl_target_gene_ids": sum(bool(ENSEMBL.fullmatch(g)) for g in all_ids),
        "unique_nonempty_target_symbols": len(all_symbols),
        "unique_guide_pair_ids": len(set().union(*pairs_by_id.values())) if pairs_by_id else 0,
        "id_overlap_with_main_panel": len(matched_ids),
        "symbol_overlap_with_main_panel": len(matched_symbols),
        "id_matches_with_exact_symbol_agreement": len(matched_ids) - len(inconsistent),
        "id_matches_without_exact_symbol_agreement": len(inconsistent),
        "symbol_matches_without_panel_id": len(symbol_only),
        "designed_guide_pairs_per_id_matched_target": {
            "min": min(counts) if counts else None,
            "median": statistics.median(counts) if counts else None,
            "max": max(counts) if counts else None,
        },
        "diagnostic_row_counts": dict(sorted(status.items())),
        "unresolved": ["observed cells per perturbation", "guide capture per cell",
                       "GEM-matched control counts", "independent biological replicates",
                       "external response-gene Ensembl overlap", "prediction performance"],
    }
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / "overlap_by_gene_id.tsv",
                ("panel_symbol", "gene_id", "library_symbols", "designed_guide_pairs", "symbol_agrees"),
                id_rows)
    write_table(out / "symbol_only_diagnostic.tsv",
                ("panel_symbol", "panel_gene_id", "library_gene_ids"), only_rows)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "parameters.json").write_text(json.dumps(params, indent=2) + "\n", encoding="utf-8")
    with (out / "SHA256SUMS").open("w", encoding="utf-8") as stream:
        for name in OUTPUTS:
            stream.write("%s  %s\n" % (digest(out / name), name))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--guide-library", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--expected-md5", default="bf1181127b314affe95e5b5ca97bdbf0")
    parser.add_argument("--expected-size", type=int, default=3648784)
    parser.add_argument("--expected-rows", type=int, default=20890)
    parser.add_argument("--task-panel", type=Path, required=True)
    parser.add_argument("--expected-tasks", type=int, default=1264)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    observed_hash = digest(args.guide_library)
    if observed_hash != args.expected_sha256.lower() or args.guide_library.stat().st_size != args.expected_size:
        raise ValueError("Guide library SHA256 or size differs from audited file")
    if md5_digest(args.guide_library) != args.expected_md5.lower():
        raise ValueError("Guide library MD5 differs from Figshare record")
    panel_by_id, panel_by_symbol = panel_ids(args.task_panel, args.expected_tasks)
    ids_to_symbols, symbols_to_ids, pairs_by_id, all_ids, all_symbols, status = audit_library(
        args.guide_library, args.expected_rows)
    params = {
        "figshare_record_id": 29190726, "figshare_version": 3,
        "figshare_file_id": 57368587, "guide_library_name": args.guide_library.name,
        "guide_library_sha256": observed_hash, "guide_library_size_bytes": args.expected_size,
        "figshare_supplied_md5": args.expected_md5.lower(),
        "task_panel_name": args.task_panel.name, "task_panel_sha256": digest(args.task_panel),
        "expected_tasks": args.expected_tasks, "expected_guide_rows": args.expected_rows,
        "python_version": sys.version.split()[0],
    }
    report(args.output_dir, params, panel_by_id, panel_by_symbol,
           ids_to_symbols, symbols_to_ids, pairs_by_id, all_ids, all_symbols, status)


if __name__ == "__main__":
    main()
