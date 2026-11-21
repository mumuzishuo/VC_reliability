#!/usr/bin/env python3
"""Compare Jiang Supplementary Table 1f *design* labels with the frozen panel.

This deliberately reads only ST1f. It does not inspect expression, DE results,
or analysis sheets in the same publisher workbook.
"""

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl


SHEET = "ST1f_all_genes_by_pathway"
EXPECTED_PATHWAY_SETS = (
    ("IFNG", "IFNB", "TGFB", "TNFA", "INS"),
    ("IFNG", "IFNB", "TGFB1", "TNFA", "INS"),
)
OUTPUTS = ("design_overlap.tsv", "summary.json", "parameters.json")


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def read_panel(path, expected_count):
    labels = {}
    ids = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"gene", "gene_id", "in_main_panel"}.issubset(reader.fieldnames or []):
            raise ValueError("Task panel lacks gene, gene_id or in_main_panel")
        for row in reader:
            flag = row["in_main_panel"]
            if flag not in ("0", "1"):
                raise ValueError("Invalid in_main_panel flag")
            if flag == "0":
                continue
            gene, gene_id = row["gene"].strip(), row["gene_id"].strip()
            if not gene or not gene_id or gene in labels or gene_id in ids:
                raise ValueError("Missing or duplicate main-panel gene or gene_id")
            labels[gene] = gene_id
            ids.add(gene_id)
    if len(labels) != expected_count:
        raise ValueError("Expected %d main tasks, found %d" % (expected_count, len(labels)))
    return labels


def read_design(path):
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        if SHEET not in book.sheetnames:
            raise ValueError("Missing design sheet " + SHEET)
        sheet = book[SHEET]
        if sheet.max_row is None or sheet.max_row > 10000:
            raise ValueError("Unexpected ST1f dimensions; inspect metadata before parsing")
        header = next(sheet.iter_rows(min_row=3, max_row=3, max_col=2, values_only=True))
        if tuple(str(v).strip() for v in header) != ("Target Gene", "Pathway"):
            raise ValueError("Unexpected ST1f header: %r" % (header,))
        pairs = []
        seen = set()
        duplicates = 0
        for row in sheet.iter_rows(min_row=4, max_row=sheet.max_row, max_col=2, values_only=True):
            if all(v is None or not str(v).strip() for v in row):
                continue
            if not all(isinstance(v, str) and v.strip() for v in row):
                raise ValueError("Incomplete ST1f target/pathway row: %r" % (row,))
            target, pathway = (v.strip() for v in row)
            key = (pathway, target)
            if key in seen:
                duplicates += 1
            else:
                seen.add(key)
                pairs.append(key)
        observed_pathways = {pathway for pathway, _ in pairs}
        pathway_order = next(
            (labels for labels in EXPECTED_PATHWAY_SETS if observed_pathways == set(labels)),
            None,
        )
        if pathway_order is None:
            raise ValueError("Unexpected ST1f pathway set: %r" % sorted(observed_pathways))
        return sorted(pairs), duplicates, sheet.max_row, pathway_order
    finally:
        book.close()


def write_outputs(outdir, parameters, pairs, duplicates, reported_rows, pathway_order, panel):
    if any((outdir / name).exists() for name in (*OUTPUTS, "SHA256SUMS")):
        raise FileExistsError("Versioned output already exists; refusing overwrite")
    outdir.mkdir(parents=True, exist_ok=True)
    grouped = defaultdict(set)
    for pathway, target in pairs:
        grouped[pathway].add(target)
    all_targets = set.union(*grouped.values())
    overlapping = all_targets & panel.keys()
    counts = Counter(target for _, target in pairs)
    summary = {
        "status": "DESIGN_LABELS_ONLY",
        "source_sheet": SHEET,
        "pathways_as_reported_in_st1f": list(pathway_order),
        "label_matching": "exact, whitespace stripped; no aliases, case conversion or inferred IDs",
        "main_panel_size": len(panel),
        "st1f_reported_rows": reported_rows,
        "design_unique_pathway_target_pairs": len(pairs),
        "duplicate_pathway_target_rows": duplicates,
        "design_unique_target_symbols": len(all_targets),
        "targets_in_multiple_pathways": sum(n > 1 for n in counts.values()),
        "main_panel_overlap_unique_symbols": len(overlapping),
        "design_unique_symbols_without_exact_panel_match": len(all_targets - panel.keys()),
        "fraction_of_main_panel_covered_by_design": len(overlapping) / len(panel),
        "per_pathway": {
            p: {
                "design_target_symbols": len(grouped[p]),
                "main_panel_overlap_symbols": len(grouped[p] & panel.keys()),
                "without_exact_panel_match": len(grouped[p] - panel.keys()),
            }
            for p in pathway_order
        },
        "unresolved": [
            "external observed perturbation cell counts", "external guide recovery",
            "batch-matched controls", "external expression Ensembl ID overlap",
            "external response or prediction performance",
        ],
    }
    with (outdir / OUTPUTS[0]).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("pathway", "design_target_symbol", "panel_gene_id"))
        for pathway, target in pairs:
            if target in panel:
                writer.writerow((pathway, target, panel[target]))
    (outdir / OUTPUTS[1]).write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (outdir / OUTPUTS[2]).write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")
    with (outdir / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for name in OUTPUTS:
            handle.write("%s  %s\n" % (digest(outdir / name), name))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--expected-workbook-sha256", required=True)
    parser.add_argument("--task-panel", required=True, type=Path)
    parser.add_argument("--expected-tasks", type=int, default=1264)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    observed_hash = digest(args.workbook)
    if observed_hash != args.expected_workbook_sha256.lower():
        raise ValueError("Supplementary workbook SHA256 mismatch")
    panel = read_panel(args.task_panel, args.expected_tasks)
    pairs, duplicates, reported_rows, pathway_order = read_design(args.workbook)
    parameters = {
        "source_article": "https://doi.org/10.1038/s41556-025-01622-z",
        "workbook_name": args.workbook.name,
        "workbook_size_bytes": args.workbook.stat().st_size,
        "workbook_sha256": observed_hash,
        "task_panel_name": args.task_panel.name,
        "task_panel_sha256": digest(args.task_panel),
        "expected_tasks": args.expected_tasks,
        "python_version": sys.version.split()[0],
        "openpyxl_version": openpyxl.__version__,
        "read_sheets": [SHEET],
    }
    write_outputs(args.output_dir, parameters, pairs, duplicates, reported_rows, pathway_order, panel)


if __name__ == "__main__":
    main()
