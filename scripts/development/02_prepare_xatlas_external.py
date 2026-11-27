#!/usr/bin/env python3
"""Project Step 29 predictions onto shared X-Atlas metadata axes."""

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "legacy_v01"))
import vc_baseline as baseline
import vc_preprocess as preprocess


RUN = Path("results/v01_four_source_calibration/job_3239681")
RUN_MANIFEST = "manifests/v01_four_source_calibration_3239681.sha256"
SUPPORT = "qc/v01_xatlas_sample_support/task_sample_support.tsv"
FEATURES = "qc/v01_xatlas_expression_feature_overlap/feature_overlap.tsv"
TARGETS = ("HCT116", "HEK293T")
SOURCES = "K562,RPE1,Jurkat,HepG2"
MIN_CELLS = 30
MIN_CONTROLS = 30
COVERAGE = 0.8
PRIMARY_RISK = "source_held_out_excess_mse"
RISK_SCORES = (PRIMARY_RISK, "source_heterogeneity", "low_source_support", "predicted_rms")
STATUS = "XATLAS_PREDICTION_CANDIDATE"
require = baseline.require


def verify(manifest, *files):
    """Verify named inputs against an earlier stage's recorded SHA256."""
    hashes = {filename: digest for digest, filename in
              (line.split(maxsplit=1) for line in Path(manifest).read_text().splitlines())}
    for file in files:
        filename = Path(file).as_posix()
        require(hashes.get(filename) == baseline.digest(file), "SHA256 mismatch: " + filename)
    return {str(Path(manifest)): baseline.digest(manifest),
            **{Path(file).as_posix(): hashes[Path(file).as_posix()] for file in files}}


def read_table(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def task_panel(path):
    by_line = {line: {} for line in TARGETS}
    for row in read_table(path):
        line, gene = row["cell_line"], row["gene_id"]
        require(line in by_line and gene and gene not in by_line[line],
                "Unexpected or repeated task row: " + line + "/" + gene)
        total = int(row["cells"])
        low = int(row["cells_in_gems_lt30_controls"])
        require(0 <= low <= total, "Invalid cell support: " + line + "/" + gene)
        by_line[line][gene] = (row["panel_symbol"], total - low)

    hct, hek = (by_line[line] for line in TARGETS)
    require(hct.keys() == hek.keys(), "Task IDs differ between guide-count tables")
    rows = []
    for gene in sorted(hct):
        symbol, n_hct = hct[gene]
        require(symbol and symbol == hek[gene][0], "Task symbol differs: " + gene)
        n_hek = hek[gene][1]
        rows.append(dict(target_gene_id=gene, source_perturbation=symbol,
                         HCT116_cells_ge30_controls=n_hct, HEK293T_cells_ge30_controls=n_hek,
                         include_in_candidate=int(min(n_hct, n_hek) >= MIN_CELLS)))
    require(len({row["source_perturbation"] for row in rows}) == len(rows),
            "Candidate identifiers do not map one to one to internal perturbations")
    return rows


def feature_panel(path):
    ids, tokens, rows = set(), set(), []
    for row in read_table(path):
        gene = row["gene_id"]
        require(gene and gene not in ids, "Empty or repeated expression gene ID: " + gene)
        ids.add(gene)
        if row["status"] == "exact_unique_token":
            token = int(row["external_gene_token_id"])
            require(token >= 0 and token not in tokens, "Invalid or repeated feature token")
            tokens.add(token)
            rows.append(dict(gene_id=gene, external_gene_token_id=token))
    return ids, sorted(rows, key=lambda row: row["gene_id"])


def axis_order(dataset, full_labels, selected_labels):
    labels = baseline.strings(dataset)
    index = {value: i for i, value in enumerate(labels)}
    require(len(index) == len(labels) and index.keys() == full_labels,
            "The archived prediction axis differs from the audited metadata")
    return np.array([index[value] for value in selected_labels], dtype=int)


def project(path, summary, task_rows, full_gene_ids, features):
    tasks = [row for row in task_rows if row["include_in_candidate"]]
    task_ids = [row["target_gene_id"] for row in tasks]
    symbols = [row["source_perturbation"] for row in tasks]
    genes = [row["gene_id"] for row in features]
    n_tasks, n_genes = summary["n_tasks"], summary["n_genes"]
    require(len(task_rows) == n_tasks and len(full_gene_ids) == n_genes,
            "Metadata axis sizes disagree with Step 29")
    require(len(tasks) > 1 and genes, "Candidate panel is empty or too small for retrieval")

    scales = summary["four_source_calibration"]
    alpha, beta = scales["paired"]["scale"], scales["common"]["scale"]
    with h5py.File(path, "r") as source:
        require(source.attrs["status"] == "INTERNAL_DEVELOPMENT_ONLY"
                and source.attrs["sources"] == SOURCES, "Unexpected source prediction file")
        require(source.attrs["alpha_paired"] == alpha and source.attrs["alpha_common"] == beta
                and 0 <= alpha <= 1 and 0 <= beta <= 1, "Scale mismatch with Step 29")
        task_ix = axis_order(source["perturbation"],
                             {r["source_perturbation"] for r in task_rows}, symbols)
        gene_ix = axis_order(source["gene_id"], full_gene_ids, genes)
        require(source["paired_raw"].shape == (n_tasks, n_genes)
                and source["common_raw"].shape == (n_genes,), "Unexpected prediction shape")
        require(all(source[name].shape == (n_tasks,) for name in (*RISK_SCORES, "tie_key")),
                "Unexpected risk or tie-key shape")
        paired = source["paired_raw"][:][np.ix_(task_ix, gene_ix)]
        common = source["common_raw"][:][gene_ix]
        scores = {key: source[key][:][task_ix] for key in RISK_SCORES}
        ties = source["tie_key"][:][task_ix]

    require(all(np.isfinite(a).all() for a in (paired, common, *scores.values())),
            "Nonfinite predictions or risk")
    require(np.unique(ties).size == len(tasks), "Duplicate tie keys")
    keep = np.zeros(len(tasks), dtype=np.uint8)
    k = math.ceil(COVERAGE * len(tasks))
    keep[baseline.risk_order(scores[PRIMARY_RISK], ties)[:k]] = 1
    arrays = dict(target_gene_id=task_ids, source_perturbation=symbols, gene_id=genes,
                  paired_raw=paired, common_raw=common, tie_key=ties, keep_primary_80=keep, **scores)
    return arrays, alpha, beta, k


def write_json(path, data):
    with path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(data, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    out = parser.parse_args().output_dir
    require(not out.exists(), "Use a new output directory: " + str(out))

    run_summary = RUN / "summary.json"
    run_hashes = RUN / "SHA256SUMS"
    source_predictions = RUN / "predictions.h5"
    input_hashes = {
        **verify(RUN_MANIFEST, run_summary, run_hashes),
        **verify(run_hashes, source_predictions),
        **verify("manifests/v01_xatlas_sample_support_public.sha256", SUPPORT),
        **verify("manifests/v01_xatlas_expression_feature_review.sha256", FEATURES),
    }
    tasks = task_panel(SUPPORT)
    full_gene_ids, features = feature_panel(FEATURES)
    summary = json.loads(run_summary.read_text(encoding="utf-8"))
    require(summary["status"] == "INTERNAL_DEVELOPMENT_ONLY", "Unexpected source run")
    arrays, alpha, beta, k = project(source_predictions, summary, tasks, full_gene_ids, features)

    out.mkdir(parents=True, exist_ok=False)
    baseline.write_table(out / "task_eligibility.tsv", tasks)
    baseline.write_table(out / "feature_axis.tsv", features)
    preprocess.write_h5(out / "predictions.h5", arrays,
                        dict(status=STATUS, sources=SOURCES, alpha_paired=alpha, alpha_common=beta))
    write_json(out / "parameters.json", dict(
        status=STATUS, created_utc=datetime.now(timezone.utc).isoformat(),
        primary_target=TARGETS[0], second_target_same_study=TARGETS[1],
        min_cells_per_context=MIN_CELLS, min_controls_per_gem=MIN_CONTROLS,
        fit_panel="full internal Step 29 panel", nominal_coverage=COVERAGE,
        retention_rule="ceil(coverage * N), ascending risk then archived tie key",
        primary_risk=PRIMARY_RISK, risk_scores=list(RISK_SCORES),
        failure_event="paired_calibrated MSE > zero_response MSE",
        input_sha256=input_hashes,
        code_sha256={Path(path).name: baseline.digest(path)
                     for path in (__file__, baseline.__file__, preprocess.__file__)},
        versions=dict(python=sys.version.split()[0], numpy=np.__version__, h5py=h5py.__version__),
    ))
    result = dict(status=STATUS, n_tasks=len(arrays["target_gene_id"]), n_genes=len(features),
                  atlas_axis_confirmed=False, alpha_paired=alpha, alpha_common=beta,
                  n_retained_primary_80=k, actual_coverage=k / len(arrays["target_gene_id"]))
    write_json(out / "summary.json", result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
