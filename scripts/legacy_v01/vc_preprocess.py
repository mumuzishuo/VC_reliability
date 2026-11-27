#!/usr/bin/env python3
"""VC reliability pilot: prepare -> smoke -> run (Python 3.9+).

Dependencies: numpy, scipy, h5py. No anndata/scanpy installation is needed.
Upload this file to $HOME/VC_reliability/vc_preprocess.py, then run:
  python vc_preprocess.py prepare
  python vc_preprocess.py smoke
Full scans belong on an allocated compute node, using the site's scheduler:
  python vc_preprocess.py run --environment K562
  python vc_preprocess.py run --environment RPE1
  python vc_preprocess.py run --environment Jurkat
  python vc_preprocess.py run --environment HepG2

Working protocol (v0.1, developmental; not a claim of optimal preprocessing):
* Respect the authors' existing cell QC; require valid, positive UMI_count.
* Use exact shared var/gene_id identifiers and an identical sorted column order.
* Use observed obs/gene labels, excluding non-targeting, with >=30 cells per
  environment by default. This is evaluation cohort eligibility, not a risk
  feature; target perturbation counts must not be fed into a deployed risk model.
* y_ig = ln(1 + 10000 * raw_count_ig / obs.UMI_count_i). The denominator is
  the original metadata total, NOT the sum over the 6642 selected features.
  These are relative-expression changes, not absolute total-RNA changes.
* delta_pg = mean(y_ig - mean_control_in_the_same_gem_group_g).
  This matches the control mixture to each perturbation's cell counts by group.
* Retain groups with >=1 control by default. Record the fraction of each
  perturbation's cells in groups with <30 controls; do not call groups biological
  replicates. --min-controls 30 provides a separately named sensitivity run.
* HepG2 outcomes go to evaluation_only; source outcomes go to sources. The
  separation is organizational, not a security barrier. Do not fit/tune models
  from evaluation_only. HepG2 controls are also saved separately as allowed input.
* Output files are ordinary HDF5 (.h5), not AnnData (.h5ad). Log-differences are
  NOT log2 fold changes. No DEG selection, batch integration or predictive model
  is performed here. Cell-level uncertainty and independent validation remain.

The raw inputs are opened read-only. Existing outputs are never overwritten.
K562/RPE1 checksums were verified in the supplied HPC transcript. HepG2/Jurkat
external checksum verification remains outstanding; a full scan validates all
matrix values but does not substitute for an external checksum.
"""

import argparse
import csv
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

FILES = {
    "K562": "K562_essential_raw_singlecell_01.h5ad",
    "RPE1": "rpe1_raw_singlecell_01.h5ad",
    "Jurkat": "GSE264667_jurkat_raw_singlecell_01.h5ad",
    "HepG2": "GSE264667_hepg2_raw_singlecell_01.h5ad",
}
CONTROL = "non-targeting"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def column(group, name):
    item = group[name]
    if isinstance(item, h5py.Group):
        values, categories = item["codes"][:], item["categories"][:]
    else:
        values = item[:]
        categories = (group["__categories"][name][:]
                      if "__categories" in group and name in group["__categories"]
                      else None)
    if categories is not None:
        codes = np.asarray(values)
        require(np.issubdtype(codes.dtype, np.integer), "Non-integer category codes")
        require(bool(((codes >= -1) & (codes < len(categories))).all()),
                "Invalid category codes")
        values = [categories[int(i)] if i >= 0 else "<missing>" for i in codes]
    return np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v)
                       for v in values])


def metadata(path):
    with h5py.File(path, "r") as f:
        require(isinstance(f["X"], h5py.Dataset), "This script requires dense X")
        genes = column(f["obs"], "gene")
        batches = column(f["obs"], "gem_group")
        ids = column(f["var"], "gene_id")
        umi = column(f["obs"], "UMI_count").astype(np.float64)
        require(f["X"].shape == (len(genes), len(ids)), "Metadata / X shape mismatch")
    require(len(batches) == len(genes) == len(umi), "Cell metadata length mismatch")
    for name, values in [("perturbation", genes), ("gem_group", batches), ("gene_id", ids)]:
        require(not np.isin(values, ["", "nan", "None", "<missing>"]).any(),
                "Missing " + name)
    require(len(set(ids)) == len(ids), "Duplicate expression gene IDs")
    require(bool((np.isfinite(umi) & (umi > 0)).all()), "Invalid UMI_count metadata")
    return genes, batches, ids, umi


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(args, work):
    path = work / "manifest.json"
    require(not path.exists(), "Manifest exists. Use smoke/run, or a new --output directory.")
    all_meta, counts, shared_ids, inputs, rows = {}, {}, None, {}, []
    for env, filename in FILES.items():
        file = args.root / "raw" / filename
        genes, batches, ids, umi = metadata(file)
        all_meta[env] = (genes, batches)
        controls = Counter(batches[genes == CONTROL])
        valid_batches = {b for b, n in controls.items() if n >= args.min_controls}
        valid = np.isin(batches, list(valid_batches))
        c = Counter(genes[valid])
        c.pop(CONTROL, None)
        counts[env] = c
        shared_ids = set(ids) if shared_ids is None else shared_ids.intersection(ids)
        stat = file.stat()
        inputs[env] = {"path": str(file.resolve()), "size_bytes": stat.st_size,
                       "mtime_ns": stat.st_mtime_ns, "n_cells": len(genes), "n_genes": len(ids)}
        rows.extend((env, b, n, controls[b], int(b in valid_batches))
                    for b, n in sorted(Counter(batches).items()))
        print(f"{env}: {len(genes):,} cells; {len(ids):,} genes; "
              f"{len(valid_batches)} usable GEM groups", flush=True)
    shared = set.intersection(*(set(c) for c in counts.values()))
    perturbations = sorted(g for g in shared
                           if all(c[g] >= args.min_cells for c in counts.values()))
    require(bool(shared_ids) and bool(perturbations), "Empty analysis cohort")
    manifest = {"protocol": "VC-pilot-v0.1", "script_sha256": sha256(Path(__file__)),
                "target": "HepG2", "control_label": CONTROL,
                "min_cells": args.min_cells, "min_controls": args.min_controls,
                "normalization": "ln(1 + 10000 * X / obs.UMI_count)",
                "gene_ids": sorted(shared_ids), "perturbations": perturbations,
                "inputs": inputs,
                "checksum_status": {"K562": "user-reported MD5 pass",
                                    "RPE1": "user-reported MD5 pass",
                                    "HepG2": "external checksum not yet verified",
                                    "Jurkat": "external checksum not yet verified"}}
    work.mkdir(parents=True, exist_ok=True)
    with (work / "batch_inventory.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["environment", "gem_group", "n_cells", "n_controls", "included"])
        w.writerows(rows)
    with (work / "perturbation_cell_counts.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["perturbation"] + list(FILES))
        w.writerows([g] + [counts[e][g] for e in FILES] for g in perturbations)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"PREPARE_OK: {len(perturbations)} perturbations x {len(shared_ids)} expression genes")
    print("Manifest:", path)


def load_manifest(work):
    manifest = json.loads((work / "manifest.json").read_text(encoding="utf-8"))
    require(manifest["script_sha256"] == sha256(Path(__file__)),
            "Script changed since prepare. Use a new --output directory.")
    for env, item in manifest["inputs"].items():
        stat = Path(item["path"]).stat()
        require(stat.st_size == item["size_bytes"] and stat.st_mtime_ns == item["mtime_ns"],
                "Input changed since prepare: " + env)
    return manifest


def validate_block(raw, umi):
    require(bool(np.isfinite(raw).all()) and bool((raw >= 0).all()),
            "Nonfinite or negative matrix values")
    require(bool(np.equal(raw, np.floor(raw)).all()), "Non-integer raw counts")
    ratios = raw.sum(axis=1, dtype=np.float64) / umi
    require(bool((ratios <= 1.0001).all()),
            "Matrix row sum exceeds UMI_count; denominator needs investigation")
    return ratios


def normalized(raw, umi, indices):
    y = raw[:, indices].astype(np.float32, copy=True)
    y *= (10000.0 / umi).astype(np.float32)[:, None]
    np.log1p(y, out=y)
    require(bool(np.isfinite(y).all()), "Nonfinite normalized values")
    return y


def accumulate(sums, codes, y):
    # Use a sparse membership matrix; never allocate cells x all perturbations.
    from scipy.sparse import csr_matrix
    unique, inverse = np.unique(codes, return_inverse=True)
    membership = csr_matrix((np.ones(len(codes), dtype=np.float64),
                             (inverse, np.arange(len(codes)))),
                            shape=(len(unique), len(codes)))
    sums[unique] += membership @ y


def smoke(manifest):
    import scipy
    print("Versions:", sys.version.split()[0], "numpy", np.__version__,
          "scipy", scipy.__version__, "h5py", h5py.__version__)
    for env, info in manifest["inputs"].items():
        genes, batches, ids, umi = metadata(info["path"])
        rows = np.unique(np.linspace(0, len(genes) - 1, min(64, len(genes)), dtype=int))
        positions = {g: i for i, g in enumerate(ids)}
        indices = np.asarray([positions[g] for g in manifest["gene_ids"]])
        with h5py.File(info["path"], "r") as f:
            raw = f["X"][rows, :]
        ratio = validate_block(raw, umi[rows])
        y = normalized(raw, umi[rows], indices)
        codes = np.arange(len(rows)) % 3
        sums = np.zeros((3, y.shape[1]), dtype=np.float64)
        accumulate(sums, codes, y)
        expected = np.stack([y[codes == k].sum(axis=0, dtype=np.float64) for k in range(3)])
        require(bool(np.allclose(sums, expected, rtol=1e-10, atol=1e-10)),
                "Aggregation smoke test failed")
        print(f"SMOKE_OK {env}: {len(rows)} cells; "
              f"stored-count / UMI_count={ratio.min():.6f}..{ratio.max():.6f}", flush=True)


def write_h5(path, arrays, attributes):
    require(not path.exists(), "Output exists; refusing overwrite: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".partial.h5")
    require(not temp.exists(), "Partial output exists; inspect it before retry: " + str(temp))
    with h5py.File(temp, "x") as h:
        for name, values in arrays.items():
            arr = np.asarray(values)
            if arr.dtype.kind in "UO":
                h.create_dataset(name, data=arr.astype(object), dtype=h5py.string_dtype("utf-8"))
            else:
                h.create_dataset(name, data=arr, compression="gzip", compression_opts=1,
                                 shuffle=True)
        for name, value in attributes.items():
            h.attrs[name] = value
    # Publish atomically without replacing any pre-existing final output.
    os.link(temp, path)
    temp.unlink()


def run_environment(manifest, env, work, chunk_size):
    info = manifest["inputs"][env]
    role = "evaluation_only" if env == manifest["target"] else "sources"
    output = work / role / (env + "_responses.h5")
    require(not output.exists(), "Output already exists: " + str(output))
    genes, batches, ids, umi = metadata(info["path"])
    perturbations = manifest["perturbations"]
    p_count, g_count = len(perturbations), len(manifest["gene_ids"])
    control_counts = Counter(batches[genes == CONTROL])
    batch_names = sorted(b for b, n in control_counts.items() if n >= manifest["min_controls"])
    b_lookup = {b: i for i, b in enumerate(batch_names)}
    p_lookup = {p: i for i, p in enumerate(perturbations)}
    positions = {g: i for i, g in enumerate(ids)}
    indices = np.asarray([positions[g] for g in manifest["gene_ids"]])
    b_codes = np.asarray([b_lookup.get(b, -1) for b in batches], dtype=np.int32)
    p_codes = np.asarray([p_lookup.get(p, -1) for p in genes], dtype=np.int32)
    is_control = genes == CONTROL
    sums = np.zeros((p_count + len(batch_names), g_count), dtype=np.float64)
    n_control = np.zeros(len(batch_names), dtype=np.int64)
    n_by_batch = np.zeros((p_count, len(batch_names)), dtype=np.int64)
    ratio_min, ratio_max = np.inf, -np.inf
    with h5py.File(info["path"], "r") as f:
        for start in range(0, len(genes), chunk_size):
            stop = min(start + chunk_size, len(genes))
            raw = f["X"][start:stop, :]
            ratios = validate_block(raw, umi[start:stop])
            ratio_min, ratio_max = min(ratio_min, ratios.min()), max(ratio_max, ratios.max())
            bc, pc, ctrl = b_codes[start:stop], p_codes[start:stop], is_control[start:stop]
            keep = (bc >= 0) & (ctrl | (pc >= 0))
            if keep.any():
                y = normalized(raw[keep], umi[start:stop][keep], indices)
                bc, pc, ctrl = bc[keep], pc[keep], ctrl[keep]
                codes = np.where(ctrl, p_count + bc, pc)
                accumulate(sums, codes, y)
                np.add.at(n_control, bc[ctrl], 1)
                np.add.at(n_by_batch, (pc[~ctrl], bc[~ctrl]), 1)
            if start == 0 or stop == len(genes) or start // chunk_size % 25 == 0:
                print(f"{env}: {stop:,}/{len(genes):,} cells scanned", flush=True)
    n_cells = n_by_batch.sum(axis=1)
    require(bool((n_cells >= manifest["min_cells"]).all()), "Insufficient perturbation cells")
    require(bool((n_control >= manifest["min_controls"]).all()), "Insufficient controls")
    means = sums[:p_count] / n_cells[:, None]
    ctrl_means = sums[p_count:] / n_control[:, None]
    weights = n_by_batch / n_cells[:, None]
    matched = weights @ ctrl_means
    delta = means - matched
    low_fraction = n_by_batch[:, n_control < 30].sum(axis=1) / n_cells
    attrs = {"protocol": manifest["protocol"], "role": role, "environment": env,
             "normalization": manifest["normalization"],
             "manifest_sha256": sha256(work / "manifest.json"),
             "stored_count_fraction_min": float(ratio_min),
             "stored_count_fraction_max": float(ratio_max),
             "warning": "GEM groups are not established independent biological replicates"}
    common = {"gene_id": manifest["gene_ids"], "gem_group": batch_names,
              "control_mean_by_batch": ctrl_means.astype(np.float32),
              "control_cell_count_by_batch": n_control}
    write_h5(output, {**common, "perturbation": perturbations, "n_cells": n_cells,
                     "n_cells_by_batch": n_by_batch,
                     "low_control_fraction": low_fraction.astype(np.float32),
                     "perturbed_mean_log1p": means.astype(np.float32),
                     "matched_control_mean_log1p": matched.astype(np.float32),
                     "delta": delta.astype(np.float32)}, attrs)
    if env == manifest["target"]:
        write_h5(work / "target_controls" / "HepG2_controls.h5", common,
                 {**attrs, "role": "allowed_target_control_input"})
    print(f"RUN_OK {env}: delta shape={delta.shape}; {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=["prepare", "smoke", "run"])
    parser.add_argument("--root", type=Path, default=Path.home() / "VC_reliability")
    parser.add_argument("--output", type=Path, help="Default: ROOT/processed/v01")
    parser.add_argument("--environment", choices=list(FILES) + ["all"])
    parser.add_argument("--min-cells", type=int, default=30, help="prepare only")
    parser.add_argument("--min-controls", type=int, default=1, help="prepare only")
    parser.add_argument("--chunk-size", type=int, default=512, help="run only")
    args = parser.parse_args()
    require(args.min_cells > 0 and args.min_controls > 0 and args.chunk_size > 0,
            "Thresholds and chunk size must be positive")
    work = args.output or args.root / "processed" / "v01"
    if args.mode == "prepare":
        prepare(args, work)
    else:
        manifest = load_manifest(work)
        if args.mode == "smoke":
            smoke(manifest)
        else:
            require(args.environment is not None, "run requires --environment")
            for env in FILES if args.environment == "all" else [args.environment]:
                run_environment(manifest, env, work, args.chunk_size)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        sys.exit(1)
