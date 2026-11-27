#!/usr/bin/env python3
"""Internal VC reliability validation; Python 3.9+, numpy/scipy/h5py.

Keep next to the original vc_preprocess.py and vc_baseline.py.
run = predict -> evaluate -> split -> stability. Existing completed stages are
verified and reused. Originals are read-only; HepG2 response/raw files are never
opened. The inherited v01 cohort DID use all four environments for eligibility.
Repeated balanced half-splits are a sampling sensitivity analysis, NOT a
bootstrap confidence interval or independent biological replication.
"""
import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import scipy
from scipy.stats import rankdata
import vc_preprocess as pp
import vc_baseline as bb

ENVS = ("K562", "RPE1", "Jurkat")
RULES = ("predicted_rms", "low_source_support", "source_heterogeneity",
         "base_rank_mean", "extended_rank_mean")
SEED = 20261102
COVERAGE = 0.8
require = bb.require


def json_write(path, data):
    temp = path.with_name(path.name + ".partial")
    with temp.open("x", encoding="utf-8") as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write("\n")
    bb.commit(temp, path)


def code_hashes():
    return {p.name: bb.digest(p) for p in
            [Path(__file__), Path(pp.__file__), Path(bb.__file__)]}


def source_path(root, env):
    require(env in ENVS, "Only the three development environments are allowed")
    return root / "processed/v01/sources" / (env + "_responses.h5")


def load_source(root, env, manifest):
    p, g, d, n, low, info = bb.read_response(source_path(root, env), "sources", env)
    require(info["manifest_sha256"] == bb.digest(root / "processed/v01/manifest.json"),
            "Source manifest mismatch")
    require(info["normalization"] == manifest["normalization"], "Normalization mismatch")
    d, n, low = bb.align(manifest["perturbations"], manifest["gene_ids"], p, g, d, n, low)
    return d, n


def make_prediction(deltas, counts):
    require(len(deltas) == len(counts) == 2, "Exactly two sources per internal fold")
    prediction = (deltas[0] + deltas[1]) / 2.0
    # For two sources the population variance is exactly (a-b)^2 / 4.
    h = np.mean((deltas[0] - deltas[1]) ** 2, axis=1) / 4.0
    a = np.sqrt(np.mean(prediction ** 2, axis=1))
    s = 1.0 / np.minimum(counts[0], counts[1])
    require(np.isfinite(prediction).all() and (np.minimum(*counts) > 0).all(),
            "Invalid source response/support")
    # Average ranks for ties. Percentile convention is rank / N, fixed in advance.
    ra, rs, rh = (rankdata(v, method="average") / len(v) for v in (a, s, h))
    return prediction, dict(zip(RULES, [a, s, h, (ra + rs) / 2, (ra + rs + rh) / 3]))


def metric_rows(error, zero, features, tie, prefix):
    n = len(error)
    k = math.ceil(COVERAGE * n)
    rows = []
    for name in RULES:
        order = bb.risk_order(features[name], tie)
        ix = order[:k]
        retained = float(error[ix].mean())
        z = float(zero[ix].mean())
        rows.append({**prefix, "rule": name, "n_tasks": n, "retained_n": k,
                     "actual_coverage": k / n, "all_mse": float(error.mean()),
                     "retained_mse": retained, "retained_zero_mse": z,
                     "gain80_vs_all": bb.ratio_gain(retained, float(error.mean())),
                     "gain_vs_zero_same_tasks": bb.ratio_gain(retained, z),
                     "fraction_better_zero_retained": float(np.mean(error[ix] < zero[ix])),
                     "AURC": float(bb.curve(error, order).mean()),
                     "spearman": bb.spearman(features[name], error)})
    return rows


def comparison(rows):
    a = next(r for r in rows if r["rule"] == "base_rank_mean")
    b = next(r for r in rows if r["rule"] == "extended_rank_mean")
    return {"extended_minus_base_retained_mse": b["retained_mse"] - a["retained_mse"],
            "relative_reduction_vs_base": bb.ratio_gain(b["retained_mse"], a["retained_mse"]),
            "extended_minus_base_AURC": b["AURC"] - a["AURC"]}


def verify_protocol(root, out):
    protocol = json.loads((out / "frozen_protocol.json").read_text())
    require(protocol["code_sha256"] == code_hashes(), "Code changed after freeze")
    require(protocol["manifest_sha256"] == bb.digest(root / "processed/v01/manifest.json"),
            "Manifest changed after freeze")
    for env, value in protocol["source_sha256"].items():
        require(bb.digest(source_path(root, env)) == value, "Source changed: " + env)
    for name, value in protocol["prediction_sha256"].items():
        require(bb.digest(out / name) == value, "Frozen prediction changed: " + name)
    return protocol


def stage_done(out, name):
    marker = out / (name + "_complete.json")
    if not marker.exists():
        return False
    data = json.loads(marker.read_text())
    require(data["protocol_sha256"] == bb.digest(out / "frozen_protocol.json"),
            "Stage protocol mismatch")
    for file, value in data["files"].items():
        require(bb.digest(out / file) == value, "Completed output changed: " + file)
    print("REUSE_OK " + name, flush=True)
    return True


def complete(out, name, files):
    json_write(out / (name + "_complete.json"),
               {"protocol_sha256": bb.digest(out / "frozen_protocol.json"),
                "files": {str(p.relative_to(out)): bb.digest(p) for p in files}})


def predict(root, out, repeats):
    if (out / "frozen_protocol.json").exists():
        p = verify_protocol(root, out)
        require(p["split_repeats"] == repeats, "Requested split repeats differ from frozen protocol")
        print("REUSE_OK predict", flush=True)
        return
    require(not out.exists(), "Incomplete/existing output directory. Inspect it; use a NEW --output path.")
    manifest_path = root / "processed/v01/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    require(manifest["script_sha256"] == bb.digest(Path(pp.__file__)),
            "Preprocessing script differs from the original manifest")
    out.mkdir(parents=True)
    sha = {e: bb.digest(source_path(root, e)) for e in ENVS}
    files = {}
    for fold, target in enumerate(ENVS):
        sources = [e for e in ENVS if e != target]
        # Only this fold's two sources are loaded before its prediction is saved.
        data = [load_source(root, e, manifest) for e in sources]
        pred, features = make_prediction([v[0] for v in data], [v[1] for v in data])
        tie = np.random.default_rng(SEED + fold).permutation(len(pred))
        path = out / (target + "_frozen.h5")
        pp.write_h5(path, {"perturbation": manifest["perturbations"],
                           "gene_id": manifest["gene_ids"], "prediction": pred,
                           "tie_key": tie, **features},
                    {"target": target, "sources": ",".join(sources),
                     "manifest_sha256": bb.digest(manifest_path)})
        files[path.name] = bb.digest(path)
        print("PREDICT_OK " + target + ": two-source prediction frozen", flush=True)
        del data, pred, features
    protocol = {"protocol": "VC-internal-v01", "created_utc": datetime.now(timezone.utc).isoformat(),
                "environments": list(ENVS), "rules": list(RULES), "seed": SEED,
                "split_repeats": repeats, "coverage": COVERAGE,
                "primary_comparison": "extended vs base retained MSE at ceil(0.8*N); lower is better",
                "secondary_metrics": ["AURC", "Spearman", "zero baseline on SAME retained tasks"],
                "combinations": "average-tie rank/N; equal weights; no fitted weights or thresholds",
                "predictor": "equal mean of TWO source responses for the same known perturbation",
                "split_design": "balanced random halves within each perturbation; controls within GEM group; independent RNG per environment; five repeats by default",
                "stability_modes": ["fixed_full_source_prediction_vs_half_target", "half_source_prediction_vs_half_target"],
                "inference_limits": ["developmental analysis motivated by revealed HepG2 pilot",
                                     "inherited cohort/genes selected using all FOUR environments",
                                     "no HepG2 response/raw file opened in this round",
                                     "folds and repeated splits overlap; not independent replicates",
                                     "half-split variation is NOT a full-sample bootstrap CI",
                                     "no causal inference or measurement-noise subtraction",
                                     "two-source findings require further three-source validation",
                                     "file freeze is computational provenance, not external preregistration"],
                "manifest_sha256": bb.digest(manifest_path), "source_sha256": sha,
                "prediction_sha256": files, "code_sha256": code_hashes(),
                "versions": {"python": sys.version.split()[0], "numpy": np.__version__,
                             "scipy": scipy.__version__, "h5py": h5py.__version__}}
    json_write(out / "frozen_protocol.json", protocol)


def read_frozen(out, target):
    with h5py.File(out / (target + "_frozen.h5"), "r") as f:
        return f["prediction"][:], {r: f[r][:] for r in RULES}, f["tie_key"][:]


def evaluate(root, out):
    verify_protocol(root, out)
    if stage_done(out, "evaluation"):
        return
    manifest = json.loads((root / "processed/v01/manifest.json").read_text())
    metrics, comparisons, details = [], [], []
    for target in ENVS:
        prediction, features, tie = read_frozen(out, target)
        truth, n = load_source(root, target, manifest)
        error = np.mean((prediction - truth) ** 2, axis=1)
        zero = np.mean(truth ** 2, axis=1)
        rows = metric_rows(error, zero, features, tie, {"target": target})
        metrics.extend(rows)
        c = {"target": target, **comparison(rows)}
        comparisons.append(c)
        for i, gene in enumerate(manifest["perturbations"]):
            details.append({"target": target, "perturbation": gene,
                            **{r: float(features[r][i]) for r in RULES},
                            "mse": float(error[i]), "zero_mse": float(zero[i]),
                            "target_n_cells_evaluation_only": int(n[i]), "tie_key": int(tie[i])})
        print("INTERNAL_OK " + target + ": extended vs base relative MSE reduction="
              + format(c["relative_reduction_vs_base"], ".4%"), flush=True)
    files = [out / n for n in ("internal_metrics.tsv", "internal_comparison.tsv", "internal_tasks.tsv")]
    for path, rows in zip(files, [metrics, comparisons, details]):
        bb.write_table(path, rows)
    complete(out, "evaluation", files)


def balanced_assignments(codes, repeats, rng):
    """Random, near-equal halves in each nonnegative stratum; -1 excluded."""
    result = np.full((repeats, len(codes)), -1, dtype=np.int8)
    valid = np.flatnonzero(codes >= 0)
    order = valid[np.argsort(codes[valid], kind="stable")]
    boundaries = np.r_[0, np.flatnonzero(np.diff(codes[order])) + 1, len(order)]
    for start, stop in zip(boundaries[:-1], boundaries[1:]):
        ids = order[start:stop]
        require(len(ids) >= 2, "A perturbation/control stratum has fewer than two cells")
        for r in range(repeats):
            perm = rng.permutation(ids)
            offset = int(rng.integers(0, 2))
            result[r, perm] = (np.arange(len(ids)) + offset) % 2
    return result


def split_environment(root, out, env, repeats, chunk_size):
    name = "split_" + env
    if stage_done(out, name):
        return
    manifest = json.loads((root / "processed/v01/manifest.json").read_text())
    info = manifest["inputs"][env]
    raw_path = Path(info["path"])
    stat = raw_path.stat()
    require(stat.st_size == info["size_bytes"] and stat.st_mtime_ns == info["mtime_ns"],
            "Raw input changed since prepare: " + env)
    labels, batches, gene_ids, umi = pp.metadata(raw_path)
    ref, full_n = load_source(root, env, manifest)
    with h5py.File(source_path(root, env), "r") as f:
        batch_names = bb.strings(f["gem_group"])
    P, G, B = len(manifest["perturbations"]), len(manifest["gene_ids"]), len(batch_names)
    K = P + B
    pmap = {v: i for i, v in enumerate(manifest["perturbations"])}
    bmap = {v: i for i, v in enumerate(batch_names)}
    gmap = {v: i for i, v in enumerate(gene_ids)}
    gene_ix = np.array([gmap[v] for v in manifest["gene_ids"]])
    pc = np.array([pmap.get(v, -1) for v in labels])
    bc = np.array([bmap.get(v, -1) for v in batches])
    control = labels == pp.CONTROL
    valid = (bc >= 0) & ((pc >= 0) | control)
    # Every control cell has ONE half assignment reused by all perturbations.
    strata = np.where(valid, np.where(control, P + bc, pc), -1)
    rng = np.random.default_rng(SEED + 1000 * (ENVS.index(env) + 1))
    half = balanced_assignments(strata, repeats, rng)
    sums = np.zeros((repeats, 2 * K, G), dtype=np.float64)
    nc = np.zeros((repeats, 2, B), dtype=np.int64)
    nb = np.zeros((repeats, 2, P, B), dtype=np.int64)
    with h5py.File(raw_path, "r") as f:
        for start in range(0, len(labels), chunk_size):
            end = min(start + chunk_size, len(labels))
            raw = f["X"][start:end, :]
            pp.validate_block(raw, umi[start:end])
            keep = valid[start:end]
            if keep.any():
                y = pp.normalized(raw[keep], umi[start:end][keep], gene_ix)
                b, p, ctrl = bc[start:end][keep], pc[start:end][keep], control[start:end][keep]
                agg = np.where(ctrl, P + b, p)
                for r in range(repeats):
                    side = half[r, start:end][keep]
                    pp.accumulate(sums[r], side * K + agg, y)
                    np.add.at(nc[r], (side[ctrl], b[ctrl]), 1)
                    np.add.at(nb[r], (side[~ctrl], p[~ctrl], b[~ctrl]), 1)
            if start == 0 or end == len(labels) or start // chunk_size % 100 == 0:
                print(f"SPLIT_SCAN {env}: {end:,}/{len(labels):,}", flush=True)
    path = out / (env + "_halves.h5")
    temp = path.with_name(path.name + ".partial")
    repeatability = []
    with h5py.File(temp, "x") as f:
        f.attrs["protocol_sha256"] = bb.digest(out / "frozen_protocol.json")
        f.attrs["environment"] = env
        for name_ids, values in [("perturbation", manifest["perturbations"]), ("gene_id", manifest["gene_ids"])]:
            f.create_dataset(name_ids, data=np.asarray(values, dtype=object), dtype=h5py.string_dtype())
        ds = f.create_dataset("delta", shape=(repeats, 2, P, G), dtype="f4",
                              chunks=(1, 1, min(P, 16), G), compression="gzip", compression_opts=1)
        n = nb.sum(axis=3)
        require((n > 0).all() and (nc > 0).all(), "Empty half-stratum")
        require(np.all(n.sum(axis=1) == full_n[None, :]), "Split cells do not reconcile")
        f.create_dataset("n_cells", data=n)
        f.create_dataset("n_cells_by_batch", data=nb, compression="gzip")
        f.create_dataset("control_n_by_batch", data=nc)
        max_diff = 0.0
        for r in range(repeats):
            halves = []
            full_sums = sums[r, :K] + sums[r, K:]
            full_ctrl = full_sums[P:] / nc[r].sum(axis=0)[:, None]
            reconstructed = full_sums[:P] / full_n[:, None] - (nb[r].sum(axis=0) / full_n[:, None]) @ full_ctrl
            discrepancy = float(np.max(np.abs(reconstructed - ref)))
            require(np.allclose(reconstructed, ref, atol=2e-6, rtol=2e-6),
                    "Raw reconstruction differs from v01 preprocessing")
            max_diff = max(max_diff, discrepancy)
            for side in range(2):
                s = sums[r, side * K:(side + 1) * K]
                ctrl = s[P:] / nc[r, side, :, None]
                d = s[:P] / n[r, side, :, None] - (nb[r, side] / n[r, side, :, None]) @ ctrl
                ds[r, side] = d.astype(np.float32)
                halves.append(d)
            corr = bb.pearson_rows(*halves)
            hmse = np.mean((halves[0] - halves[1]) ** 2, axis=1)
            for i, p in enumerate(manifest["perturbations"]):
                repeatability.append({"environment": env, "repeat": r + 1, "perturbation": p,
                                      "n_half0": int(n[r, 0, i]), "n_half1": int(n[r, 1, i]),
                                      "half_response_pearson": float(corr[i]),
                                      "half_difference_mse": float(hmse[i])})
        f.attrs["full_reconstruction_max_abs_difference"] = max_diff
    bb.commit(temp, path)
    table = out / (env + "_repeatability.tsv")
    bb.write_table(table, repeatability)
    complete(out, name, [path, table])
    print(f"SPLIT_OK {env}: {repeats} paired splits; full-response reconstruction verified", flush=True)


def split(root, out, chunk_size):
    protocol = verify_protocol(root, out)
    for env in ENVS:
        split_environment(root, out, env, protocol["split_repeats"], chunk_size)


def stability(root, out):
    protocol = verify_protocol(root, out)
    if stage_done(out, "stability"):
        return
    for env in ENVS:
        require(stage_done(out, "split_" + env), "Run split first")
    rows, comps = [], []
    for target in ENVS:
        full_prediction, full_features, tie = read_frozen(out, target)
        sources = [e for e in ENVS if e != target]
        for r in range(protocol["split_repeats"]):
            for side in range(2):
                data = {}
                for env in sources:
                    with h5py.File(out / (env + "_halves.h5"), "r") as f:
                        data[env] = (f["delta"][r, side].astype(np.float64), f["n_cells"][r, side])
                hp, hf = make_prediction([data[e][0] for e in sources], [data[e][1] for e in sources])
                with h5py.File(out / (target + "_halves.h5"), "r") as f:
                    truth = f["delta"][r, side].astype(np.float64)
                zero = np.mean(truth ** 2, axis=1)
                for mode, prediction, features in [("fixed_sources", full_prediction, full_features),
                                                     ("half_sources", hp, hf)]:
                    error = np.mean((prediction - truth) ** 2, axis=1)
                    prefix = {"target": target, "repeat": r + 1, "half": side, "mode": mode}
                    rr = metric_rows(error, zero, features, tie, prefix)
                    rows.extend(rr)
                    comps.append({**prefix, **comparison(rr),
                                  "H_rank_correlation_to_full_sources": bb.spearman(features["source_heterogeneity"], full_features["source_heterogeneity"])})
                del data, hp, hf, truth
    files = [out / "stability_metrics.tsv", out / "stability_comparison.tsv", out / "stability_summary.tsv"]
    summary = []
    for target in ENVS:
        for mode in ("fixed_sources", "half_sources"):
            x = [v for v in comps if v["target"] == target and v["mode"] == mode]
            diff = np.array([v["extended_minus_base_retained_mse"] for v in x])
            summary.append({"target": target, "mode": mode, "n_half_evaluations": len(x),
                            "n_extended_better": int((diff < 0).sum()),
                            "difference_min": float(diff.min()), "difference_median": float(np.median(diff)),
                            "difference_max": float(diff.max()),
                            "interpretation": "descriptive half-split sensitivity; NOT a confidence interval"})
    for file, rr in zip(files, [rows, comps, summary]):
        bb.write_table(file, rr)
    complete(out, "stability", files)
    print("STABILITY_OK: paired half-split sensitivity completed; no independence/CI claim", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["predict", "evaluate", "split", "stability", "run"])
    p.add_argument("--root", type=Path, default=Path.home() / "VC_reliability")
    p.add_argument("--output", type=Path)
    p.add_argument("--split-repeats", type=int, default=5, help="Frozen at predict; default 5")
    p.add_argument("--chunk-size", type=int, default=512)
    args = p.parse_args()
    require(1 <= args.split_repeats <= 10 and args.chunk_size > 0, "Invalid resource parameters")
    root = args.root.resolve()
    out = args.output.resolve() if args.output else root / "results/internal_v01"
    if args.mode in ("predict", "run"):
        predict(root, out, args.split_repeats)
    if args.mode in ("evaluate", "run"):
        evaluate(root, out)
    if args.mode in ("split", "run"):
        split(root, out, args.chunk_size)
    if args.mode in ("stability", "run"):
        stability(root, out)
    if args.mode == "run":
        print("ALL_OK: " + str(out), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR " + type(exc).__name__ + ": " + str(exc), file=sys.stderr, flush=True)
        sys.exit(1)
