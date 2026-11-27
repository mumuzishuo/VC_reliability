#!/usr/bin/env python3
"""Retrospective diagnostics for the three existing VC-internal-v01 folds."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
from scipy.stats import rankdata, spearmanr

import vc_baseline as bb
import vc_internal as vi
import vc_preprocess as pp

SEED = 20261105
PERMUTATIONS = 200
BOOTSTRAPS = 1000
require = bb.require


def ms(x):
    return np.mean(x * x, axis=1)


def summarize(e, z):
    return dict(n_tasks=len(e), mse=float(e.mean()), zero_mse=float(z.mean()),
                mean_D=float((e - z).mean()), fraction_worse_zero=float(np.mean(e > z)),
                fraction_better_zero=float(np.mean(e < z)))


def raw_audit(root, env, manifest, truth, n):
    """Independent float64 normalization and direct batch-matched averaging."""
    info = manifest["inputs"][env]
    raw_path = root / "raw" / Path(info["path"]).name
    stat = raw_path.stat()
    require(stat.st_size == info["size_bytes"] and stat.st_mtime_ns == info["mtime_ns"],
            f"{env}: raw file differs from manifest")
    labels, batches, ids, umi = pp.metadata(raw_path)
    perturbations = np.asarray(manifest["perturbations"])
    genes = manifest["gene_ids"]
    index = {g: i for i, g in enumerate(ids)}
    cols = np.asarray([index[g] for g in genes])
    ctrl = labels == manifest["control_label"]
    names, count = np.unique(batches[ctrl], return_counts=True)
    eligible = names[count >= manifest["min_controls"]]
    raw_ctrl_n = dict(zip(names, count))
    valid = np.isin(batches, eligible)
    with h5py.File(vi.source_path(root, env), "r") as h:
        batch_names = bb.strings(h["gem_group"])
        pindex = {p: i for i, p in enumerate(bb.strings(h["perturbation"]))}
        gindex = {g: i for i, g in enumerate(bb.strings(h["gene_id"]))}
        stored_cols = [gindex[g] for g in genes]
        require(set(batch_names) == set(eligible), f"{env}: eligible GEM groups differ")
        require(np.array_equal(h["control_cell_count_by_batch"][:],
                               [raw_ctrl_n[b] for b in batch_names]),
                f"{env}: control counts differ")
        by_batch = h["n_cells_by_batch"][:][[pindex[p] for p in perturbations]]
        control_ref = h["control_mean_by_batch"][:][:, stored_cols]
        perturb_ref = h["perturbed_mean_log1p"][:]
        matched_ref = h["matched_control_mean_log1p"][:]
    low = np.asarray([raw_ctrl_n[b] < 30 for b in batch_names])
    low_fraction = by_batch[:, low].sum(axis=1) / n
    # Fixed sample: extremes and middle of cell support and response size,
    # plus the tasks with the largest low-control share and single-GEM share.
    selected = set()
    for value in (n, ms(truth)):
        order = np.argsort(value, kind="stable")
        selected.update(order[np.linspace(0, len(n) - 1, 5, dtype=int)])
    selected.update((np.argmax(low_fraction), np.argmax(by_batch.max(axis=1) / n)))
    selected = sorted(selected)
    needed = batch_names[(by_batch[selected].sum(axis=0) > 0)]
    rows_needed = (valid & np.isin(labels, perturbations[selected])) | (ctrl & np.isin(batches, needed))
    control_sums = {b: np.zeros(len(genes)) for b in needed}
    perturb_sums = np.zeros((len(selected), len(genes)))
    observed_n = np.zeros((len(selected), len(batch_names)), dtype=int)
    observed_ctrl = dict.fromkeys(needed, 0)
    ratio_min, ratio_max = np.inf, -np.inf
    with h5py.File(raw_path, "r") as raw_file:
        x = raw_file["X"]
        for start in range(0, len(labels), 512):
            end = min(start + 512, len(labels))
            keep = rows_needed[start:end]
            if not keep.any():
                continue
            raw = x[start:end, :][keep]
            total = umi[start:end][keep]
            ratio = raw.sum(axis=1, dtype=np.float64) / total
            require(np.isfinite(raw).all() and (raw >= 0).all()
                    and np.equal(raw, np.floor(raw)).all() and (ratio <= 1.0001).all(),
                    f"{env}: invalid raw count or UMI denominator")
            ratio_min, ratio_max = min(ratio_min, ratio.min()), max(ratio_max, ratio.max())
            y = np.log1p(10000.0 * raw[:, cols].astype(np.float64) / total[:, None])
            lb, bt = labels[start:end][keep], batches[start:end][keep]
            for batch in np.unique(bt):
                ix = bt == batch
                control_rows = ix & (lb == manifest["control_label"])
                if control_rows.any():
                    control_sums[batch] += y[control_rows].sum(axis=0)
                    observed_ctrl[batch] += int(control_rows.sum())
                j = int(np.where(batch_names == batch)[0][0])
                for i, p in enumerate(perturbations[selected]):
                    pr = ix & (lb == p)
                    if pr.any():
                        perturb_sums[i] += y[pr].sum(axis=0)
                        observed_n[i, j] += int(pr.sum())
    require(np.array_equal(observed_n, by_batch[selected]), f"{env}: task GEM counts differ")
    require(np.array_equal(observed_n.sum(axis=1), n[selected]), f"{env}: task cell counts differ")
    require(all(observed_ctrl[b] == raw_ctrl_n[b] for b in needed), f"{env}: raw controls differ")
    control_means = np.zeros_like(control_ref, dtype=np.float64)
    used = np.isin(batch_names, needed)
    for j in np.flatnonzero(used):
        b = batch_names[j]
        control_means[j] = control_sums[b] / observed_ctrl[b]
    require(np.allclose(control_means[used], control_ref[used], atol=2e-6, rtol=2e-6),
            f"{env}: control means differ")
    out = []
    for i, p in enumerate(perturbations[selected]):
        row = selected[i]
        mean = perturb_sums[i] / n[row]
        matched = (observed_n[i] / n[row]) @ control_means
        pairs = {"perturb_mean": (mean, perturb_ref[pindex[p], stored_cols]),
                 "matched_control": (matched, matched_ref[pindex[p], stored_cols]),
                 "delta": (mean - matched, truth[row])}
        differences = {}
        for name, (value, reference) in pairs.items():
            require(np.allclose(value, reference, atol=2e-6, rtol=2e-6),
                    f"{env}/{p}: {name} differs")
            differences[f"{name}_max_abs_diff"] = float(np.max(np.abs(value - reference)))
        out.append(dict(target=env, perturbation=p, n_cells=int(n[row]),
                        gem_groups=int(np.count_nonzero(observed_n[i])),
                        low_control_fraction=float(low_fraction[row]),
                        stored_counts_over_UMI_min=float(ratio_min),
                        stored_counts_over_UMI_max=float(ratio_max), **differences))
    print(f"RAW_AUDIT_OK {env}: {len(selected)} tasks, {int(used.sum())} GEMs; "
          f"low-control GEMs (<30)={int(low.sum())}; "
          f"max delta difference={max(row['delta_max_abs_diff'] for row in out):.3g}", flush=True)
    return out, low_fraction


def read_halves(internal, env, manifest, pred=None):
    with h5py.File(internal / f"{env}_halves.h5", "r") as f:
        require(f.attrs["environment"] == env and
                f.attrs["protocol_sha256"] == bb.digest(internal / "frozen_protocol.json"),
                f"{env}: half-split provenance mismatch")
        for key, expected in (("perturbation", manifest["perturbations"]),
                              ("gene_id", manifest["gene_ids"])):
            require(np.array_equal(bb.strings(f[key]), expected), f"{env}: half {key} order differs")
        d = f["delta"]
        corr, half_mse, half_d = [], [], []
        for r in range(d.shape[0]):
            a, b = d[r].astype(np.float64)
            require(np.isfinite(a).all() and np.isfinite(b).all(), f"{env}: nonfinite halves")
            corr.append(bb.pearson_rows(a, b))
            half_mse.append(ms(a - b))
            if pred is not None:
                half_d.extend([ms(pred - a) - ms(a), ms(pred - b) - ms(b)])
    with np.errstate(invalid="ignore"):
        finite = np.isfinite(corr).sum(axis=0)
    median = np.ma.median(np.ma.masked_invalid(corr), axis=0).filled(np.nan)
    return median, np.median(half_mse, axis=0), finite, np.asarray(half_d)


def diagnostic(pred, truth):
    a, b = np.sqrt(ms(pred)), np.sqrt(ms(truth))
    inner = np.einsum("ij,ij->i", pred, truth) / truth.shape[1]
    ab = a * b
    cosine = np.divide(inner, ab, out=np.full_like(a, np.nan), where=ab > 0)
    amplitude = (a - b) ** 2
    direction = np.maximum(2 * (ab - inner), 0)
    e, z = ms(pred - truth), ms(truth)
    require(np.allclose(amplitude + direction, e, atol=1e-12, rtol=1e-10),
            "MSE decomposition failed")
    return a, b, cosine, amplitude, direction, e, z


def permutation_rows(target, pred, truth, e, z, rng):
    # MSE after permutation from the pairwise inner products. The marginal
    # prediction and truth energies are unchanged by any permutation.
    dot = (pred @ truth.T) / truth.shape[1]
    total = ms(pred).mean() + z.mean()
    order = np.argsort(ms(pred), kind="stable")
    groups = np.array_split(order, 5)
    rows = []
    for r in range(PERMUTATIONS):
        shuffle = rng.permutation(len(e))
        within = np.arange(len(e))
        for ix in groups:
            within[ix] = rng.permutation(ix)
        for kind, perm in (("unrestricted", shuffle), ("within_predicted_rms_quintile", within)):
            mse = total - 2 * dot[perm, np.arange(len(e))].mean()
            rows.append(dict(target=target, permutation=r + 1, scheme=kind,
                             mse=float(mse), mse_minus_paired=float(mse - e.mean())))
    return rows


def classifier_metrics(score, failure):
    """AUROC with average ties, AUPRC as non-interpolated average precision."""
    positive = int(failure.sum())
    negative = len(failure) - positive
    if not positive or not negative:
        return np.nan, np.nan
    ranks = rankdata(score)
    auc = (ranks[failure].sum() - positive * (positive + 1) / 2) / (positive * negative)
    order = np.argsort(-score, kind="stable")
    s, y = score[order], failure[order]
    end = np.r_[np.flatnonzero(s[1:] != s[:-1]), len(s) - 1]
    tp = np.cumsum(y)[end]
    ap = np.sum(np.diff(np.r_[0, tp]) * tp / (end + 1)) / positive
    return float(auc), float(ap)


def risk_metrics(target, scores, tie, e, z, half_d, rng):
    k = math.ceil(0.8 * len(e))
    d = e - z
    failure = d > 0
    boot = rng.integers(0, len(e), size=(BOOTSTRAPS, len(e)))
    rows = []
    for name, score in scores.items():
        order = bb.risk_order(score, tie)
        retained = order[:k]
        use = np.zeros(len(e), dtype=bool)
        use[retained] = True
        policy = np.where(use, e, z)
        # Fixed actions; only the evaluated tasks are resampled. These are
        # conditional descriptive intervals, not independent-screen inference.
        vs_zero = policy - z
        lower, upper = np.quantile(vs_zero[boot].mean(axis=1), [0.025, 0.975])
        auc, ap = classifier_metrics(score, failure)
        hd = half_d[:, use].sum(axis=1) / len(e)
        rows.append(dict(target=target, rule=name, coverage=k / len(e),
                         n_tasks=len(e), n_retained=k, failure_prevalence=float(failure.mean()),
                         spearman_E=float(spearmanr(score, e)[0]),
                         spearman_D=float(spearmanr(score, d)[0]), auroc_D_gt_0=auc, auprc_AP=ap,
                         AURC=float(bb.curve(e, order).mean()),
                         retained_mse=float(e[retained].mean()),
                         retained_zero_mse=float(z[retained].mean()),
                         retained_failure_rate=float(failure[retained].mean()),
                         all_prediction_mse=float(e.mean()), all_zero_mse=float(z.mean()),
                         policy_mse=float(policy.mean()),
                         random_same_coverage_expected_mse=float(z.mean() + use.mean() * d.mean()),
                         policy_minus_all_prediction=float((policy - e).mean()),
                         policy_minus_all_zero=float(vs_zero.mean()),
                         policy_minus_zero_boot_conditional_lo=float(lower),
                         policy_minus_zero_boot_conditional_hi=float(upper),
                         half_policy_minus_zero_min=float(hd.min()),
                         half_policy_minus_zero_max=float(hd.max())))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.home() / "VC_reliability")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root
    internal = root / "results/internal_v01"
    out = args.output or root / "results/diagnostics_v01"
    require(not out.exists(), f"Output exists: {out}. Choose a new --output.")
    protocol = vi.verify_protocol(root, internal)
    for env in vi.ENVS:
        require(vi.stage_done(internal, "split_" + env), f"{env}: missing completed half splits")
    manifest = json.loads((root / "processed/v01/manifest.json").read_text())
    cohorts = manifest["perturbations"]
    data = {env: vi.load_source(root, env, manifest) for env in vi.ENVS}
    raw_rows, low_fractions = [], {}
    for env, (truth, n) in data.items():
        rows, low_fractions[env] = raw_audit(root, env, manifest, truth, n)
        raw_rows.extend(rows)
    tasks, strata, models, permutations, risks = [], [], [], [], []
    rng = np.random.default_rng(SEED)
    repeatability = {}
    # Read target halves once per environment and reuse the resulting summaries.
    for env in vi.ENVS:
        pred, _, _ = vi.read_frozen(internal, env)
        repeatability[env] = read_halves(internal, env, manifest, pred)
    for target in vi.ENVS:
        pred, scores, tie = vi.read_frozen(internal, target)
        truth, n = data[target]
        src = [env for env in vi.ENVS if env != target]
        src_arrays = [data[env][0] for env in src]
        require(np.array_equal(pred, (src_arrays[0] + src_arrays[1]) / 2),
                f"{target}: frozen prediction differs from the source mean")
        a, b, cosine, amp, angle, e, z = diagnostic(pred, truth)
        corr, half_mse, finite, half_d = repeatability[target]
        source_noise_scale = sum(repeatability[s][1] / 4 for s in src) / 4
        strength_order = bb.risk_order(a, tie)
        other_scores = {**scores, "exploratory_low_rms_as_failure_risk": -a}
        risks.extend(risk_metrics(target, other_scores, tie, e, z, half_d, rng))
        fold_permutations = permutation_rows(target, pred, truth, e, z, rng)
        permutations.extend(fold_permutations)
        model_vectors = [(src[0], src_arrays[0]), (src[1], src_arrays[1]),
                         ("paired_source_mean", pred),
                         ("common_source_mean", np.broadcast_to(pred.mean(axis=0), pred.shape)),
                         ("zero", np.zeros_like(pred))]
        common_mse = np.nan
        for label, vector in model_vectors:
            pa, _, cs, da, dd, pe, pz = diagnostic(vector, truth)
            good_cosine = cs[np.isfinite(cs)]
            models.append(dict(target=target, model=label, **summarize(pe, pz),
                               median_predicted_rms=float(np.median(pa)),
                               median_observed_rms=float(np.median(b)),
                               median_cosine=float(np.median(good_cosine)) if len(good_cosine) else np.nan,
                               mean_amplitude_term=float(da.mean()), mean_direction_term=float(dd.mean())))
            if label == "common_source_mean":
                common_mse = pe.mean()
        for i, p in enumerate(cohorts):
            tasks.append(dict(target=target, perturbation=p, target_cells=int(n[i]),
                              min_source_cells=min(data[s][1][i] for s in src),
                              target_low_control_fraction=low_fractions[target][i],
                              predicted_rms=a[i], observed_rms=b[i], cosine=cosine[i],
                              amplitude_term=amp[i], direction_term=angle[i],
                              mse=e[i], zero_mse=z[i], D=e[i] - z[i],
                              half_pearson_median=corr[i], half_mse_median=half_mse[i],
                              n_defined_half_pearson=int(finite[i]),
                              half_D_fraction_positive=float(np.mean(half_d[:, i] > 0)),
                              source_heterogeneity=scores["source_heterogeneity"][i],
                              source_half_noise_scale_heuristic=source_noise_scale[i]))
        defined = np.flatnonzero(np.isfinite(corr))

        def add_stratum(factor, group, ix):
            if len(ix):
                valid_corr = corr[ix][np.isfinite(corr[ix])]
                by_h = bb.risk_order(scores["source_heterogeneity"][ix], tie[ix])
                top = ix[by_h[:math.ceil(0.8 * len(ix))]]
                strata.append(dict(target=target, factor=factor, group=group,
                                   **summarize(e[ix], z[ix]), median_cells=float(np.median(n[ix])),
                                   median_observed_rms=float(np.median(b[ix])),
                                   median_half_pearson=float(np.median(valid_corr)) if len(valid_corr) else np.nan,
                                   median_half_mse=float(np.median(half_mse[ix])),
                                   median_source_heterogeneity=float(np.median(scores["source_heterogeneity"][ix])),
                                   median_source_half_noise_scale=float(np.median(source_noise_scale[ix])),
                                   retained80_n=len(top), retained80_mse_ranked_by_H=float(e[top].mean()),
                                   retained80_mse_minus_random=float(e[top].mean() - e[ix].mean())))
        for q, ix in enumerate(np.array_split(defined[np.lexsort((tie[defined], corr[defined]))], 4), 1):
            add_stratum("target_half_pearson", f"Q{q}", ix)
        add_stratum("target_half_pearson", "undefined", np.flatnonzero(~np.isfinite(corr)))
        for q, ix in enumerate(np.array_split(strength_order, 5), 1):
            add_stratum("predicted_rms", f"Q{q}", ix)
        add_stratum("target_low_control", "absent", np.flatnonzero(low_fractions[target] == 0))
        add_stratum("target_low_control", "present", np.flatnonzero(low_fractions[target] > 0))
        support_order = bb.risk_order(-scores["low_source_support"], tie)
        for q, ix in enumerate(np.array_split(support_order, 4), 1):
            add_stratum("min_source_cells", f"Q{q}", ix)
        for cell_high in (False, True):
            for rms_high in (False, True):
                ix = defined[((n[defined] >= np.median(n)) == cell_high)
                             & ((b[defined] >= np.median(b)) == rms_high)]
                ranked = ix[np.lexsort((tie[ix], corr[ix]))]
                for high, part in enumerate(np.array_split(ranked, 2)):
                    add_stratum("target_cells_observed_rms_half_pearson",
                                f"cellsHigh{int(cell_high)}_rmsHigh{int(rms_high)}_pearsonHigh{high}",
                                part)
        print(f"DIAGNOSTICS_OK {target}: source-mean MSE={e.mean():.8f}, "
              f"zero={z.mean():.8f}, worse={np.mean(e > z):.2%}", flush=True)
        shuffled = [r["mse"] for r in fold_permutations
                    if r["scheme"] == "within_predicted_rms_quintile"]
        print(f"CONTROLS_OK {target}: common-source MSE={common_mse:.8f}; "
              f"strength-matched shuffled median MSE={np.median(shuffled):.8f}", flush=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".vc_diagnostics_", dir=out.parent) as tmp:
        stage = Path(tmp) / out.name
        stage.mkdir()
        for name, rows in (("raw_audit", raw_rows), ("task_diagnostics", tasks),
                           ("strata", strata), ("model_comparison", models),
                           ("permutations", permutations), ("risk_metrics", risks)):
            bb.write_table(stage / (name + ".tsv"), rows)
        (stage / "protocol.json").write_text(json.dumps({
            "scope": "retrospective internal K562, RPE1, Jurkat",
            "frozen_protocol_sha256": bb.digest(internal / "frozen_protocol.json"),
            "diagnostics_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "seed": SEED, "permutations_per_scheme": PERMUTATIONS,
            "bootstrap_replicates": BOOTSTRAPS, "split_repeats": protocol["split_repeats"],
            "conditional_bootstrap": "fixed actions, tasks resampled as if exchangeable; shared controls can invalidate nominal coverage",
            "half_splits": "sampling sensitivity, not independent replicates",
            "source_half_noise_scale_heuristic": "sum of median source half-difference MSEs / 16; equal-half independent noise approximation only",
            "utility_score_orientation": "original high-risk orientations plus exploratory reversed predicted RMS",
            "raw_audit": "float64 comparison on deterministic stress cases; no alternate normalization fitted",
        }, indent=2) + "\n")
        os.rename(stage, out)
    print(f"ALL_OK {out}", flush=True)


if __name__ == "__main__":
    main()
