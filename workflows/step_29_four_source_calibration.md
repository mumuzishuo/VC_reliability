# Step 29: calibrate four internal sources

This development run uses the audited v01 responses from K562 essential,
RPE1, Jurkat and HepG2 on the internal panel of 1,264 tasks and 6,642
expression genes. It verifies the response and preprocessing manifest hashes
before reading the four matrices. The HepG2 file keeps its `evaluation_only`
role from the old experiment. Its use as a development source is specific
to the new cross-study prediction. The run has no X-Atlas input.

## Predictions and calibration

The paired predictor averages the four responses for the same perturbation.
The common predictor averages those vectors across all source perturbations
and assigns the same response profile to every task. Each source contributes
equally. A zero vector provides the reference.

For each source context `s`, build a predictor `m[-s]` from the remaining
three and compare it with response `y[s]`. For each predictor separately,
the constrained least-squares coefficient is:

```text
scale = clip(sum_s <m[-s], y[s]> / sum_s <m[-s], m[-s]>, 0, 1)
```

The inner product sums over the same tasks and expression genes in each
fold. For a zero denominator, use scale zero. Apply the two fitted scales
to the four-source paired and common predictors, respectively. With two
sources this formula matches the archived two-direction calibration.

To evaluate the calibration internally, hold out each context and repeat
the fit among the other three. Each inner calibration fold uses two sources
to predict the third. The outer holdout contributes to scoring. Its response
and cell support are absent from the outer prediction, scale and risk inputs.
The three-to-one internal scores are retrospective and describe the inherited
panel. The four-to-one performance remains unmeasured.

## Risk diagnostics

Compare the three existing scores (source heterogeneity, inverse minimum
source cell count, and raw paired response RMS) with one candidate score.
The candidate is a task's mean excess MSE over zero in held-out source
contexts. For each of those contexts, estimate the predictor and its scale
from the remaining sources. A higher value is a proposed signal for failure.

The event is `D > 0`, with `D = paired_calibrated MSE - zero_response MSE`.
Report AUROC, average precision and failure prevalence. Retain the lowest
scores for `ceil(0.8 * n_tasks)` tasks with a fixed random tie order. Report
MSE for that subset and MSE across all tasks with zero response on rejected
tasks. Compare the latter with always predicting, always using zero and random
selection at the same coverage. This run reports every candidate without
selecting an external risk rule. Split-half and retrieval will be part of the
locked external evaluation plan.

## Run on the HPC

```bash
cd "$HOME/VC_reliability_repo" &&
git pull --ff-only origin dev &&
mkdir -p logs/v01_four_source_calibration &&
sbatch slurm/step_29_four_source_calibration.slurm
```

The job uses account `zhenwei01_bio214`, partition `cpu8358`, QOS `26cores`,
2 CPUs and 8 GB. Allow a few minutes after it starts, with a 30 minute limit.
It reads the previously audited summaries, about 340 MB of H5 files. The
step runs with the HPC Python 3.9 module.

After completion, set `vc_job_id` to the numeric Slurm job ID. From the
repository root, check the scheduler and output checksums:

```bash
vc_job_id=YOUR_JOB_ID
vc_out="results/v01_four_source_calibration/job_$vc_job_id"
sacct -X -j "$vc_job_id" --format=JobID,State,ExitCode,Elapsed
sha256sum --status -c "$vc_out/SHA256SUMS" && cat "$vc_out/summary.json"
```

Look for `COMPLETED 0:0` and `FOUR_SOURCE_CALIBRATION_COMPLETE` in the
scheduler and stdout, respectively. If either is missing, inspect the
matched `.out` and `.err` files in `logs/v01_four_source_calibration/`.

When the checks pass, archive the small result tables, parameters and log:

```bash
(
  set -euo pipefail
  git pull --ff-only origin dev
  sha256sum --status -c "$vc_out/SHA256SUMS"
  vc_files=("$vc_out"/{summary.json,parameters.json,code_commit.txt,SHA256SUMS}
            "$vc_out"/internal_{model_metrics,risk_metrics}.tsv
            "logs/v01_four_source_calibration/step29_$vc_job_id.out")
  vc_public="manifests/v01_four_source_calibration_${vc_job_id}.sha256"
  sha256sum -- "${vc_files[@]}" > "$vc_public"
  sha256sum --status -c "$vc_public"
  git add -- "${vc_files[@]}" "$vc_public"
  git -c core.whitespace=cr-at-eol diff --cached --check -- "${vc_files[@]}" "$vc_public"
  git commit --only -m "results: archive four-source development calibration" -- \
    "${vc_files[@]}" "$vc_public"
  git push origin dev
)
```

The full local run retains `predictions.h5`, per-task errors and a SHA256
manifest for all outputs. The H5 stores unscaled paired and common vectors,
the scales, source cell counts and risk scores. The full internal panel is
for development. Before external evaluation, fix the external task and
feature axes from metadata, freeze the evaluation and risk rules, and save
prediction hashes for those axes before accessing the external responses.

## HPC run 3239681

The user reported `COMPLETED 0:0` in 21 seconds, `SHA256_OK` for the
run's output manifest, and `FOUR_SOURCE_CALIBRATION_COMPLETE 3239681`
from commit `7e64c3d`. The run used 1,264 tasks and 6,642 expression
genes. The four-source scale was **0.6019758093** for the paired predictor
and **0.7550802413** for the common predictor. These scales describe the
full inherited panel. External task and feature axes remain to be fixed.

| Held-out context | Raw paired MSE (reconstructed) | Calibrated paired MSE | Zero MSE | Calibrated gain vs zero | Calibrated tasks better than zero |
| --- | ---: | ---: | ---: | ---: | ---: |
| K562 | 0.004115 | 0.002675 | 0.003219 | 16.89% | 600/1,264 (47.47%) |
| RPE1 | 0.007629 | 0.008246 | 0.010121 | 18.52% | 925/1,264 (73.18%) |
| Jurkat | 0.005094 | 0.004166 | 0.005081 | 18.00% | 758/1,264 (59.97%) |
| HepG2 | 0.004972 | 0.004799 | 0.006185 | 22.41% | 859/1,264 (67.96%) |

All four nested internal folds have lower mean calibrated MSE than zero.
In K562, 664 tasks are worse than zero and 600 improve. The net reduction
in mean MSE therefore coexists with a failure rate above one half.
Across the four equally weighted internal folds, the reconstructed raw
paired MSE is 0.00545279 and the calibrated paired MSE is 0.00497150.
Calibration reduces MSE in three folds. RPE1's calibrated MSE increases 8.08% relative to its unscaled paired
predictor, while remaining below zero-response MSE.

Calibrated and zero MSEs come directly from the HPC summary. Counts are
recovered from the reported task fractions and the known denominator.
The summary also allows reconstruction of raw paired MSE. Let `L` be
the squared norm of the sum of the four response matrices divided by
`P*G`, where `P` is the task count and `G` is the gene count.
`L = sum_e(z_e) + 3 * paired_cross_product / (P*G)`, with `z_e` the
fold's zero MSE. For the raw prediction `m_e` in fold `e`, let
`A_e = <m_e,m_e>/(P*G)` and `B_e = <m_e,truth_e>/(P*G)`.
Solve `9*A_e + 6*B_e = L - z_e` and
`a_e^2*A_e - 2*a_e*B_e = calibrated_MSE_e - z_e`, then compute
`raw_MSE_e = z_e + A_e - 2*B_e`. Here `a_e` is the reported
paired scale for the fold. The mean reconstructed raw MSE agrees
with the pooled raw MSE calculated from the reported cross product and
predictor sum of squares to within `1e-15`.

Archive the actual run files using the command above with `vc_job_id=3239681`.
Check the reconstructed values against `internal_model_metrics.tsv`.
Use that table for the comparison with `common_calibrated`. Use
`internal_risk_metrics.tsv` to assess `D > 0` discrimination and the
all-task MSE of the 80% policy. The pasted summary contains neither of
those comparisons. External X-Atlas responses remain sealed.
