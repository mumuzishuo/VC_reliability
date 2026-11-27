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
