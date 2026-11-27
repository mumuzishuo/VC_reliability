# Step 30 · Prepare X-Atlas external predictions

This job projects the audited four-source predictions from job `3239681`
onto a candidate X-Atlas panel. It records the model coefficients, risk
ranking and eligibility rules before the target responses are available.
The inputs are the original Step 29 prediction H5 and two archived metadata
tables. The code has no external expression input.

## Candidate cohort

HCT116 is the primary external target. Its guide audit found no GEM with
fewer than 30 non-targeting controls. Its published expression file is also
smaller than the HEK293T file. HEK293T is a second test context within the
same study. Report its results separately without changing the method after
HCT116 is scored. The two cell lines represent two test contexts.

- Use the same Ensembl target IDs in both contexts. A task must have at least
  30 guide-validated cells in GEMs with at least 30 non-targeting cells in
  each context. The current guide metadata give 1,110 HCT116 tasks, 1,130
  HEK293T tasks and a shared candidate panel of 1,088. The shared panel
  omits 22 HCT116-only and 42 HEK293T-only eligible tasks so the main
  comparison uses the same tasks in both contexts.
- Use the 6,638 exact unique expression IDs from the audited feature
  dictionary. Verify the actual H5AD `var` IDs before scoring. The four
  missing dictionary IDs remain documented in the feature overlap audit.
- Preserve observed guide pair counts for interpretation and subgroup
  checks. The 104 tasks with at least two observed pairs in both contexts
  form a diagnostic subgroup. Their guide-specific cell support and response
  concordance still require verification.

## Predictions and decision rule

Carry forward the Step 29 fits on the full internal panel. The paired
prediction uses the four-source mean for the same perturbation and the
common prediction uses the same vector for every target. Both raw vectors
are sliced to the candidate axes. Multiplying them by the archived scalars
gives the two calibrated baselines. The coefficients are `0.6019758092849216`
and `0.7550802413023521`, respectively. Include the zero baseline.

The primary risk score is `source_held_out_excess_mse` from Step 29. Report
the three legacy scores in the same direction as in the development run.
The primary policy retains the lowest `ceil(0.8 * N)` scores with the
unchanged Step 29 tie keys. It uses paired calibrated predictions for those
tasks and zero for the remainder. For 1,088 tasks this retains 871. The
external task counts only set cohort eligibility. The model and scores
retain their original values when projected to the candidate axes.

## Requirements before evaluating target expression

1. Confirm unique expression IDs and reconcile full H5AD `.obs` barcodes,
   sample and guide labels with the pinned guide CSVs. Verify source file
   version and checksum. Document any discrepancy before opening expression
   values. Use GEMs for control matching and report the status of independent
   biological replicate provenance.
2. Verify the candidate axes against that metadata audit. Fix the actual
   axes and interpretation rules, and hash the prediction and evaluation
   artifacts before computing or viewing target responses. If the axes
   differ, record the metadata reason and create a new numbered bundle
   before opening expression values.
3. Use mean `log1p(10000 * raw_count / original_total_UMI)` by task, minus the
   weighted mean of its GEM controls. Verify the count layer and original
   denominator. Use GEMs with at least 30 controls as primary and repeat
   with all GEMs with a control on the same tasks if applicable. Use the
   HCT116 and HEK293T barcode/guide joins to reconstruct actual cell support.

Evaluation will compare per-task and aggregate error for all five
baselines, including the fraction with `D = MSE - MSE_zero > 0`. Report
Pearson, Spearman and cosine on defined rows and retrieval on a common
eligible candidate set with its actual chance levels. Use label permutations
with a fixed seed, split-half reproducibility with matched controls, and
cell-support and perturbation-strength strata. Report AUROC for `D > 0`,
average precision with its prevalence reference, retained-task error and
the all-task zero-fallback policy error. Compare the policy with always
predicting, always using zero and random selection at the same coverage.
Report split-half and batch sensitivity with the unit of resampling and
the dependence from shared controls. The final evaluation implementation,
metadata confirmation and protocol tag are pending.

## Run the candidate projection

```bash
cd "$HOME/VC_reliability_repo" &&
git pull --ff-only origin dev &&
mkdir -p logs/v01_xatlas_prediction_candidate &&
sbatch slurm/step_30_prepare_xatlas_external.slurm
```

This job uses account `zhenwei01_bio214`, partition `cpu8358`, QOS `26cores`,
1 CPU and 4 GB of memory with a 10 minute limit. Expect under one minute
after scheduling. It reads the small local Step 29 H5 and metadata tables.

Check `sacct`, then the files from the submitted job:

```bash
vc_job_id=YOUR_JOB_ID
vc_out="results/v01_xatlas_prediction_candidate/job_$vc_job_id"
sacct -X -j "$vc_job_id" --format=JobID,State,ExitCode,Elapsed
sha256sum --quiet -c "$vc_out/SHA256SUMS" && cat "$vc_out/summary.json"
```

Keep `predictions.h5` on the HPC. Preserve `task_eligibility.tsv`,
`feature_axis.tsv`, `summary.json`, `parameters.json`, the code commit,
logs and SHA256 entries in the audit trail. Archive the small outputs on
`dev` after the job completes.

## Verified HPC run 3239864

The user reported Slurm `COMPLETED 0:0` in eight seconds, `SHA256_OK`,
and `XATLAS_PREDICTION_CANDIDATE_COMPLETE 3239864` for commit `5c7492c`.

| Project result | Value |
| --- | ---: |
| Candidate perturbations | 1,088 |
| Candidate expression genes | 6,638 |
| Paired calibration coefficient | 0.6019758092849216 |
| Common calibration coefficient | 0.7550802413023521 |
| Retained predictions | 871 |
| Retained fraction | 80.0551% |

The run status is `XATLAS_PREDICTION_CANDIDATE`, with
`atlas_axis_confirmed: false`. The counts describe the projection onto
the candidate metadata axes. The exact H5AD gene-axis audit follows in
Step 31. The full cell barcode and guide checks follow at intake.

## Archive the actual HPC outputs

```bash
(
  set -euo pipefail
  cd "$HOME/VC_reliability_repo"
  git pull --ff-only origin dev
  vc_out=results/v01_xatlas_prediction_candidate/job_3239864
  vc_manifest=manifests/v01_xatlas_prediction_candidate_3239864.sha256
  vc_files=("$vc_out"/{task_eligibility.tsv,feature_axis.tsv,summary.json,parameters.json,code_commit.txt,SHA256SUMS}
            logs/v01_xatlas_prediction_candidate/step30_3239864.out)
  sha256sum --quiet -c "$vc_out/SHA256SUMS"
  sha256sum -- "${vc_files[@]}" > "$vc_manifest"
  sha256sum --quiet -c "$vc_manifest"
  git add -- "${vc_files[@]}" "$vc_manifest"
  git -c core.whitespace=cr-at-eol diff --cached --check -- "${vc_files[@]}" "$vc_manifest"
  git commit --only -m "results: archive X-Atlas prediction candidate job 3239864" -- \
    "${vc_files[@]}" "$vc_manifest"
  git push origin dev
)
```

The public manifest covers the archived reports and the original run
checksum list. The original list also records the hash of the local H5.

[Step 31](step_31_xatlas_gene_axis.md)
