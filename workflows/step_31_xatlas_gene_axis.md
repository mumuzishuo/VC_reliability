# Step 31 · Match the actual X-Atlas gene axis

Use the verified feature table from Step 30 (`job_3239864`) and the
pinned [Figshare v3 record](https://doi.org/10.25452/figshare.plus.29190726.v3).
This audit visits the `var` index and available identifier fields by HTTPS
ranges. It does not decode `X`, per-gene expression summaries, or target
responses. The archived Step 26 Range reader limits each file to 8 MiB
and 256 requests. The full expression file is not downloaded in this step.

## Run HCT116 first

After [archiving Step 30](step_30_prepare_xatlas_external.md#archive-the-actual-hpc-outputs),
run from the HPC repository root:

```bash
cd "$HOME/VC_reliability_repo" &&
git pull --ff-only origin dev &&
mkdir -p logs/v01_xatlas_gene_axis &&
sbatch slurm/step_31_xatlas_gene_axis.slurm
```

The job uses account `zhenwei01_bio214`, `cpu8358`, QOS `26cores`,
one CPU and 2 GB. Allow several minutes after scheduling, with a 30
minute time limit. Inspect the run output if the remote server reaches a
limit or returns an HTTP error. With the same script, the HEK293T job
can be submitted with `sbatch slurm/step_31_xatlas_gene_axis.slurm HEK293T`
after the HCT116 result is checked.

## Verify and interpret the output

```bash
vc_job_id=YOUR_JOB_ID
vc_out="qc/v01_xatlas_gene_axis/job_$vc_job_id"
sacct -X -j "$vc_job_id" --format=JobID,State,ExitCode,Elapsed
sha256sum --quiet -c "$vc_out/SHA256SUMS" && cat "$vc_out/summary.json"
```

`EXACT_VAR_AXIS_MATCH` requires 6,638 unique exact Ensembl IDs with no
disagreement between the examined identifier fields. The mapping table
gives zero-based matrix column indices in Step 30 feature order.
`VAR_AXIS_REVIEW_NEEDED` records missing, duplicate or conflicting IDs.
Review the `field_diagnostics` and TSV before changing the panel. A
symbol-only `var` index requires a documented crosswalk. The script never
accepts symbol-only or version-stripped matches automatically.

Preserve the JSON summary, mapping, code commit, logs and checksums. A
successful probe establishes feature correspondence for the pinned
release. The full-file MD5 check and one-to-one `.obs` barcode and guide
label reconciliation remain required before reading target responses.
Any later axis changes must be based on metadata and recorded before
target response evaluation. The final evaluation protocol tag remains
pending.
