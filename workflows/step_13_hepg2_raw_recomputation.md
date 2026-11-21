# Step 13: direct HepG2 raw-count audit

Purpose: close the direct raw-count verification gap for the **target** HepG2
matrix. The older 35-task direct audit covered only three **source** contexts.
This is an arithmetic quality check, not a model evaluation or external test.

## Inputs and selection

- Legacy `processed/v01/manifest.json`, `vc_preprocess.py`, raw HepG2 `.h5ad`
  and `processed/v01/evaluation_only/HepG2_responses.h5`.
- `prepare` inspects only manifest, response axes, perturbation cell/GEM counts,
  and control cell counts. It does not read raw `X` or stored `delta` values.
- It selects five deterministic cell-support quantiles, the task concentrated
  in one GEM, and the top two exposed tasks for **each** GEM with fewer than
  30 controls. All selection reasons and input hashes are written before `audit`.
- `audit` independently normalizes raw counts in float64 using the original
  UMI_count, accumulates perturbed and GEM control means, and computes the
  matched mixture. It compares control means, perturbed means, matched controls,
  delta and low-control fractions with saved results at absolute tolerance
  `2e-6`. The legacy metadata reader decodes categorical labels; its SHA256
  must match the original preprocessing manifest.

## Run once on HPC

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
mkdir -p logs/v01_hepg2_raw_recomputation
bash -n slurm/step_13_recompute_hepg2_raw_counts.slurm
sbatch --account="$VC_ACCOUNT" --partition="$VC_PARTITION" --qos="$VC_QOS" \
  slurm/step_13_recompute_hepg2_raw_counts.slurm
```

Slurm requests one job, two CPUs, 8 GB memory and one hour of wall time.
Set `VC_ACCOUNT`, `VC_PARTITION`, and `VC_QOS` locally to the allocation
confirmed for your HPC account before submission. Keep site-specific account
details out of the public repository.
The job refuses to replace `qc/v01_hepg2_raw_recomputation` if it already
exists. Never remove that directory just to rerun; inspect a failure first and
change to a new versioned directory in the code and Slurm script if warranted.

After the job finishes, replace `JOB_ID` below with the submitted job ID:

```bash
sacct -X -j JOB_ID --format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem
tail -n 30 logs/v01_hepg2_raw_recomputation/step13_JOB_ID.out
tail -n 30 logs/v01_hepg2_raw_recomputation/step13_JOB_ID.err
sha256sum --status -c qc/v01_hepg2_raw_recomputation/SHA256SUMS
python3 -m json.tool qc/v01_hepg2_raw_recomputation/audit_summary.json
```

The job saves the frozen `selection.tsv` and `selection.json`, `per_task.tsv`,
`audit_summary.json`, relative-path SHA256 input list, script and job SHA256,
Python version, commit and Slurm job metadata. Scheduler success alone is not
evidence of an arithmetic pass; inspect `audit_summary.json: status` and the
per-task errors. If the audit fails, preserve logs and hashes for diagnosis.

## Interpretation

A pass confirms this deterministic subset of HepG2 responses can be recovered
from the raw counts within tolerance. It does not check every task, establish
independent biological replicates, prove all controls unbiased, or validate a
model's cross-context risk score. Selection uses target cell/GEM metadata for
**QC only**; neither it nor target responses may train a prospective risk rule.
