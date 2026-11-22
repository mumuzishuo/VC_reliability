# Step 26: probe the X-Atlas expression cell index

[Steps 24 and 25](step_25_xatlas_hek293t_observed.md) establish 1,088 shared
tasks with guide cell support. Before counting the available expression
cells, check how the exact Figshare v3 H5ADs encode the cell index and GEM.
The published HCT116 and HEK293T H5ADs have file IDs 55021257 and 55074802,
respectively. The two files total 559.52 GB (521.09 GiB).
[Figshare v3](https://doi.org/10.25452/figshare.plus.29190726.v3)

The audit sends HTTPS range requests. It reads the length of the H5AD
`obs` index, the first 16 barcodes and sample labels, and the available
`obs` field names. For each cell line, it compares the preview with the first
4,096 guide CSV records and the verified guide CSV total.

The reader accepts only a matching HTTP 206 response and caps bytes read
at 8 MiB per H5AD. The probe uses the Anaconda Python 3.9 and h5py environment.
Counts, response matrices, and individual barcodes are excluded from the output.
The probe does not access the expression matrix.

```bash
cd "$HOME/VC_reliability_repo" &&
git pull --ff-only origin dev &&
mkdir -p logs/v01_xatlas_expression_index_probe &&
sbatch slurm/step_26_xatlas_expression_index_probe.slurm
```

Slurm uses account `zhenwei01_bio214`, partition `cpu8358`, QOS `26cores`,
one CPU, 2 GB RAM, and a 15 minute limit. The probe is expected to finish
in a few minutes once scheduled. It requires HTTPS access from the compute node.

Once the job prints `XATLAS_OBS_INDEX_PROBE_OK`, verify and read its summary:

```bash
sha256sum --status -c qc/v01_xatlas_expression_index_probe/SHA256SUMS &&
cat qc/v01_xatlas_expression_index_probe/summary.json
```

Save the Slurm logs, summary, code commit and SHA256 manifest. If a range
request fails, use the Slurm error log to diagnose the transport. On success,
use the confirmed index and sample fields to plan the full cohort join.
Identical row totals and matching preview barcodes still require an exact
full join. Distinct order can explain a preview mismatch. A GEM or sample
label gives a capture batch and cannot establish independent replicates.
Freeze the model, calibration and risk rules before reading target responses.
