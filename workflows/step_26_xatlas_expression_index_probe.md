# Step 26: probe the X-Atlas expression cell index

**Current status:** Slurm job `3220676` stopped while opening the HCT116 H5AD.
Run the [login-node range check](#login-node-range-check) before retrying the
index probe.

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

## Slurm attempt 3220676

The job ran from commit `e699a01be51a43b919e4e91c174cc680c11dd155`.
Slurm reported `FAILED`, exit code `1:0`, after 23 seconds. The last stdout
line was `Checking HCT116 .obs metadata from Figshare file 55021257`.
While `h5py` was opening the remote file, `urllib` followed an HTTP 302
redirect and the subsequent HTTPS connection returned `Errno 104`
(connection reset by peer). The probe did not produce a verified index
summary or a SHA256 manifest. The available log does not identify the
component that reset the connection. Keep the original stdout and stderr
under `logs/v01_xatlas_expression_index_probe/`.

## Login node range check

Run this on the login node. It uses the same curl environment as the guide
CSV intake and requests bytes 0–7 from each pinned Figshare v3 H5AD. Curl
has a 25 second timeout per file, a 2,048 byte size limit when the length
is known, and a 1 KiB/s transfer rate limit. Response bodies are discarded.

```bash
git pull --ff-only origin dev &&
bash scripts/external_metadata/09_probe_xatlas_range.sh
```

`RANGE_OK` requires curl exit 0, HTTP 206, 8 received bytes, and a
`Content-Range` that exactly matches the requested range and recorded H5AD
size. `RANGE_UNCONFIRMED` retains the HTTP and curl diagnostics. This checks
transfer support. Cohort membership still requires the H5AD index probe.

Logs and their SHA256 manifests go to
`logs/v02_xatlas_expression_index_probe/`. The log contains response statuses
and the returned range. Run the full index probe after a route is verified
and the failed run has been reviewed.
