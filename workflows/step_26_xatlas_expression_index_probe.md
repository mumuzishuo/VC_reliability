# Step 26: probe the X-Atlas expression cell index

**Current status:** Slurm v05 completed its bounded `.obs` index probe for both
cell lines. H5AD row counts and the first 16 barcode positions match the guide
CSVs; full barcode membership and independent replicate structure remain open.
The v05 outputs are under `qc/v05_xatlas_expression_index_probe/job_3237540/`.

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

The Anaconda Python 3.9 and h5py job uses curl for each HTTPS range.
It first reads eight bytes per file to verify the route on the compute node.
It accepts only matching HTTP 206 responses and caps the experiment at
8 MiB of accepted range payload per H5AD.
Counts, response matrices, and individual barcodes are excluded from the output.
The probe does not access the expression matrix.

```bash
cd "$HOME/VC_reliability_repo" &&
git pull --ff-only origin dev &&
mkdir -p logs/v05_xatlas_expression_index_probe &&
sbatch slurm/step_26_xatlas_expression_index_probe.slurm
```

Slurm uses account `zhenwei01_bio214`, partition `cpu8358`, QOS `26cores`,
one CPU, 2 GB RAM, and a 15 minute limit. The Anaconda module supplies
curl 7.84 or newer for the response header fields. The probe is expected to finish
in a few minutes once scheduled. It requires HTTPS access from the compute node.

Once the job prints `XATLAS_OBS_INDEX_PROBE_OK`, verify and read its summary:

```bash
vc_job=1234567  # Replace with the job ID printed by sbatch.
vc_out="qc/v05_xatlas_expression_index_probe/job_${vc_job}"
sha256sum --status -c "$vc_out/SHA256SUMS" &&
cat "$vc_out/summary.json"
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
and the returned range.

### Verified login node result

At 2026-11-22 14:12:24 UTC, the login-node check at commit
`ace6e89f10723ad076b05964bf8a1fe88f0514cb` observed two HTTP 302 redirects,
a final HTTP 206, and eight bytes for each file. Curl 7.84.0 exited with
status 0 for both files.

| File | Final Content-Range |
| --- | --- |
| HCT116 | `bytes 0-7/209354246272` |
| HEK293T | `bytes 0-7/350164035901` |

The source is the HPC log
`logs/v02_xatlas_expression_index_probe/range_GqsiQv.log`.
Its checksum manifest was saved alongside it. This proves a valid partial
transfer from the login node for each file. A successful v05 job confirms the first
read on its compute node with `COMPUTE_RANGE_OK` and reports the H5AD
index and batch checks in the summary. If it stops before that marker,
inspect its stderr for the requested range, curl exit code, and HTTP
status. This test does not isolate the node from the HTTP client in the
earlier failure.

## Slurm attempt 3237360

The v03 job exited `FAILED 1:0` after five seconds. Input-hash checks
completed, then Python printed the curl version and stopped before reading
any H5AD bytes: `qc/v03_xatlas_expression_index_probe` already existed and
`mkdir(..., exist_ok=False)` raised `FileExistsError`. The user checked
that the directory contained no files or SHA256 manifest. The v03 logs and
empty directory remain available for audit. The v04 job uses a separate
`job_<Slurm ID>` directory for each attempt.

## Slurm attempt 3237464

The v04 job exited `FAILED 1:0` after 16 seconds. It printed
`COMPUTE_RANGE_OK HCT116 bytes=8`, then the subsequent HDF5 `.obs` read
received HTTP 403 on bytes 209354245576–209354245751 (176 bytes).
The v04 code reused the final redirected URL from the first request. The
observed sequence is consistent with a redirected URL refusing reuse; it does
not prove that is the sole cause. In v05 each range starts at the original
Figshare API URL and follows a fresh HTTPS redirect. Keep the v04 stdout,
stderr and job-specific partial directory; report whether v05 reaches the
HCT116 and HEK293T summaries before drawing a conclusion about transport.

Before publishing v03, the range reader was tested with a local HTTPS
server for two redirects, exact byte reads, HTTP errors, oversized bodies,
inconsistent ranges, and changed ETags. The Slurm script passed Bash syntax
checking. The live `.obs` read still needs the HPC run.

## Verified v05 metadata probe (job 3237540)

The user reported Slurm `COMPLETED 0:0` in 3 min 21 sec, an empty stderr,
`XATLAS_OBS_INDEX_PROBE_OK 3237540`, and `sha256sum -c` passing for the
summary, code commit, both scripts and all referenced source QC files.
These are current-project HPC calculations, not author-reported expression
performance.

| Metadata diagnostic | HCT116 | HEK293T |
| --- | ---: | ---: |
| Guide CSV unique cells | 3,409,169 | 4,534,299 |
| H5AD `.obs` rows | 3,409,169 | 4,534,299 |
| First 16 H5AD barcodes matched among first 4,096 guide rows | 16 | 16 |
| Same barcode order in first 16 | Yes | Yes |
| First 16 `sample` values match barcode suffix | 16 | 16 |
| HTTPS range requests | 36 | 36 |
| Remote bytes accepted | 146,448 | 115,520 |

Both `.obs` tables expose `gene_target`, `guide_target`, `num_features`,
`pass_guide_filter`, `sample`, `total_counts`, `total_counts_mt`,
`n_genes_by_counts` and `pct_counts_mt`. Only the index preview and
the schema were read; `X` and target perturbation responses were not read.
Equal row totals and 16 ordered matches strongly support consistent export
but do not prove full barcode membership, guide agreement for all cells,
independent biological replicates, or external predictive performance.

Next: archive the exact HPC summary and a portable SHA256 manifest in `dev`.
Predeclare how to handle unmatched barcodes and missing/low-support GEM
controls, then perform a full metadata barcode join as the first operation
when the complete expression cohort is ingested. Check replicate provenance
separately; GEMs are capture batches and must not be counted as independent
biological repeats. Freeze prediction, calibration and risk rules before
using any external target response for evaluation.
