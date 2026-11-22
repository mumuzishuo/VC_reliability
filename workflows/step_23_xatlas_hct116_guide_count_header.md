# Step 23: inspect X-Atlas HCT116 per-cell guide-count metadata

Goal: determine the schema and join keys of the **guide-count** file before
estimating observed per-perturbation support or non-targeting controls.

The archived Figshare+ record 29190726, version 3, lists:

| Field | Pinned value |
| --- | --- |
| File | HCT116_filtered_guide_calls_per_cell.csv.gz |
| File ID | 59490731 |
| Bytes | 69606079 |
| MD5 | 2f6a6c751e14c600dcc833964e90b1b2 |
| Record DOI | 10.25452/figshare.plus.29190726.v3 |

These values come from the archived database record and
`manifests/v01_xatlas_guide_count_source.tsv`.
The companion HEK293T file is listed there but is not downloaded in this step.
The author describes these files as sgRNA count data. Their cell barcode key,
sample/GEM join, control encoding and relationship to filtered expression cells
have not been verified.

## Transport diagnosis reported from the HPC on 2026-11-22

| Request | User-reported result |
| --- | --- |
| Original GET through ndownloader.figshare.com | curl 35, SSL connection reset |
| Range GET through figshare.com/ndownloader | HTTP 202, 0 bytes, 0 redirects, curl exit 0 |
| Range GET through api.figshare.com/v2/file/download/59490731 | HTTP 206, 1024 bytes, 2 redirects, curl exit 0 |

These are results pasted by the user, not a new test run elsewhere.
The last request establishes that the API route returned a byte range on that
attempt. It does not establish full-file integrity or identify the cause of
the other responses. Figshare documents the public file download endpoint at
<https://docs.figshare.com/v2/>.
The original published download URL is retained in the source manifest;
the API request URL used for this attempt is recorded in the new parameters file.

## Download, verify, and preview

The script downloads the pinned HCT116 guide-count file, verifies its exact
byte size and official MD5, checks gzip integrity, and previews the CSV header
and three rows. It reuses the earlier CSV preview logic.
This is file intake on the login node using Python 3.9; no Slurm job is
submitted. The expression H5AD is not accessed.

```bash
(
  set -e -o pipefail
  cd "$HOME/VC_reliability_repo"
  git pull --ff-only origin dev
  module load anaconda3/2022.10
  bash -n scripts/external_metadata/04_xatlas_hct116_intake.sh
  bash scripts/external_metadata/04_xatlas_hct116_intake.sh
)
```

The new run directory is `logs/v02_xatlas_guide_calls`.
An existing run directory causes a stop, preserving the previous attempt.
An existing final input is verified rather than downloaded again. An existing
partial transfer also causes a stop for diagnosis. A newly downloaded file
is published from its Git-ignored `.partial` path only after verification
and preview succeed; an existing final path is never replaced.

The transfer prints a start message, then may be quiet until completion.
Each curl attempt has a 900-second limit and at most two retries.
A checksum, size, gzip, or parsing error stops execution and preserves the
candidate file. Do not delete evidence or rerun into an existing output
directory without first diagnosing the failure.

`gzip -t` scans the complete compressed stream for integrity without saving
an uncompressed copy. CSV parsing itself reads only the header and three rows.

## Outputs and interpretation

- `logs/v02_xatlas_guide_calls/run.log`: UTC start time, code commit,
  Python/curl versions, transfer outcome, integrity results and CSV preview.
- `logs/v02_xatlas_guide_calls/parameters.tsv`: pinned file, source record,
  API URL, transfer limits and preview parameters.
- `logs/v02_xatlas_guide_calls/SHA256SUMS`: hashes of the final input, parameters,
  run log, intake script and original source manifest, using repository-relative
  paths.
- `metadata/v01_external_catalog/HCT116_filtered_guide_calls_per_cell.csv.gz`:
  verified source file, kept Git-ignored.

Success ends with `HCT116_METADATA_INTAKE_OK`. To verify the saved manifest,
run from the repository root:

```bash
sha256sum --status -c logs/v02_xatlas_guide_calls/SHA256SUMS
```

The script's shell syntax was checked before publication. The full HPC intake
has not yet been confirmed. Review the resulting schema before writing the
streaming cell-count audit or committing the intake logs and manifest.

A guide count is not a unique cell, and a non-targeting guide is not proof of
sufficient matched control cells. This step cannot establish eligible task
counts, independent biological replicates, or prediction performance.
Keep external target responses unseen.
