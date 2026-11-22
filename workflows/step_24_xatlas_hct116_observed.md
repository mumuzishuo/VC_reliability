# Step 24: count observed HCT116 guide support and GEM controls

Run this on the HPC after the [HCT116 intake](step_23_xatlas_hct116_guide_count_header.md).
The job verifies the SHA256 of the downloaded guide calls, the released guide
library, and the existing 1,264-task panel. The CSV contains guide calls and
guide UMI counts. The job reads those three metadata files.

```bash
cd "$HOME/VC_reliability_repo" &&
git pull --ff-only origin dev &&
mkdir -p logs/v01_xatlas_hct116_observed &&
sbatch slurm/step_24_xatlas_hct116_observed.slurm
```

Slurm uses account `zhenwei01_bio214`, partition `cpu8358`, QOS `26cores`,
one CPU, 4 GB memory, and a 20-minute wall time limit. The scan is expected
to take a few minutes after the job starts. The actual time will be available
from `sacct`.

## How the counts are made

1. A full cell barcode identifies one cell. The script stops if it encounters
   a repeated barcode. The `HCT116_…` suffix provides the GEM identifier.
2. A counted cell passes `pass_guide_filter`, has two guide features, and has
   an exact guide pair in the released library. Pair order can be reversed.
   The script checks that the cell gene label agrees with the library label.
   The summary counts exclusions by the first reason that applies.
3. The guide library gives the Ensembl ID for matching to the main panel.
   It also identifies non-targeting pairs. Count the observed cells carrying
   those pairs as controls in each GEM.
4. Report tasks with at least 30 guide-validated cells. Also count tasks with
   at least 30 cells after restricting to GEMs with any non-targeting controls
   or with at least 30 controls. These are support diagnostics. Freeze the
   final analysis thresholds before reading target expression responses.

The [author dataset card](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion)
defines the barcode suffix as the sample, identifies the sample as GEM, and
describes the guide filter. The
[Figshare record](https://doi.org/10.25452/figshare.plus.29190726.v3) provides
the guide calls and library. These counts describe the guide CSV cohort.
Membership in the expression cohort requires a barcode join. A GEM identifies
a capture batch. Independent biological replication needs experimental metadata.

## Verify and record

Wait for `HCT116_GUIDE_AUDIT_OK` in the Slurm stdout log. Then run from the
repository root:

```bash
sha256sum -c qc/v01_xatlas_hct116_observed/SHA256SUMS &&
cat qc/v01_xatlas_hct116_observed/summary.json
```

| File in `qc/v01_xatlas_hct116_observed/` | Contents |
| --- | --- |
| `summary.json` | Unique cells, exclusion reasons, GEMs, controls, observed panel coverage and threshold counts |
| `gem_control_counts.tsv` | Non-targeting cell and observed control pair counts by GEM |
| `task_support.tsv` | Task-level support, guide pairs and cells in GEMs with controls |
| `task_gem_counts.tsv` | Per task and GEM counts with the matching control count |
| `task_guide_counts.tsv` | Observed guide pair support by panel target |
| `parameters.json`, `code_commit.txt`, `SHA256SUMS` | Rules, code version and checksums |

The input hashes are in `manifests/v01_xatlas_hct116_observed_inputs.sha256`.
Keep the Slurm stdout/stderr logs with the run. After verification, sync the
aggregate QC outputs and run records to GitHub. The source guide files remain
on the HPC.

## Verified HPC result

Slurm job `3220603` used commit `72ccd2595b081344d11f6ed9bcd8ee25f3e27a2a`.
The user ran `sha256sum -c qc/v01_xatlas_hct116_observed/SHA256SUMS`
on the HPC and reported that every listed file passed. The counts below are
from that job's `summary.json`.

| HCT116 guide metadata measure | Project result |
| --- | ---: |
| Unique full barcodes | 3,409,169 |
| Excluded cells | 0 |
| GEMs | 109 |
| Non-targeting cells | 165,777 |
| Minimum non-targeting cells in one GEM | 868 |
| GEMs with zero or fewer than 30 non-targeting cells | 0 |
| Targeting cells | 3,243,392 |
| Unique observed target Ensembl IDs | 18,293 |
| Main panel targets in the design | 1,231 of 1,264 |
| Main panel targets with observed cells | 1,230 of 1,264 |
| Cells assigned to main panel targets | 179,096 |
| Main panel targets with at least 30 cells | 1,110 of 1,264 (87.8%) |
| Main panel targets with at least 30 cells from GEMs with controls | 1,110 |
| Main panel targets with at least 30 cells from GEMs with at least 30 controls | 1,110 |

This passes the guide data support check for 1,110 HCT116 candidate tasks.
The next step repeats these counts in HEK293T and computes their intersection.
The check of expression cohort membership and independent biological
replicates remains open. [Step 25](step_25_xatlas_hek293t_observed.md)
