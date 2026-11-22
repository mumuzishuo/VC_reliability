# Step 25: count HEK293T guide support and compare both contexts

Start after the HCT116 guide audit passes its SHA256 check. This step uses
HEK293T guide calls, the same guide library and task panel, and the verified
HCT116 task table. The HEK293T file has ID 59490734, 90,886,568 bytes, and
MD5 `f2106fd8fdca7350160e9e5a91c67be8`. Its source is the
[Figshare v3 record](https://doi.org/10.25452/figshare.plus.29190726.v3).

Run these commands on the HPC login node:

```bash
cd "$HOME/VC_reliability_repo" &&
git pull --ff-only origin dev &&
bash scripts/external_metadata/06_xatlas_hek293t_intake.sh &&
sbatch slurm/step_25_xatlas_hek293t_observed.slurm
```

The intake runs over HTTPS on the login node. It verifies the byte count,
official MD5 and gzip integrity before publishing the downloaded file.
It stores the requested URL and check results in
`logs/v01_xatlas_hek293t_observed/intake_*.log`, and writes the file's SHA256
to `logs/v01_xatlas_hek293t_observed/guide_calls.sha256`. If a partial
file is left after an error, inspect its intake log before retrying.

The job runs on account `zhenwei01_bio214`, partition `cpu8358`, QOS
`26cores`, with one CPU, 6 GB memory, and a 20 minute limit.
The scan is expected to finish within a few minutes of starting.

## Compare support using the same rules

The job checks the input hashes, uses the same guide audit as Step 24,
and joins the task tables by exact Ensembl ID. It counts tasks with at least
30 guide validated cells in HCT116, HEK293T, and both contexts. It repeats
that calculation for cells in GEMs with any non-targeting controls and in
GEMs with at least 30 controls.

These are guide data diagnostics. A later gate must verify expression cell
membership by barcode and check the experimental replicate structure.
Target expression stays sealed until the external protocol is frozen.
The author defines samples as GEMs and describes the cell barcode suffix
and guide filter in the
[X-Atlas dataset card](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion).

## Verify and save

The Slurm job ends with `XATLAS_PAIRED_GUIDE_AUDIT_OK` in its stdout log.
After the job completes, run this from the repository root:

```bash
sha256sum --status -c qc/v01_xatlas_hek293t_observed/SHA256SUMS &&
cat qc/v01_xatlas_hek293t_observed/{summary,paired_support_summary}.json
```

`qc/v01_xatlas_hek293t_observed/` contains the HEK293T task, GEM, guide pair
and control tables, plus `paired_task_support.tsv` and
`paired_support_summary.json`. Keep the parameters, code commit,
SHA256 manifests and Slurm logs with the QC outputs. The two source
guide CSV files stay on the HPC and are excluded from Git.

## Verified HPC result

On 2026-11-22, the intake at commit `0c5c1c996135f4a1340e1d66f156aabe67143c0f`
downloaded 90,886,568 bytes from Figshare file 59490734.
The byte size, published MD5, and gzip integrity checks passed. The verified
SHA256 is pinned in `manifests/v01_xatlas_hek293t_guide_calls.sha256`.
The user verified the Step 25 output manifest for job `3220633` on the HPC
and reported the following summary. HCT116 values come from Step 24.

| Guide CSV metric | HCT116 | HEK293T |
| --- | ---: | ---: |
| Unique cell barcodes | 3,409,169 | 4,534,299 |
| Excluded guide CSV records | 0 | 0 |
| GEMs | 109 | 223 |
| Non-targeting cells | 165,777 | 218,838 |
| Minimum non-targeting cells in one GEM | 868 | 12 |
| GEMs with fewer than 30 non-targeting cells | 0 | 1 |
| Targeting cells | 3,243,392 | 4,315,461 |
| Observed target Ensembl IDs | 18,293 | 18,311 |
| Main panel targets with observed cells | 1,230 | 1,229 |
| Cells assigned to main panel targets | 179,096 | 237,422 |
| Main panel targets with at least 30 guide-validated cells | 1,110 | 1,130 |

Across the 1,264 main panel targets, **1,088 (86.1%)** reach the cell cutoff
in both contexts. There are 22 that reach it in HCT116 only, 42 in HEK293T
only, and 112 in neither. Restricting the counts to cells from GEMs with any
controls, or at least 30 controls, yields the same qualifying task counts and
intersection. The HEK293T GEM with 12 controls does not change which tasks
reach the cutoff. Its effect on expression estimates remains untested.

The 7,943,468 total barcodes are project counts from the two guide CSVs.
The author-reported total is approximately 8 million cells. The archived
Figshare description uses the same guide filter for its guide CSVs and H5ADs.
Exact barcode membership in the expression cohort still needs to be checked.
Both GEM counts are capture batch counts. Biological replicate provenance
remains open. [Step 26](step_26_xatlas_expression_index_probe.md)
