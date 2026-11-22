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
