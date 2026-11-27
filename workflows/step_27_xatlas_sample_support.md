# Step 27: inspect X-Atlas task support across GEM batches

The v05 `.obs` probe and both guide metadata audits have passed. The
author's dataset card defines `sample` as a GEM batch. A GEM is useful for
control matching, but is not established as an independent biological
replicate. This step uses the already computed guide tables; it reads no
external expression matrix or target response.

Run from the HPC repository root:

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
module load anaconda3/2022.10
sha256sum --status -c qc/v01_xatlas_hct116_observed/SHA256SUMS
sha256sum --status -c qc/v01_xatlas_hek293t_observed/SHA256SUMS

python3 scripts/external_metadata/10_xatlas_sample_support.py \
  --hct116-dir qc/v01_xatlas_hct116_observed \
  --hek293t-dir qc/v01_xatlas_hek293t_observed \
  --output-dir qc/v01_xatlas_sample_support

sha256sum -- qc/v01_xatlas_sample_support/{task_sample_support.tsv,summary.json,parameters.json} \
  > qc/v01_xatlas_sample_support/SHA256SUMS
sha256sum -c qc/v01_xatlas_sample_support/SHA256SUMS
cat qc/v01_xatlas_sample_support/summary.json
```

The output directory must not exist before the run. This small metadata
calculation uses Python 3.9 and should finish on the login node without a
Slurm allocation. It does not reopen the multi-gigabyte guide CSV files.
If an input or reconciliation check fails, inspect that error before
choosing a new versioned output directory.

`task_sample_support.tsv` reports, for each target and cell line, the
observed cells, guide pairs, GEMs, fraction in the largest GEM, effective
GEM count `1 / sum((cells_gem / cells_task)^2)`, cells in GEMs with fewer
than 30 controls, and the minimum control count in any used GEM.
`summary.json` counts targets meeting the existing 30-cell diagnostic
cutoff, and those with at least two GEMs or guide pairs. The script checks
these counts against the previous audits, checks control counts for every
used GEM, and saves input SHA256 values in `parameters.json`.

These are **metadata diagnostics**, not a new eligibility rule. Multiple
GEMs or guide pairs do not establish independent biological replication,
guide concordance, knockdown efficiency, or predictive accuracy. The next
gate must document experimental replicate provenance and verify full
barcode and label membership between expression and guide files. Freeze
the external prediction and evaluation protocol before reading target
expression responses.

Source: [Xaira X-Atlas/Orion dataset card](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion)
and [Figshare v3 record](https://doi.org/10.25452/figshare.plus.29190726.v3).
