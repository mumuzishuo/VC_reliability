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

## Verified project metadata result

The user ran the script after pulling commit `c3ae519` and reported that
`summary.json`, `parameters.json` and `task_sample_support.tsv` all
passed the output SHA256 check. The results are current-project calculations
from guide metadata, not author-reported performance results:

| Diagnostic | HCT116 | HEK293T | Both contexts |
| --- | ---: | ---: | ---: |
| Tasks with at least 30 guide-validated cells | 1,110 | 1,130 | 1,088 |
| Those present in at least two GEMs | 1,110 | 1,130 | 1,088 |
| Those with at least two observed guide **pairs** | 104 | 104 | 104 |

Thus 984 of the 1,088 cell-supported tasks lack a second observed guide
pair in both contexts. The two members of a single dual-guide pair are
not two independent guide-pair measurements. The GEM threshold is a
technical distribution diagnostic, not biological replication.

The aggregate JSON does not describe the largest-GEM fractions, effective
GEM counts or exposure to GEMs with fewer than 30 controls. Inspect the
per-task TSV after it is archived before drawing conclusions about those
quantities. Archive the three derived outputs, their local SHA256 list and
a portable manifest on `dev`; keep the guide CSVs outside Git. No target
expression response was accessed in this step.

## Distribution review from the archived task table

Input: `qc/v01_xatlas_sample_support/task_sample_support.tsv`, SHA256
`20c6d38894d11c88bb1fe6216ffed29ff10f2b11abe686558e4e1bf5f9ba3544`
(pinned in `manifests/v01_xatlas_sample_support_public.sha256`). The following
are current-project calculations restricted to tasks with at least 30
validated cells **in that context**:

| Diagnostic | HCT116 (1,110 tasks) | HEK293T (1,130 tasks) |
| --- | ---: | ---: |
| Minimum distinct GEMs for any qualifying task | 24 | 22 |
| Median distinct GEMs | 75 | 106 |
| Median fraction of task cells in its largest GEM | 3.70% | 3.57% |
| Maximum fraction of task cells in its largest GEM | 17.65% | 21.88% |
| Qualifying tasks exposed to a GEM with <30 controls | 0 | 15 |
| Total qualifying perturbation cells in such GEMs | 0 | 26 |
| Minimum controls in a used GEM | 868 | 12 |

All 15 affected HEK293T tasks are among the 1,088 qualifying in both
contexts. The largest low-control-GEM contribution to any one qualifying
HEK293T task is 3/121 cells (2.48%). These small exposures do not establish
how much the eventual response changes if those cells are removed. Keep
both the complete and the >=30-control-GEM response constructions as a
predeclared sensitivity comparison; do not choose after seeing target MSE.

A compact reproduction of the summary (Python 3.9, run from repository root):

```python
import csv
from statistics import median

path = "qc/v01_xatlas_sample_support/task_sample_support.tsv"
with open(path, encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

for line in ("HCT116", "HEK293T"):
    selected = [r for r in rows if r["cell_line"] == line and int(r["cells"]) >= 30]
    gems = [int(r["observed_gems"]) for r in selected]
    largest = [float(r["largest_gem_cell_fraction"]) for r in selected]
    low = [int(r["cells_in_gems_lt30_controls"]) for r in selected]
    controls = [int(r["minimum_controls_in_used_gem"]) for r in selected]
    print(line, len(selected), min(gems), median(gems),
          median(largest), max(largest), sum(x > 0 for x in low),
          sum(low), min(controls))
```

The median and range describe technical capture distribution, not a count
of independent experimental repeats. The released `.obs` schema has a
`sample` field, but the available metadata do not yet establish a mapping
from its GEM labels to independent biological experiments. Full barcode and
label matching, replicate provenance, and the locked evaluation protocol
remain necessary before reading external response data.
