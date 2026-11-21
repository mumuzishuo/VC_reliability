# Step 15: Jiang ST1f design label overlap

The input is the publisher's Supplementary Tables XLSX (local to the HPC),
SHA256 `2037dcde3ce88656e334f9c498c3f2571c05ae08727ea7fea69b1c7e29f070af`,
2,962,962 bytes. The repository ignores the publisher workbook. The script
reads only the ST1f *design-target* sheet, not DE, MultiCCA, clustering,
expression or observed target responses. It avoids the reported 1,048,576
rows in ST1a. A match means exact equality of target symbols after trimming
whitespace; aliases are not assumed. The target panel was derived from four
existing contexts with at least 30 cells per context.

Run on the login node with Python 3.9; this scans a 3 MB workbook and the
existing panel only, so no Slurm allocation is needed:

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
module load anaconda3/2022.10
set -o pipefail
mkdir -p logs/v01_external_jiang_design
python3 scripts/external_metadata/01_jiang_design_overlap.py \
  --workbook metadata/v01_external_catalog/41556_2025_1622_MOESM3_ESM.xlsx \
  --expected-workbook-sha256 2037dcde3ce88656e334f9c498c3f2571c05ae08727ea7fea69b1c7e29f070af \
  --task-panel qc/v01_task_panel/shared_perturbation_tasks.tsv \
  --expected-tasks 1264 \
  --output-dir qc/v01_external_jiang_design \
  2>&1 | tee logs/v01_external_jiang_design/retry_tgfb1.log
(cd qc/v01_external_jiang_design && sha256sum --status -c SHA256SUMS)
```

Review `summary.json`, `parameters.json`, `design_overlap.tsv` and
`SHA256SUMS` before committing the small derived outputs. The script refuses
to overwrite an existing versioned result. Its synthetic five-pathway fixture
and hash verification passed before publication. The design-target overlap
is not evidence of observed per-guide cells, matching controls, independent
biological replicates, shared expression features, or predictive performance.
Do not select a target based on external response or DE statistics; external
validation requires freezing the model, calibration, risk rule and evaluation
protocol first.

## Diagnosed first-run failure

The first HPC run stopped before creating the output directory because the
actual ST1f pathway cell for AKT1 is `TGFB1`, while the initial script required
`TGFB`. The updated script accepts either complete set of five original
ST1f pathway labels and reports the **observed raw label**; it never changes
target gene symbols. An unexpected sixth or missing pathway still raises an
error. The prior traceback remains in `run.log`; the retry uses a separate
`retry_tgfb1.log` file. No result from the first run was published.
