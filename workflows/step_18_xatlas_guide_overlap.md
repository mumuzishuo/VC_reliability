# Step 18: X-Atlas guide-design overlap with the 1,264 tasks

Input: the verified local Figshare+ v3 `guide_library.csv`, file ID 57368587,
3,648,784 bytes, MD5 `bf1181127b314affe95e5b5ca97bdbf0`, SHA256
`1841320684c72f4a6c5f51a98ebaab0054711367a9ee81cc2167dd3fbbd6381d`.
The first exploratory reader mistakenly used `target_gene`; the actual CSV
header is `target_gene_name`. Therefore its printed zero unique symbols and
20,890 missing symbols are **not** data quality findings. Actual labels and
missingness must be recomputed with the correct column.

The new script matches `target_gene_id` to the current panel's Ensembl IDs
**exactly** and reports symbol matches and disagreements separately. It does
not normalize aliases, strip Ensembl versions, infer controls from gene IDs,
read single-cell expression or use external response outcomes. Guide pair
counts are designed guide pair counts, not observed cell support.

Run on the login node with Python 3.9; the input is 3.65 MB (no Slurm needed):

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
module load anaconda3/2022.10
set -o pipefail
sha256sum --status -c qc/v01_task_panel/SHA256SUMS
mkdir -p logs/v01_xatlas_guide_design
python3 scripts/external_metadata/02_xatlas_guide_overlap.py --guide-library metadata/v01_external_catalog/guide_library.csv --expected-sha256 1841320684c72f4a6c5f51a98ebaab0054711367a9ee81cc2167dd3fbbd6381d --expected-md5 bf1181127b314affe95e5b5ca97bdbf0 --expected-size 3648784 --expected-rows 20890 --task-panel qc/v01_task_panel/shared_perturbation_tasks.tsv --expected-tasks 1264 --output-dir qc/v01_xatlas_guide_design 2>&1 | tee logs/v01_xatlas_guide_design/run.log
(cd qc/v01_xatlas_guide_design && sha256sum -c SHA256SUMS)
```

Review `summary.json` plus `overlap_by_gene_id.tsv` and
`symbol_only_diagnostic.tsv` before publishing. The script also writes
`parameters.json` and `SHA256SUMS`; it refuses to overwrite an existing
version. Its synthetic fixture checks exact ID matching, symbol mismatches,
non-Ensembl IDs, repeated guide-pair identifiers and output checksums.
A target design match alone does not establish matching controls, guide
capture, observed cell support, biological replication or expression-gene
coverage. No large H5ADs should be downloaded at this stage.
