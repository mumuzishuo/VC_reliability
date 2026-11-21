# Step 22: X-Atlas expression-feature review and archive

Current project calculation from the verified author
`metadata/gene_metadata.parquet` (38,606 nonempty, unique Ensembl IDs)
and the audited four-context 6,642-feature panel:

- 6,638 / 6,642 (99.94%) have an exact, unique Ensembl ID to
  `gene_token_id` mapping in the external **feature dictionary**.
- Four panel IDs have no exact external ID:
  `ENSG00000112096`, `ENSG00000130723`,
  `ENSG00000167747`, and `ENSG00000196312`.
- There were no duplicate external IDs, reused tokens, missing tokens,
  versioned Ensembl IDs, or versionless candidates in this check.
- `ENSG00000112096` also occurs among the earlier guide-design
  **symbol-only** matches under `SOD2`, where the external guide-library
  target ID differs. Do not replace this expression ID from the symbol;
  neither the guide nor expression ID is reconciled here.

These are *feature-dictionary* findings, not cell-line-specific expression
measurements. The provisional external evaluation feature set is the exact
6,638-ID intersection marked `exact_unique_token` in
`feature_overlap.tsv`. Freeze this exact-ID rule and the derived set before
reading external target responses. Once external data are eligible, compute
all predictions, zero-response baselines, metrics, and uncertainty on the
same fixed feature set. Keep the four unmatched IDs and all label exceptions
out of the exact-match set; record any ID-resolution analyses separately.

On HPC, verify the files before archiving. The older feature SHA256SUMS
contains paths from the repository root; the new output SHA256SUMS contains
filenames relative to its own directory.

```bash
(
  set -e -o pipefail
  cd "$HOME/VC_reliability_repo"
  git pull --ff-only origin dev
  sha256sum --status -c qc/v01_expression_features/SHA256SUMS
  (cd qc/v01_xatlas_expression_feature_overlap && sha256sum --status -c SHA256SUMS)

  vc_files=(
    qc/v01_xatlas_expression_feature_overlap/feature_overlap.tsv
    qc/v01_xatlas_expression_feature_overlap/summary.json
    qc/v01_xatlas_expression_feature_overlap/parameters.json
    qc/v01_xatlas_expression_feature_overlap/SHA256SUMS
    logs/v01_xatlas_expression_features/metadata_check.log
    logs/v01_xatlas_expression_features/overlap_run.log
    logs/v01_xatlas_expression_features/unmatched_review.log
  )
  vc_manifest=manifests/v01_xatlas_expression_feature_review.sha256
  mkdir -p manifests
  if [ ! -e "$vc_manifest" ]; then
    sha256sum -- "${vc_files[@]}" > "$vc_manifest"
  fi
  sha256sum --status -c "$vc_manifest"
  git add -- "${vc_files[@]}" "$vc_manifest"
  git diff --cached --check
  git diff --cached --stat
  git commit --only -m "metadata: archive X-Atlas expression feature audit" \
    -- "${vc_files[@]}" "$vc_manifest"
  git push origin dev
)
```

Do not add the original Parquet, Python wheel, unrelated QC directories,
raw H5ADs, or target outcomes. The next gate is observed cell/guide counts
and non-targeting controls by sample/GEM, independently in HCT116 and
HEK293T; source guide-library design counts are not substitutes.
