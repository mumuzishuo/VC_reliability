# Step 20: X-Atlas expression-gene metadata intake

Primary author source: [X-Atlas/Orion gene_metadata.parquet](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/blob/main/metadata/gene_metadata.parquet).
The author reports an 881,915-byte Parquet file with SHA256
`33900a0dbaeafb84601c7c8ba7c6642c440774fc3f49613dab11cf1f2bf0fdbf`.
Its dataset card describes `ensembl_id`, `gene_name`, and `gene_token_id`.
This is a *feature dictionary*. Its ID overlap with the existing 6,642-gene
panel does not establish expression or effective counts in either cell line.

Download and inspect this small file on the HPC login node. The repository
ignores the source Parquet and partial transfer; retain them locally. Do not
download any large expression matrix yet. Python 3.9 is used.

```bash
(
  set -e -o pipefail
  cd "$HOME/VC_reliability_repo"
  git pull --ff-only origin dev
  module load anaconda3/2022.10
  vc_gene=metadata/v01_external_catalog/gene_metadata.parquet
  vc_sha=33900a0dbaeafb84601c7c8ba7c6642c440774fc3f49613dab11cf1f2bf0fdbf
  vc_url=https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/resolve/main/metadata/gene_metadata.parquet
  mkdir -p logs/v01_xatlas_expression_features
  if [ ! -e "$vc_gene" ]; then
    if [ -e "$vc_gene.partial" ]; then
      echo "PARTIAL_EXISTS: inspect transfer before retry" >&2
      exit 1
    fi
    curl --fail --location --retry 2 --max-time 120 --max-filesize 2000000 \
      "$vc_url" -o "$vc_gene.partial"
    test "$(wc -c < "$vc_gene.partial")" -eq 881915
    printf '%s  %s\n' "$vc_sha" "$vc_gene.partial" | sha256sum --status -c -
    mv -- "$vc_gene.partial" "$vc_gene"
  fi
  test "$(wc -c < "$vc_gene")" -eq 881915
  printf '%s  %s\n' "$vc_sha" "$vc_gene" | sha256sum --status -c -
  {
    git rev-parse --short HEAD
    python3 --version
    wc -c -- "$vc_gene"
    sha256sum -- "$vc_gene"
    python3 - <<'PY'
import importlib.util
for module in ("pyarrow", "fastparquet", "duckdb"):
    print("READER_AVAILABLE", module, importlib.util.find_spec(module) is not None)
if importlib.util.find_spec("pyarrow") is not None:
    import pyarrow.parquet as pq
    pf = pq.ParquetFile("metadata/v01_external_catalog/gene_metadata.parquet")
    print("PARQUET_ROWS", pf.metadata.num_rows)
    print("PARQUET_COLUMNS", pf.schema.names)
PY
  } 2>&1 | tee logs/v01_xatlas_expression_features/metadata_check.log
)
```

Stop if the size or SHA256 differs: the remote file may have changed or a
download may have returned a pointer/error document. Next use the audited
feature dictionary and the existing
`qc/v01_expression_features/common_expression_features.tsv` to calculate
**exact Ensembl ID overlap**, detect duplicate/missing IDs and report any
versioned IDs separately. Compute cell-line-specific detectability only after
cell-level metadata is available. Keep external responses unseen and all
predictor/calibration/risk rules frozen before external evaluation.
