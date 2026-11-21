# Step 21: X-Atlas feature dictionary overlap

Read-only audit of X-Atlas/Orion `gene_metadata.parquet` against the four-context
6,642-gene feature panel. Uses exact `ensembl_id` only for the primary count.
Duplicate IDs, missing tokens, reused tokens, and versionless candidates are
reported separately. The feature dictionary alone cannot show whether a gene
is expressed in HCT116 or HEK293T.

Prerequisites: step 20 downloaded and SHA256-verified the author Parquet; the
current four-context feature table and `SHA256SUMS` are present on HPC.
Python 3.9.13's Anaconda module has no Parquet reader, so install a pinned
prebuilt PyArrow 16.1.0 wheel to a Git-ignored, project-local directory.
This version publishes a CPython 3.9 Linux x86-64 wheel; no system packages
or shared Anaconda installation are changed. Record its version in the output
parameters and log. Install and compute only on small metadata on the login
node; Slurm is unnecessary for this step.

```bash
(
  set -e -o pipefail
  cd "$HOME/VC_reliability_repo"
  git pull --ff-only origin dev
  module load anaconda3/2022.10
  (cd qc/v01_expression_features && sha256sum --status -c SHA256SUMS)

  vc_pkg=.local_packages/pyarrow16
  if [ ! -d "$vc_pkg/pyarrow" ]; then
    if [ -e "$vc_pkg" ]; then
      echo "PACKAGE_DIR_EXISTS: inspect failed install before retry" >&2
      exit 1
    fi
    mkdir -p .local_packages
    python3 -m pip install --disable-pip-version-check --only-binary=:all: \
      --no-deps --target "$vc_pkg" 'pyarrow==16.1.0'
  fi
  export PYTHONPATH="$PWD/$vc_pkg${PYTHONPATH:+:$PYTHONPATH}"
  mkdir -p logs/v01_xatlas_expression_features

  {
    git rev-parse HEAD
    python3 - <<'PY'
import sys, numpy, pyarrow
print("PYTHON", sys.version.split()[0])
print("NUMPY", numpy.__version__)
print("PYARROW", pyarrow.__version__)
assert sys.version_info[:2] == (3, 9)
assert pyarrow.__version__ == "16.1.0"
PY
    python3 scripts/external_metadata/03_xatlas_expression_overlap.py \
      --gene-metadata metadata/v01_external_catalog/gene_metadata.parquet \
      --expected-sha256 33900a0dbaeafb84601c7c8ba7c6642c440774fc3f49613dab11cf1f2bf0fdbf \
      --expected-size 881915 \
      --features qc/v01_expression_features/common_expression_features.tsv \
      --expected-features 6642 \
      --output-dir qc/v01_xatlas_expression_feature_overlap
    (cd qc/v01_xatlas_expression_feature_overlap && sha256sum -c SHA256SUMS)
  } 2>&1 | tee logs/v01_xatlas_expression_features/overlap_run.log
)
```

The script refuses to overwrite an existing output directory. If pip or
Parquet import fails, preserve and diagnose the output before changing the
command. Outputs are `feature_overlap.tsv`, `summary.json`,
`parameters.json`, and `SHA256SUMS`; keep the local metadata file and
reader package outside Git. After inspecting the results, archive derived QC
and logs on `dev` with a SHA256 manifest. Only then consider downloading
per-cell guide calls to establish usable perturbation cell counts and
sample/GEM control structure. Do not read external response results yet.
