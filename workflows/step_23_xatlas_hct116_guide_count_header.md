# Step 23: inspect X-Atlas HCT116 per-cell guide-count metadata

Goal: determine the schema and join keys of the **guide-count** file before
estimating observed per-perturbation support or non-targeting controls.
The Figshare+ v3 record 29190726 lists:
`HCT116_filtered_guide_calls_per_cell.csv.gz`, file ID 59490731,
69,606,079 bytes, MD5 `2f6a6c751e14c600dcc833964e90b1b2`.
The companion HEK293T file is listed in
`manifests/v01_xatlas_guide_count_source.tsv` but is not downloaded yet.
The author describes these as sgRNA count data. The structure, cell barcode
key, sample/GEM join, control encoding and relationship to filtered expression
cells have **not** been verified.

Download and preview only the HCT116 metadata file; do not download the
209 GB HCT116 expression H5AD. This is a metadata read on the login node
(Python 3.9); no Slurm allocation is needed for the header and three rows.
The downloaded file and any partial transfer are Git-ignored. Stop if size,
official MD5, or gzip/CSV parsing differs from expectations.

```bash
(
  set -e -o pipefail
  cd "$HOME/VC_reliability_repo"
  git pull --ff-only origin dev
  module load anaconda3/2022.10

  vc_file=metadata/v01_external_catalog/HCT116_filtered_guide_calls_per_cell.csv.gz
  vc_part="$vc_file.partial"
  vc_md5=2f6a6c751e14c600dcc833964e90b1b2
  mkdir -p logs/v01_xatlas_guide_calls

  if [ ! -e "$vc_file" ]; then
    if [ -e "$vc_part" ]; then
      echo "PARTIAL_EXISTS: inspect the partial transfer before retry" >&2
      exit 1
    fi
    curl --fail --location --retry 2 --max-time 900 --max-filesize 100000000 \
      "https://ndownloader.figshare.com/files/59490731" -o "$vc_part"
    test "$(wc -c < "$vc_part")" -eq 69606079
    printf '%s  %s\n' "$vc_md5" "$vc_part" | md5sum --status -c -
    mv -- "$vc_part" "$vc_file"
  fi

  test "$(wc -c < "$vc_file")" -eq 69606079
  printf '%s  %s\n' "$vc_md5" "$vc_file" | md5sum --status -c -
  {
    git rev-parse --short HEAD
    python3 --version
    wc -c -- "$vc_file"
    md5sum -- "$vc_file"
    sha256sum -- "$vc_file"
    python3 - <<'PY'
import csv
import gzip
from itertools import islice

path = "metadata/v01_external_catalog/HCT116_filtered_guide_calls_per_cell.csv.gz"
with gzip.open(path, mode="rt", encoding="utf-8-sig", newline="") as handle:
    reader = csv.reader(handle)
    header = next(reader)
    print("HEADER_FIELD_COUNT:", len(header))
    print("HEADER_FIRST_20:", repr(header[:20]))
    for number, row in enumerate(islice(reader, 3), 1):
        print("ROW", number, "FIELD_COUNT:", len(row),
              "FIRST_12:", repr(row[:12]))
PY
  } 2>&1 | tee logs/v01_xatlas_guide_calls/header_HCT116.log
)
```

This step reads only the header and three rows, not the entire decompressed
file. After reviewing the schema, write a versioned streaming QC script
before counting cells by perturbation or controls by sample/GEM. A guide
count is not a unique cell, and a non-targeting *guide* is not proof of
sufficient matched control *cells*. Keep target responses unseen.
