#!/usr/bin/env bash
# Pinned Figshare v3 metadata intake; run after loading Python 3.9.
set -e -o pipefail
cd "$HOME/VC_reliability_repo"

vc_file=metadata/v01_external_catalog/HCT116_filtered_guide_calls_per_cell.csv.gz
vc_part="$vc_file.partial"
vc_run=logs/v02_xatlas_guide_calls
vc_url=https://api.figshare.com/v2/file/download/59490731
vc_size=69606079
vc_md5=2f6a6c751e14c600dcc833964e90b1b2
vc_script=scripts/external_metadata/04_xatlas_hct116_intake.sh
vc_source=manifests/v01_xatlas_guide_count_source.tsv

mkdir -p logs metadata/v01_external_catalog
if [ -e "$vc_run" ]; then
    echo "RUN_DIRECTORY_EXISTS: inspect $vc_run before retry" >&2
    exit 1
fi
mkdir "$vc_run"

{
    date -u '+START_UTC=%Y-%m-%dT%H:%M:%SZ'
    git rev-parse HEAD
    python3 - <<'PY'
import sys
print("PYTHON:", sys.version.split()[0])
if sys.version_info[:2] != (3, 9):
    raise SystemExit("Python 3.9 required: load anaconda3/2022.10 first")
PY
    curl --version

    printf '%s\t%s\n' \
      parameter value \
      record_id 29190726 \
      figshare_version 3 \
      file_id 59490731 \
      file_name HCT116_filtered_guide_calls_per_cell.csv.gz \
      expected_size_bytes "$vc_size" \
      expected_md5 "$vc_md5" \
      requested_download_url "$vc_url" \
      source_manifest "$vc_source" \
      preview_rows 3 \
      max_time_seconds_per_attempt 900 \
      retries 2 \
      connect_timeout_seconds 20 \
      curl_max_filesize_bytes 100000000 \
      output_file "$vc_file" \
      > "$vc_run/parameters.tsv"

    if [ -e "$vc_file" ]; then
        vc_candidate="$vc_file"
        echo "EXISTING_FILE: verifying pinned size, MD5 and gzip integrity"
    else
        if [ -e "$vc_part" ]; then
            echo "PARTIAL_EXISTS: inspect $vc_part before retry" >&2
            exit 1
        fi
        vc_candidate="$vc_part"
        echo "DOWNLOADING: HCT116 guide-count metadata (69,606,079 expected bytes)"
        echo "The transfer may be quiet until it finishes or reports an error."
        curl --silent --show-error --fail --location \
          --connect-timeout 20 --max-time 900 --retry 2 \
          --max-filesize 100000000 \
          --write-out 'HTTP=%{http_code} BYTES=%{size_download} REDIRECTS=%{num_redirects}\n' \
          "$vc_url" --output "$vc_part"
    fi

    vc_bytes=$(wc -c < "$vc_candidate")
    printf 'SIZE_BYTES=%s EXPECTED_BYTES=%s\n' "$vc_bytes" "$vc_size"
    if [ "$vc_bytes" -ne "$vc_size" ]; then
        echo "SIZE_MISMATCH: preserve the candidate file and diagnose" >&2
        exit 1
    fi
    printf '%s  %s\n' "$vc_md5" "$vc_candidate" | md5sum -c -
    gzip -t -- "$vc_candidate"
    echo "SIZE_MD5_GZIP_OK"
    sha256sum -- "$vc_candidate"

    python3 - "$vc_candidate" <<'PY'
import csv
import gzip
import sys
from itertools import islice

with gzip.open(sys.argv[1], mode="rt", encoding="utf-8-sig", newline="") as handle:
    reader = csv.reader(handle)
    header = next(reader)
    print("HEADER_FIELD_COUNT:", len(header))
    print("HEADER_FIRST_20:", repr(header[:20]))
    for number, row in enumerate(islice(reader, 3), 1):
        print("ROW", number, "FIELD_COUNT:", len(row),
              "FIRST_12:", repr(row[:12]))
PY

    # Publish without replacing an existing final file; both paths share a filesystem.
    if [ "$vc_candidate" != "$vc_file" ]; then
        ln -- "$vc_candidate" "$vc_file"
        rm -- "$vc_candidate"
    fi
    echo "VALIDATION_AND_PREVIEW_OK"
} 2>&1 | tee "$vc_run/run.log"

# Paths in this manifest are relative to the repository root.
sha256sum -- "$vc_file" "$vc_run/parameters.tsv" "$vc_run/run.log" \
  "$vc_script" "$vc_source" > "$vc_run/SHA256SUMS"
sha256sum --status -c "$vc_run/SHA256SUMS"
echo "HCT116_METADATA_INTAKE_OK"
