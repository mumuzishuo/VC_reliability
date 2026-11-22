#!/usr/bin/env bash
# Download and verify the Figshare v3 HEK293T guide-call metadata.
set -euo pipefail
cd "$HOME/VC_reliability_repo"
module load anaconda3/2022.10

file="metadata/v01_external_catalog/HEK293T_filtered_guide_calls_per_cell.csv.gz"
partial="$file.partial"
log_dir=logs/v01_xatlas_hek293t_observed
url="https://api.figshare.com/v2/file/download/59490734"
expected_bytes=90886568
expected_md5=f2106fd8fdca7350160e9e5a91c67be8

mkdir -p "$log_dir" metadata/v01_external_catalog
log_file=$(mktemp "$log_dir/intake_XXXXXX.log")
echo "Intake log: $log_file"

{
    date -u '+START_UTC=%Y-%m-%dT%H:%M:%SZ'
    git rev-parse HEAD
    printf 'URL=%s EXPECTED_BYTES=%s EXPECTED_MD5=%s\n' "$url" "$expected_bytes" "$expected_md5"

    candidate="$file"
    if [ ! -e "$file" ]; then
        if ! (set -C; : > "$partial"); then
            echo "PARTIAL_EXISTS: inspect $partial before retrying" >&2
            exit 1
        fi
        candidate="$partial"
        echo "Downloading HEK293T guide metadata (90,886,568 bytes)"
        curl --silent --show-error --fail --location --retry 2 \
            --connect-timeout 20 --max-time 900 --max-filesize 100000000 \
            --write-out 'HTTP=%{http_code} BYTES=%{size_download} REDIRECTS=%{num_redirects}\n' \
            --output "$partial" "$url"
    fi

    bytes=$(wc -c < "$candidate")
    if [ "$bytes" -ne "$expected_bytes" ]; then
        echo "SIZE_MISMATCH: found $bytes bytes; expected $expected_bytes" >&2
        exit 1
    fi
    printf '%s  %s\n' "$expected_md5" "$candidate" | md5sum -c -
    gzip -t -- "$candidate"

    if [ "$candidate" = "$partial" ]; then
        ln -- "$partial" "$file"
        rm -- "$partial"
    fi
    echo "HEK293T_SIZE_MD5_GZIP_OK"
    sha256sum -- "$file" | tee "$log_dir/guide_calls.sha256"
} 2>&1 | tee "$log_file"
