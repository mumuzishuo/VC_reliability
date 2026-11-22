#!/usr/bin/env bash
# Check the first 8 bytes of the pinned Figshare v3 H5ADs from the login node.
set -euo pipefail
cd "$(dirname -- "$0")/../.."
module load anaconda3/2022.10

source_manifest=manifests/v01_xatlas_figshare_article.sha256
sha256sum --status -c "$source_manifest"

log_dir=logs/v02_xatlas_expression_index_probe
mkdir -p "$log_dir"
log_file=$(mktemp "$log_dir/range_XXXXXX.log")

probe() {
    local cell_line=$1 file_id=$2 file_size=$3 report curl_status=0
    echo "FILE $cell_line ID=$file_id EXPECTED_RANGE=bytes 0-7/$file_size"

    # Keep only status and range headers. Redirect URLs may contain temporary signatures.
    report=$(curl --disable --silent --show-error --location \
        --proto '=https' --proto-redir '=https' --max-redirs 4 \
        --range 0-7 --max-filesize 2048 --limit-rate 1K \
        --connect-timeout 8 --max-time 25 --retry 0 \
        --dump-header - --output /dev/null \
        --write-out 'HTTP=%{http_code} BYTES=%{size_download} REDIRECTS=%{num_redirects}\n' \
        "https://api.figshare.com/v2/file/download/$file_id" |
        awk '
            { sub(/\r$/, ""); header = tolower($0) }
            /^HTTP\// { range = ""; print }
            header ~ /^content-range:/ {
                sub(/^content-range:[[:space:]]*/, "content-range: ", header)
                sub(/[[:space:]]+$/, "", header)
                range = header
            }
            /^HTTP=/ { if (range != "") print range; print }
        ') || curl_status=$?

    printf '%s\nCURL_EXIT=%s\n' "$report" "$curl_status"
    if [[ "$curl_status" == 0 &&
          "$report" == *"content-range: bytes 0-7/$file_size"$'\n'* &&
          "$report" == *"HTTP=206 BYTES=8 REDIRECTS="* ]]; then
        echo "RANGE_OK $cell_line"
    else
        echo "RANGE_UNCONFIRMED $cell_line"
    fi
}

{
    date -u '+START_UTC=%Y-%m-%dT%H:%M:%SZ'
    git rev-parse HEAD
    curl --version | sed -n '1p'
    probe HCT116 55021257 209354246272
    probe HEK293T 55074802 350164035901
} 2>&1 | tee "$log_file"

sha256sum -- "$log_file" scripts/external_metadata/09_probe_xatlas_range.sh \
    "$source_manifest" > "$log_file.sha256"
echo "Saved $log_file and $log_file.sha256"
