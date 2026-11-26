#!/usr/bin/env python3
"""Probe only the cell index and sample metadata in the Figshare v3 H5ADs."""

import argparse
import csv
import gzip
import io
import json
import subprocess
import sys
from tempfile import TemporaryFile
from itertools import islice
from pathlib import Path
from urllib.parse import urlsplit

import h5py

CATALOG = Path("metadata/v01_external_catalog")
CELL_LINES = ("HCT116", "HEK293T")
PROBE_ROWS = 16
GUIDE_WINDOW = 4096
MAX_BYTES = 8 * 1024 * 1024
MAX_REQUESTS = 256
CONNECT_TIMEOUT = 8
REQUEST_TIMEOUT = 25
CURL_WRITEOUT = (
    "%{stderr}%{http_code}\n%header{content-range}\n"
    "%header{content-encoding}\n%header{etag}\n%{url_effective}\n"
)


def curl_range(url, start, count):
    """Return one bounded range and its final response metadata."""
    end = start + count - 1
    byte_range = f"{start}-{end}"
    command = [
        "curl", "--disable", "--silent", "--fail", "--location", "--globoff",
        "--proto", "=https", "--proto-redir", "=https", "--max-redirs", "4",
        "--connect-timeout", str(CONNECT_TIMEOUT), "--max-time", str(REQUEST_TIMEOUT),
        "--range", byte_range, "--max-filesize", str(max(count, 2048)),
        "--header", "Accept-Encoding: identity",
        "--output", "-", "--write-out", CURL_WRITEOUT, "--", url,
    ]
    with TemporaryFile() as metadata:
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=metadata) as process:
            data = process.stdout.read(count + 1)
            if len(data) > count:
                process.kill()
                raise OSError(f"Response exceeded the requested range {byte_range}")
            exit_code = process.wait()
        metadata.seek(0)
        fields = metadata.read().decode("utf-8", errors="replace").splitlines()

    status = fields[0] if fields else "000"
    if exit_code != 0 or status != "206":
        raise OSError(f"Range {byte_range}: curl exit {exit_code}, HTTP {status}, received {len(data)} bytes")
    if len(fields) != 5:
        raise OSError(f"Range {byte_range}: missing curl response metadata")
    _, content_range, encoding, etag, effective_url = fields
    return data, content_range, encoding.lower(), etag or None, effective_url


class RangeReader(io.RawIOBase):
    """Expose exact HTTPS byte ranges to h5py as a read-only file."""

    def __init__(self, url, size):
        super().__init__()
        self.url = url
        self.size = size
        self.position = 0
        self.requests = 0
        self.bytes_read = 0
        self.etag = None

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=io.SEEK_SET):
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self.position, io.SEEK_END: self.size}
        position = base[whence] + offset
        if position < 0:
            raise ValueError("Negative file position")
        self.position = position
        return position

    def readinto(self, buffer):
        start = self.position
        count = min(len(buffer), self.size - start)
        if count <= 0:
            return 0
        if self.bytes_read + count > MAX_BYTES or self.requests >= MAX_REQUESTS:
            raise OSError("Metadata probe exceeded its byte or request limit")

        end = start + count - 1
        expected_range = f"bytes {start}-{end}/{self.size}"
        self.requests += 1
        data, actual_range, encoding, etag, effective_url = curl_range(self.url, start, count)
        if actual_range != expected_range or encoding not in ("", "identity"):
            raise OSError(f"Incorrect range or encoding for bytes {start}-{end}")
        if urlsplit(effective_url).scheme != "https":
            raise OSError("Expected an HTTPS download")
        if self.etag is not None and etag != self.etag:
            raise OSError("The remote file changed during the probe")
        self.etag = etag
        # Start every range at the official API URL. Reused redirects may reject
        # later requests, as observed in the v04 HCT116 probe.
        if len(data) != count:
            raise OSError(f"Short HTTP range: expected {count} bytes, read {len(data)}")
        buffer[:count] = data
        self.position += count
        self.bytes_read += count
        return count


def sample_labels(obs, count):
    sample = obs["sample"]
    if isinstance(sample, h5py.Group):
        codes = sample["codes"][:count]
        categories = sample["categories"].asstr()[:]
    elif "categories" in sample.attrs:
        codes = sample[:count]
        categories = obs.file[sample.attrs["categories"]].asstr()[:]
    else:
        return sample.asstr()[:count].tolist()
    return [categories[code] if code >= 0 else None for code in codes]


def probe_cell_line(cell_line, files):
    expression = files[f"{cell_line}_filtered_dual_guide_cells.h5ad"]
    guide = files[f"{cell_line}_filtered_guide_calls_per_cell.csv.gz"]
    qc = Path(f"qc/v01_xatlas_{cell_line.lower()}_observed/summary.json")
    guide_summary = json.loads(qc.read_text(encoding="utf-8"))
    if guide_summary["cell_line"] != cell_line or guide_summary["status"] != "GUIDE_METADATA_ONLY":
        raise ValueError(f"Unexpected guide audit summary in {qc}")
    guide_rows = guide_summary["unique_cells_in_guide_csv"]
    with gzip.open(CATALOG / guide["name"], "rt", encoding="utf-8-sig", newline="") as stream:
        guide_barcodes = [
            row["cell_barcode"] for row in islice(csv.DictReader(stream), GUIDE_WINDOW)
        ]

    url = f"https://api.figshare.com/v2/file/download/{expression['id']}"
    print(f"Checking {cell_line} .obs metadata from Figshare file {expression['id']}", flush=True)
    with RangeReader(url, expression["size"]) as source:
        source.read(8)
        print(f"COMPUTE_RANGE_OK {cell_line} bytes=8", flush=True)
        source.seek(0)
        with h5py.File(source, "r") as handle:
            obs = handle["obs"]
            index_name = obs.attrs["_index"]
            if isinstance(index_name, bytes):
                index_name = index_name.decode("utf-8")
            index = obs[index_name]
            obs_rows = len(index)
            barcodes = index.asstr()[:PROBE_ROWS].tolist()
            if not barcodes or not guide_barcodes:
                raise ValueError(f"Empty barcode preview for {cell_line}")
            samples = sample_labels(obs, len(barcodes))
            columns = sorted(name for name in obs.keys() if name not in (index_name, "__categories"))

    guide_prefix = set(guide_barcodes)
    return {
        "cell_line": cell_line,
        "h5ad_file_id": expression["id"],
        "h5ad_size_bytes": expression["size"],
        "h5ad_source_url": url,
        "guide_csv_file_id": guide["id"],
        "guide_csv_unique_cells": guide_rows,
        "h5ad_obs_rows": obs_rows,
        "same_cohort_size": obs_rows == guide_rows,
        "obs_index": index_name,
        "obs_columns": columns,
        "obs_barcodes_checked": len(barcodes),
        "guide_prefix_rows_checked": len(guide_barcodes),
        "obs_prefix_barcodes_found_in_guide_prefix": sum(b in guide_prefix for b in barcodes),
        "same_barcode_order_in_preview": barcodes == guide_barcodes[:len(barcodes)],
        "sample_matches_barcode_suffix_in_preview": sum(
            barcode.partition("-")[2] == sample for barcode, sample in zip(barcodes, samples)
        ),
        "range_requests": source.requests,
        "remote_bytes_read": source.bytes_read,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    curl_version = subprocess.check_output(["curl", "--version"], text=True).splitlines()[0]
    version = tuple(int(part) for part in curl_version.split()[1].split(".")[:3])
    if version < (7, 84, 0):
        raise RuntimeError("Load the Anaconda curl 7.84.0 environment before running the probe")
    print(curl_version, flush=True)
    record = json.loads((CATALOG / "xatlas_figshare_article.json").read_text(encoding="utf-8"))
    if record["id"] != 29190726 or record["version"] != 3:
        raise ValueError("Expected the archived Figshare+ v3 record")
    files = {entry["name"]: entry for entry in record["files"]}
    args.output_dir.mkdir(parents=True, exist_ok=False)

    summary = {
        "status": "OBS_INDEX_PROBE_ONLY",
        "figshare_version_doi": record["doi"],
        "transport": "curl_https_range",
        "curl_version": curl_version,
        "python_version": sys.version.split()[0],
        "h5py_version": h5py.__version__,
        "obs_preview_rows": PROBE_ROWS,
        "guide_prefix_rows": GUIDE_WINDOW,
        "max_remote_bytes_per_file": MAX_BYTES,
        "max_range_requests_per_file": MAX_REQUESTS,
        "connect_timeout_seconds": CONNECT_TIMEOUT,
        "request_timeout_seconds": REQUEST_TIMEOUT,
        "cohorts": [probe_cell_line(cell_line, files) for cell_line in CELL_LINES],
    }
    report = json.dumps(summary, indent=2)
    (args.output_dir / "summary.json").write_text(report + "\n", encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
