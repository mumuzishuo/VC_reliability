# Step 17: obtain and verify X-Atlas guide design table

The archived Figshare+ v3 article JSON (ID 29190726, DOI
`10.25452/figshare.plus.29190726.v3`) identifies exactly one
`guide_library.csv`: file ID 57368587, 3,648,784 bytes, reported MD5
`bf1181127b314affe95e5b5ca97bdbf0`. The documented download URL is
`https://ndownloader.figshare.com/files/57368587`. The repository
ignores the downloaded CSV; keep it on the HPC and publish only provenance
and small derived overlap outputs after review. This table describes guide
pair design and does not establish actual assigned cells, GEM-matched
controls or independent biological replicates.

On the login node (small file; Slurm is unnecessary):

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
python3 - <<'PY'
import hashlib, json
path = 'metadata/v01_external_catalog/xatlas_figshare_article.json'
with open(path, 'rb') as handle:
    raw = handle.read()
assert hashlib.sha256(raw).hexdigest() == 'effabb3539f95f22481352e88cb39aa6468b9ae92a50bc8edf4d23ec97bf79f7'
record = json.loads(raw)
assert record['id'] == 29190726 and record['version'] == 3
matches = [item for item in record['files'] if item['name'] == 'guide_library.csv']
assert len(matches) == 1
item = matches[0]
assert (item['id'], item['size'], item['supplied_md5']) == (57368587, 3648784, 'bf1181127b314affe95e5b5ca97bdbf0')
assert item['download_url'] == 'https://ndownloader.figshare.com/files/57368587'
print('FILE_METADATA_OK', item['name'], item['size'], item['download_url'])
PY
vc_guide=metadata/v01_external_catalog/guide_library.csv
if test -e "$vc_guide"; then printf '%s\n' 'GUIDE_FILE_EXISTS: inspect before retry'; else curl --fail --location --retry 2 --max-time 180 --max-filesize 5000000 'https://ndownloader.figshare.com/files/57368587' -o "$vc_guide"; fi
wc -c -- "$vc_guide"
printf '%s  %s\n' bf1181127b314affe95e5b5ca97bdbf0 "$vc_guide" | md5sum -c -
sha256sum -- "$vc_guide"
python3 - <<'PY'
import csv
path = 'metadata/v01_external_catalog/guide_library.csv'
with open(path, encoding='utf-8-sig', newline='') as handle:
    rows = csv.reader(handle)
    print('COLUMNS:', next(rows))
    print('DATA_ROWS:', sum(1 for _ in rows))
PY
```

If curl fails or the size/MD5 differs, preserve the output and diagnose the
file transfer before any analysis. Do not use a partially downloaded table.
Do not download the 209 GB or 350 GB H5ADs at this gate. The separate small
Hugging Face gene metadata will be evaluated after guide-target overlap.
