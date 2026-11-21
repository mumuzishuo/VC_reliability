# Step 16: X-Atlas Figshare+ metadata gate

The authors' [Figshare+ article](https://plus.figshare.com/articles/dataset/Processed_data_for_X-Atlas_Orion_Genome-wide_Perturb-seq_Datasets_via_a_Scalable_Fix-Cryopreserve_Platform_for_Training_Dose-Dependent_Biological_Foundation_Models/29190726)
is version 3 and describes filtered dual-guide HCT116 and HEK293T H5ADs,
`guide_library.csv` (target symbols and Ensembl IDs), and per-cell guide-call
CSVs. The page lists CC BY-NC-SA 4.0. The author's
[Hugging Face card](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/blob/main/README.md)
also lists a separate small `metadata/gene_metadata.parquet` mapping gene
tokens to Ensembl IDs. The full article download is hundreds of GB; start
with the public metadata API only. Figshare's
[API documentation](https://docs.figshare.com/v2/) defines
`GET /v2/articles/{article_id}` as the public article details endpoint.

On the HPC login node, fetch the JSON listing only (no Slurm needed):

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
mkdir -p metadata/v01_external_catalog
curl --fail --location --retry 2 --max-time 60 --max-filesize 5000000 "https://api.figshare.com/v2/articles/29190726" -o metadata/v01_external_catalog/xatlas_figshare_article.json
python3 -m json.tool metadata/v01_external_catalog/xatlas_figshare_article.json >/dev/null
sha256sum -- metadata/v01_external_catalog/xatlas_figshare_article.json
python3 - <<'PY'
import json
p = 'metadata/v01_external_catalog/xatlas_figshare_article.json'
with open(p, encoding='utf-8') as handle:
    record = json.load(handle)
print('RECORD_ID:', record.get('id'))
print('DOI:', record.get('doi'))
print('VERSION:', record.get('version'))
print('TITLE:', record.get('title'))
print('LICENSE:', record.get('license'))
files = record.get('files', [])
print('FILE_COUNT:', len(files))
for item in files:
    print('FILE:', item.get('name'), 'SIZE_BYTES:', item.get('size'),
          'CHECKSUM:', item.get('supplied_md5', item.get('computed_md5')))
PY
```

If the API is unavailable, preserve the HTTP status and diagnose before
trying other official endpoints. Do not download H5ADs or guide-call tables
until metadata establishes the relevant small-file sizes and valid controls.
The guide library lists designed targets, not observed per-batch cell support.
A dataset card cell count is not a biological replicate count.
