# Step 2. AnnData HDF5 structure audit

Run the h5py-only audit from the code repository. The four raw files remain in the separate raw-data directory.

```bash
cd "$HOME/VC_reliability_repo"
module load anaconda3/2022.10

mkdir -p qc/v01_raw_audit logs/v01_raw_audit

python3 scripts/data_audit/01_inspect_h5ad_h5py.py \
  --input "HepG2=$HOME/VC_reliability/raw/GSE264667_hepg2_raw_singlecell_01.h5ad" \
  --input "Jurkat=$HOME/VC_reliability/raw/GSE264667_jurkat_raw_singlecell_01.h5ad" \
  --input "K562_essential=$HOME/VC_reliability/raw/K562_essential_raw_singlecell_01.h5ad" \
  --input "RPE1=$HOME/VC_reliability/raw/rpe1_raw_singlecell_01.h5ad" \
  --output qc/v01_raw_audit/h5ad_structure.json \
  --top-n 10

python3 -m json.tool \
  qc/v01_raw_audit/h5ad_structure.json | sed -n '1,220p'

ls -lh qc/v01_raw_audit/h5ad_structure.json
```

This audit reads HDF5 metadata and selected observation columns. It does not load the full count matrix or modify raw files.

Record the JSON output, the file size and any error message.
