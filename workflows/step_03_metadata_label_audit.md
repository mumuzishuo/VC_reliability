# Step 3. Perturbation label and metadata audit

Update the code repository before running this step.

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
mkdir -p qc/v01_metadata_audit
```

Run the metadata profiling script.

```bash
module load anaconda3/2022.10

python3 scripts/data_audit/02_profile_obs_fields.py \
  --input "HepG2=$HOME/VC_reliability/raw/GSE264667_hepg2_raw_singlecell_01.h5ad" \
  --input "Jurkat=$HOME/VC_reliability/raw/GSE264667_jurkat_raw_singlecell_01.h5ad" \
  --input "K562_essential=$HOME/VC_reliability/raw/K562_essential_raw_singlecell_01.h5ad" \
  --input "RPE1=$HOME/VC_reliability/raw/rpe1_raw_singlecell_01.h5ad" \
  --output qc/v01_metadata_audit/obs_field_profile.json \
  --top-n 20

python3 -m json.tool \
  qc/v01_metadata_audit/obs_field_profile.json | sed -n '1,280p'
```

The output checks whether `sgID_AB` maps consistently to `gene` and `gene_id`. It also reports label counts, GEM-group counts, UMI summaries and mitochondrial-fraction summaries.

Do not construct the perturbation response matrix until the label mapping has been checked.
