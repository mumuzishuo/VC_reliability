# Step 4. Export guide-to-gene mappings

Update the repository.

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
module load anaconda3/2022.10
mkdir -p metadata/v01_label_mapping
```

Run the decoded mapping export.

```bash
python3 scripts/data_audit/03_export_label_mapping.py \
  --input "HepG2=$HOME/VC_reliability/raw/GSE264667_hepg2_raw_singlecell_01.h5ad" \
  --input "Jurkat=$HOME/VC_reliability/raw/GSE264667_jurkat_raw_singlecell_01.h5ad" \
  --input "K562_essential=$HOME/VC_reliability/raw/K562_essential_raw_singlecell_01.h5ad" \
  --input "RPE1=$HOME/VC_reliability/raw/rpe1_raw_singlecell_01.h5ad" \
  --output-dir metadata/v01_label_mapping
```

Check the generated files.

```bash
ls -lh metadata/v01_label_mapping

head -n 8 metadata/v01_label_mapping/guide_to_gene_all_environments.tsv
head -n 8 metadata/v01_label_mapping/gene_summary_all_environments.tsv
head -n 12 metadata/v01_label_mapping/cross_environment_gene_id_check.tsv
```

The mapping files use decoded gene symbols and Ensembl IDs. They preserve guide-level counts, gene-level counts and GEM-group counts for each environment.
