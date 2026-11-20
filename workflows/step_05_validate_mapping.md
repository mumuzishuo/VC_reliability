# Step 5. Validate cross-environment mappings

Update the repository.

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
module load anaconda3/2022.10
mkdir -p qc/v01_mapping_validation
```

Run the mapping validation.

```bash
python3 scripts/data_audit/04_validate_mapping.py \
  --mapping-dir metadata/v01_label_mapping \
  --output qc/v01_mapping_validation/mapping_summary.json
```

Print the summary.

```bash
python3 -m json.tool \
  qc/v01_mapping_validation/mapping_summary.json
```

The main checks are:

- non-control gene counts per environment
- guide counts per environment
- control cell counts
- coverage distribution across one, two, three and four environments
- gene ID consistency
- number of shared non-control genes across all four environments

The shared-gene count must be compared with the current project value of 2,053. Do not change the task panel until the comparison is explained.
