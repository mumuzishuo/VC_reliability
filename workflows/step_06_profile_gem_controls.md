# Step 6. Profile GEM-level controls and perturbations

Update the repository.

```bash
cd "$HOME/VC_reliability_repo"
git pull --ff-only origin dev
module load anaconda3/2022.10
mkdir -p qc/v01_gem_control
```

Run the GEM-level profiling.

```bash
python3 scripts/data_audit/05_profile_gem_controls.py \
  --input "HepG2=$HOME/VC_reliability/raw/GSE264667_hepg2_raw_singlecell_01.h5ad" \
  --input "Jurkat=$HOME/VC_reliability/raw/GSE264667_jurkat_raw_singlecell_01.h5ad" \
  --input "K562_essential=$HOME/VC_reliability/raw/K562_essential_raw_singlecell_01.h5ad" \
  --input "RPE1=$HOME/VC_reliability/raw/rpe1_raw_singlecell_01.h5ad" \
  --output-dir qc/v01_gem_control
```

Inspect the main output files.

```bash
ls -lh qc/v01_gem_control

head -n 12 qc/v01_gem_control/gem_control_summary_all_environments.tsv

head -n 12 qc/v01_gem_control/perturbation_cell_counts_by_gem_all_environments.tsv
```

The outputs provide GEM-specific control counts for later matched-control response calculations. This stage does not yet calculate expression responses.
