#!/bin/bash -l
#SBATCH --job-name=vc_hct116_guides
#SBATCH --account=zhenwei01_bio214
#SBATCH --partition=cpu8358
#SBATCH --qos=26cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/v01_xatlas_hct116_observed/step24_%j.out
#SBATCH --error=logs/v01_xatlas_hct116_observed/step24_%j.err

set -euo pipefail
cd "$HOME/VC_reliability_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd

out=qc/v01_xatlas_hct116_observed
inputs=manifests/v01_xatlas_hct116_observed_inputs.sha256
script=scripts/external_metadata/05_xatlas_observed_support.py
sha256sum -c "$inputs"

python3 -u "$script" \
    --cell-line HCT116 \
    --guide-calls metadata/v01_external_catalog/HCT116_filtered_guide_calls_per_cell.csv.gz \
    --guide-library metadata/v01_external_catalog/guide_library.csv \
    --task-panel qc/v01_task_panel/shared_perturbation_tasks.tsv \
    --output-dir "$out"

git rev-parse HEAD > "$out/code_commit.txt"
sha256sum -- "$out"/*.tsv "$out"/*.json "$out/code_commit.txt" "$inputs" \
    "$script" scripts/external_metadata/02_xatlas_guide_overlap.py \
    slurm/step_24_xatlas_hct116_observed.slurm > "$out/SHA256SUMS"
sha256sum --status -c "$out/SHA256SUMS"
echo "HCT116_GUIDE_AUDIT_OK $SLURM_JOB_ID"
