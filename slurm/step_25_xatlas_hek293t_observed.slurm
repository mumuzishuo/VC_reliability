#!/bin/bash -l
#SBATCH --job-name=vc_hek293t_guides
#SBATCH --account=zhenwei01_bio214
#SBATCH --partition=cpu8358
#SBATCH --qos=26cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --output=logs/v01_xatlas_hek293t_observed/step25_%j.out
#SBATCH --error=logs/v01_xatlas_hek293t_observed/step25_%j.err

set -euo pipefail
cd "$HOME/VC_reliability_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd

hct116=qc/v01_xatlas_hct116_observed
out=qc/v01_xatlas_hek293t_observed
intake_hash=logs/v01_xatlas_hek293t_observed/guide_calls.sha256
shared_inputs=manifests/v01_xatlas_hct116_observed_inputs.sha256
audit=scripts/external_metadata/05_xatlas_observed_support.py
compare=scripts/external_metadata/07_compare_xatlas_observed_support.py

sha256sum -c "$intake_hash" "$shared_inputs" "$hct116/SHA256SUMS"

python3 -u "$audit" \
    --cell-line HEK293T \
    --guide-calls "metadata/v01_external_catalog/HEK293T_filtered_guide_calls_per_cell.csv.gz" \
    --guide-library metadata/v01_external_catalog/guide_library.csv \
    --task-panel qc/v01_task_panel/shared_perturbation_tasks.tsv \
    --output-dir "$out"

python3 -u "$compare" "$hct116" "$out"

git rev-parse HEAD > "$out/code_commit.txt"
sha256sum -- "$out"/*.tsv "$out"/*.json "$out/code_commit.txt" \
    "$intake_hash" "$shared_inputs" "$hct116/SHA256SUMS" \
    "$audit" "$compare" scripts/external_metadata/02_xatlas_guide_overlap.py \
    scripts/external_metadata/06_xatlas_hek293t_intake.sh \
    slurm/step_25_xatlas_hek293t_observed.slurm > "$out/SHA256SUMS"
sha256sum --status -c "$out/SHA256SUMS"
echo "XATLAS_PAIRED_GUIDE_AUDIT_OK $SLURM_JOB_ID"
