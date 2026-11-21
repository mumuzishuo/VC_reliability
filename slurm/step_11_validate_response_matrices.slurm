#!/bin/bash -l
#SBATCH --job-name=vc_matrix_qc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/v01_response_matrix_validation/step11_%j.out
#SBATCH --error=logs/v01_response_matrix_validation/step11_%j.err

set -euo pipefail
vc_repo="$HOME/VC_reliability_repo"
vc_old="$HOME/VC_reliability"
cd "$vc_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

vc_out=qc/v01_response_matrix_validation
mkdir -p "$vc_out" logs/v01_response_matrix_validation
if [ -e "$vc_out/response_matrix_validation.json" ] || \
   [ -e "$vc_out/response_matrix_validation_by_task.tsv" ]; then
    echo "Validation outputs exist: inspect them before rerunning; refusing overwrite" >&2
    exit 1
fi

git rev-parse HEAD > "$vc_out/code_commit.txt"
python3 --version > "$vc_out/python_version.txt" 2>&1
scontrol show job "$SLURM_JOB_ID" > "$vc_out/slurm_job.txt"
cat > "$vc_out/run_parameters.json" <<'JSON'
{
  "script": "scripts/data_audit/10_validate_response_matrices.py",
  "arithmetic_tolerance": 0.000005,
  "low_control_fraction_tolerance": 0.000001,
  "low_control_diagnostic_cutoff": 30,
  "chunk_size": 64,
  "cpus_per_task": 2,
  "memory": "8G",
  "time_limit": "01:00:00"
}
JSON

(cd qc/v01_gem_control && sha256sum --status -c SHA256SUMS)
sha256sum --status -c qc/v01_task_panel/SHA256SUMS
sha256sum --status -c qc/v01_expression_features/SHA256SUMS
sha256sum --status -c qc/v01_existing_responses/SHA256SUMS
sha256sum --status -c qc/v01_response_provenance/SHA256SUMS
(cd "$vc_old" && sha256sum --status -c "$vc_repo/qc/v01_existing_responses/legacy_inputs.sha256")
(cd "$vc_old" && sha256sum --status -c "$vc_repo/qc/v01_response_provenance/method_inputs.sha256")

python3 -u scripts/data_audit/10_validate_response_matrices.py \
    --features qc/v01_expression_features/common_expression_features.tsv \
    --tasks qc/v01_task_panel/shared_perturbation_tasks.tsv \
    --manifest "$vc_old/processed/v01/manifest.json" \
    --provenance qc/v01_response_provenance/response_provenance.json \
    --input "HepG2=$vc_old/processed/v01/evaluation_only/HepG2_responses.h5" \
    --input "Jurkat=$vc_old/processed/v01/sources/Jurkat_responses.h5" \
    --input "K562_essential=$vc_old/processed/v01/sources/K562_responses.h5" \
    --input "RPE1=$vc_old/processed/v01/sources/RPE1_responses.h5" \
    --output-dir "$vc_out"

sha256sum -- \
    "$vc_out/response_matrix_validation.json" \
    "$vc_out/response_matrix_validation_by_task.tsv" \
    "$vc_out/run_parameters.json" \
    "$vc_out/code_commit.txt" \
    "$vc_out/python_version.txt" \
    "$vc_out/slurm_job.txt" \
    scripts/data_audit/10_validate_response_matrices.py \
    slurm/step_11_validate_response_matrices.slurm \
    > "$vc_out/SHA256SUMS"
sha256sum --status -c "$vc_out/SHA256SUMS"
echo "MATRIX_QC_COMPLETE $SLURM_JOB_ID"
