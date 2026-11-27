#!/bin/bash -l
#SBATCH --job-name=vc_four_source
#SBATCH --account=zhenwei01_bio214
#SBATCH --partition=cpu8358
#SBATCH --qos=26cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/v01_four_source_calibration/step29_%j.out
#SBATCH --error=logs/v01_four_source_calibration/step29_%j.err

set -euo pipefail
cd "$HOME/VC_reliability_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2

vc_script=scripts/development/01_four_source_calibration.py
vc_job=slurm/step_29_four_source_calibration.slurm
vc_out="results/v01_four_source_calibration/job_${SLURM_JOB_ID}"
echo "OUTPUT_DIR=$vc_out"

sha256sum --status -c manifests/v01_legacy_code.sha256
sha256sum --status -c qc/v01_existing_responses/SHA256SUMS
sha256sum --status -c qc/v01_response_provenance/SHA256SUMS
sha256sum --status -c qc/v01_response_matrix_validation/SHA256SUMS

python3 -u "$vc_script" \
  --legacy-root "$HOME/VC_reliability" \
  --response-checksums qc/v01_existing_responses/legacy_inputs.sha256 \
  --provenance qc/v01_response_provenance/response_provenance.json \
  --output-dir "$vc_out"

git rev-parse HEAD > "$vc_out/code_commit.txt"
sha256sum -- "$vc_out"/{summary.json,parameters.json,predictions.h5,code_commit.txt} \
  "$vc_out"/internal_{model_metrics,task_errors,risk_metrics}.tsv \
  "$vc_script" "$vc_job" manifests/v01_legacy_code.sha256 > "$vc_out/SHA256SUMS"
sha256sum --status -c "$vc_out/SHA256SUMS"
echo "FOUR_SOURCE_CALIBRATION_COMPLETE $SLURM_JOB_ID"
