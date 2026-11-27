#!/bin/bash -l
#SBATCH --job-name=vc_xatlas_prepare
#SBATCH --account=zhenwei01_bio214
#SBATCH --partition=cpu8358
#SBATCH --qos=26cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/v01_xatlas_prediction_candidate/step30_%j.out
#SBATCH --error=logs/v01_xatlas_prediction_candidate/step30_%j.err

set -euo pipefail
cd "$HOME/VC_reliability_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

vc_script=scripts/development/02_prepare_xatlas_external.py
vc_job=slurm/step_30_prepare_xatlas_external.slurm
vc_out="results/v01_xatlas_prediction_candidate/job_${SLURM_JOB_ID}"
echo "OUTPUT_DIR=$vc_out"

sha256sum --quiet -c manifests/v01_legacy_code.sha256
python3 -u "$vc_script" --output-dir "$vc_out"

git rev-parse HEAD > "$vc_out/code_commit.txt"
sha256sum -- "$vc_out"/{task_eligibility.tsv,feature_axis.tsv,predictions.h5,parameters.json,summary.json,code_commit.txt} \
  "$vc_script" "$vc_job" manifests/v01_legacy_code.sha256 > "$vc_out/SHA256SUMS"
sha256sum --quiet -c "$vc_out/SHA256SUMS"
echo "XATLAS_PREDICTION_CANDIDATE_COMPLETE $SLURM_JOB_ID"
