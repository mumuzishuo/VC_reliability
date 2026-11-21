#!/bin/bash -l
#SBATCH --job-name=vc_risk_audit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/v01_risk_alignment/step12_%j.out
#SBATCH --error=logs/v01_risk_alignment/step12_%j.err

set -euo pipefail
vc_repo="$HOME/VC_reliability_repo"
vc_old="$HOME/VC_reliability"
cd "$vc_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

vc_out=qc/v01_risk_alignment
if [ -e "$vc_out" ]; then
    echo "Risk audit output exists: inspect it before rerunning; refusing overwrite" >&2
    exit 1
fi
sha256sum --status -c qc/v01_response_matrix_validation/SHA256SUMS

python3 -u scripts/data_audit/11_audit_risk_alignment.py \
    --per-task "$vc_old/results/baseline_v01/prediction_level_results.tsv" \
    --risk-table "$vc_old/results/baseline_v01/risk_predictions.tsv" \
    --frozen-protocol "$vc_old/results/baseline_v01/frozen_protocol.json" \
    --baseline-code "$vc_old/vc_baseline.py" \
    --prediction-h5 "$vc_old/results/baseline_v01/predictions.h5" \
    --output-dir "$vc_out"

git rev-parse HEAD > "$vc_out/code_commit.txt"
python3 --version > "$vc_out/python_version.txt" 2>&1
scontrol show job "$SLURM_JOB_ID" > "$vc_out/slurm_job.txt"
(cd "$vc_old" && sha256sum -- \
    vc_baseline.py \
    results/baseline_v01/prediction_level_results.tsv \
    results/baseline_v01/risk_predictions.tsv \
    results/baseline_v01/frozen_protocol.json \
    results/baseline_v01/predictions.h5) > "$vc_out/legacy_inputs.sha256"

sha256sum -- \
    "$vc_out/risk_alignment_summary.tsv" \
    "$vc_out/risk_alignment_summary.json" \
    "$vc_out/legacy_inputs.sha256" \
    "$vc_out/code_commit.txt" \
    "$vc_out/python_version.txt" \
    "$vc_out/slurm_job.txt" \
    scripts/data_audit/11_audit_risk_alignment.py \
    slurm/step_12_risk_alignment.slurm \
    > "$vc_out/SHA256SUMS"
sha256sum --status -c "$vc_out/SHA256SUMS"
echo "RISK_ALIGNMENT_QC_COMPLETE $SLURM_JOB_ID"
