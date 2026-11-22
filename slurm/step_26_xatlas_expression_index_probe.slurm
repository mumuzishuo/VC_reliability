#!/bin/bash -l
#SBATCH --job-name=vc_xatlas_obs
#SBATCH --account=zhenwei01_bio214
#SBATCH --partition=cpu8358
#SBATCH --qos=26cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=logs/v03_xatlas_expression_index_probe/step26_%j.out
#SBATCH --error=logs/v03_xatlas_expression_index_probe/step26_%j.err

set -euo pipefail
cd "$HOME/VC_reliability_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd

out=qc/v03_xatlas_expression_index_probe
script=scripts/external_metadata/08_probe_xatlas_expression_index.py
job=slurm/step_26_xatlas_expression_index_probe.slurm
figshare=manifests/v01_xatlas_figshare_article.sha256
hct_qc=qc/v01_xatlas_hct116_observed
hek_qc=qc/v01_xatlas_hek293t_observed
guide_inputs=manifests/v01_xatlas_hct116_observed_inputs.sha256
hek_input=manifests/v01_xatlas_hek293t_guide_calls.sha256

sha256sum --status -c "$figshare" "$guide_inputs" "$hek_input" \
    "$hct_qc/SHA256SUMS" "$hek_qc/SHA256SUMS"
python3 -u "$script" "$out"

git rev-parse HEAD > "$out/code_commit.txt"
sha256sum -- "$out/summary.json" "$out/code_commit.txt" "$script" "$job" \
    "$figshare" "$guide_inputs" "$hek_input" \
    "$hct_qc/summary.json" "$hek_qc/summary.json" \
    "$hct_qc/SHA256SUMS" "$hek_qc/SHA256SUMS" > "$out/SHA256SUMS"
sha256sum --status -c "$out/SHA256SUMS"
echo "XATLAS_OBS_INDEX_PROBE_OK $SLURM_JOB_ID"
