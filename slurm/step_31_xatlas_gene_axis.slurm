#!/bin/bash -l
#SBATCH --job-name=vc_xatlas_var
#SBATCH --account=zhenwei01_bio214
#SBATCH --partition=cpu8358
#SBATCH --qos=26cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/v01_xatlas_gene_axis/step31_%j.out
#SBATCH --error=logs/v01_xatlas_gene_axis/step31_%j.err

set -euo pipefail
cd "$HOME/VC_reliability_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd

vc_cell_line=${1:-HCT116}
vc_out="qc/v01_xatlas_gene_axis/job_$SLURM_JOB_ID"
vc_candidate=results/v01_xatlas_prediction_candidate/job_3239864
vc_script=scripts/external_metadata/11_probe_xatlas_gene_axis.py
vc_range=scripts/external_metadata/08_probe_xatlas_expression_index.py
vc_job=slurm/step_31_xatlas_gene_axis.slurm
vc_catalog=manifests/v01_xatlas_figshare_article.sha256

echo "OUTPUT_DIR=$vc_out"
sha256sum --quiet -c "$vc_candidate/SHA256SUMS" "$vc_catalog" manifests/v01_legacy_code.sha256
python3 -u "$vc_script" \
    --cell-line "$vc_cell_line" \
    --features "$vc_candidate/feature_axis.tsv" \
    --output-dir "$vc_out"

git rev-parse HEAD > "$vc_out/code_commit.txt"
sha256sum -- "$vc_out"/{feature_h5ad_mapping.tsv,summary.json,code_commit.txt} \
    "$vc_script" "$vc_range" "$vc_job" "$vc_catalog" \
    manifests/v01_legacy_code.sha256 "$vc_candidate/feature_axis.tsv" \
    > "$vc_out/SHA256SUMS"
sha256sum --quiet -c "$vc_out/SHA256SUMS"
echo "XATLAS_VAR_PROBE_COMPLETE $SLURM_JOB_ID"
