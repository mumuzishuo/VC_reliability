#!/bin/bash -l
#SBATCH --job-name=vc_hepg2_raw_qc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/v01_hepg2_raw_recomputation/step13_%j.out
#SBATCH --error=logs/v01_hepg2_raw_recomputation/step13_%j.err

set -euo pipefail
vc_repo="$HOME/VC_reliability_repo"
vc_old="$HOME/VC_reliability"
cd "$vc_repo"
module load anaconda3/2022.10-gcc-8.5.0-4dp3trd
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

vc_out=qc/v01_hepg2_raw_recomputation
vc_script=scripts/data_audit/12_recompute_hepg2_raw_counts.py
vc_job=slurm/step_13_recompute_hepg2_raw_counts.slurm
vc_manifest="$vc_old/processed/v01/manifest.json"
vc_target="$vc_old/processed/v01/evaluation_only/HepG2_responses.h5"
vc_preprocess="$vc_old/vc_preprocess.py"

if [ -e "$vc_out" ]; then
    echo "Audit output exists; inspect it and use a new version before rerunning" >&2
    exit 1
fi
sha256sum --status -c qc/v01_response_matrix_validation/SHA256SUMS
sha256sum --status -c qc/v01_risk_alignment/SHA256SUMS

# This writes and hashes the fixed labels before reading raw X or target delta.
python3 -u "$vc_script" prepare \
    --manifest "$vc_manifest" \
    --response "$vc_target" \
    --preprocess-code "$vc_preprocess" \
    --output-dir "$vc_out"
(cd "$vc_out" && sha256sum -- selection.tsv selection.json > SHA256SUMS.prepare)
(cd "$vc_out" && sha256sum --status -c SHA256SUMS.prepare)

set +e
python3 -u "$vc_script" audit \
    --manifest "$vc_manifest" \
    --response "$vc_target" \
    --preprocess-code "$vc_preprocess" \
    --output-dir "$vc_out"
vc_audit_exit=$?
set -e

git rev-parse HEAD > "$vc_out/code_commit.txt"
python3 --version > "$vc_out/python_version.txt" 2>&1
scontrol show job "$SLURM_JOB_ID" > "$vc_out/slurm_job.txt"
(cd "$vc_old" && sha256sum -- \
    raw/GSE264667_hepg2_raw_singlecell_01.h5ad \
    processed/v01/manifest.json \
    processed/v01/evaluation_only/HepG2_responses.h5 \
    vc_preprocess.py) > "$vc_out/input_sha256.txt"

vc_files=("$vc_out/selection.tsv" "$vc_out/selection.json"
          "$vc_out/SHA256SUMS.prepare" "$vc_out/input_sha256.txt"
          "$vc_out/code_commit.txt" "$vc_out/python_version.txt"
          "$vc_out/slurm_job.txt" "$vc_script" "$vc_job")
if [ -f "$vc_out/audit_summary.json" ] && [ -f "$vc_out/per_task.tsv" ]; then
    vc_files+=("$vc_out/audit_summary.json" "$vc_out/per_task.tsv")
fi
sha256sum -- "${vc_files[@]}" > "$vc_out/SHA256SUMS"
sha256sum --status -c "$vc_out/SHA256SUMS"
if [ "$vc_audit_exit" -ne 0 ]; then
    echo "HepG2 audit failed (exit $vc_audit_exit); inspect diagnostics and log" >&2
    exit "$vc_audit_exit"
fi
echo "HEPG2_RAW_QC_COMPLETE $SLURM_JOB_ID"
