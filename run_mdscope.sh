#!/bin/bash
#SBATCH --job-name=mdscope_synd4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH -A bevanlab
#SBATCH --array=0-5
#SBATCH --output=slurm-out/mdscope/%x-%A_%a.out
#SBATCH --error=slurm-out/mdscope/%x-%A_%a.err

set -euo pipefail

# Run from the directory with YAML
ROOT="${SLURM_SUBMIT_DIR}"
YAML="${1:-${ROOT}/synd4_baseline.yaml}"

mkdir -p "${ROOT}/slurm-out/mdscope"
cd "${ROOT}"

# Activate the environment containing MDscope and its dependencies.
source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate mdanalysis
MDSCOPE_ROOT="/projects/bevanlab/students/marionlopresti/hepsulf_MD"
export PYTHONPATH="${MDSCOPE_ROOT}:${PYTHONPATH:-}"

echo "============================================================"
echo "MDscope replicate array"
echo "Job ID:       ${SLURM_ARRAY_JOB_ID}"
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Config:       ${YAML}"
echo "Host:         $(hostname)"
echo "Started:      $(date)"
echo "============================================================"

python -m MDscope.analysis.cli replicate \
    "${YAML}" \
    --index "${SLURM_ARRAY_TASK_ID}"

echo "============================================================"
echo "Finished task ${SLURM_ARRAY_TASK_ID}: $(date)"
echo "============================================================"
