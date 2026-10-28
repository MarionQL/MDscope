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

# SLURM keeps stdout and stderr in separate files. Batch progress writes to
# stderr; combine streams so progress and errors are visible in the .out log.
exec 2>&1

# Run from the directory with YAML
ROOT="${SLURM_SUBMIT_DIR}"
YAML="${1:-${ROOT}/synd4_full_mechanisms.yaml}"

mkdir -p "${ROOT}/slurm-out/mdscope"
cd "${ROOT}"

# Activate the environment containing MDscope and its dependencies.
source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate mdanalysis
export PYTHONUNBUFFERED=1
export MDSCOPE_PROGRESS_EVERY_PERCENT="${MDSCOPE_PROGRESS_EVERY_PERCENT:-1}"
# Enable lightweight trajectory/analysis timing during SLURM runs by default.
export MDSCOPE_PROFILE="${MDSCOPE_PROFILE:-1}"
MDSCOPE_ROOT="/projects/bevanlab/students/marionlopresti/hepsulf_MD"
export PYTHONPATH="${MDSCOPE_ROOT}:${PYTHONPATH:-}"

python - <<'PY'
import MDscope.analysis.cli as cli
import MDscope.analysis.progress as progress
print(f"MDscope CLI source: {cli.__file__}", flush=True)
print(f"Progress reporter source: {progress.__file__}", flush=True)
PY
if git -C "${ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "MDscope commit: $(git -C "${ROOT}" rev-parse --short HEAD)"
fi

echo "============================================================"
echo "MDscope replicate array"
echo "Job ID:       ${SLURM_ARRAY_JOB_ID}"
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Config:       ${YAML}"
echo "Host:         $(hostname)"
echo "Started:      $(date)"
echo "============================================================"

python -u -m MDscope.analysis.cli replicate \
    "${YAML}" \
    --index "${SLURM_ARRAY_TASK_ID}"

echo "============================================================"
echo "Finished task ${SLURM_ARRAY_TASK_ID}: $(date)"
echo "============================================================"
