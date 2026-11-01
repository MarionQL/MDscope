#!/bin/bash
#SBATCH --job-name=mdscope_finalize
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH -A bevanlab
#SBATCH --output=slurm-out/mdscope/%x-%j.out
#SBATCH --error=slurm-out/mdscope/%x-%j.err

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
YAML="${1:-${ROOT}/synd4_full_mechanisms.yaml}"

mkdir -p "${ROOT}/slurm-out/mdscope"
cd "${ROOT}"
if [[ ! -f "${YAML}" ]]; then
    echo "ERROR: YAML not found: ${YAML}" >&2
    exit 2
fi

# Finalization reads completed replicate outputs and generates the
# cross-replicate statistics and figures; it does not rerun trajectories.
source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate mdanalysis
export PYTHONUNBUFFERED=1
python -m pip install --no-deps -e "${ROOT}"

echo "============================================================"
echo "MDscope finalize"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Config:       ${YAML}"
echo "Host:         $(hostname)"
echo "Started:      $(date)"
echo "============================================================"

python -u -m mdscope.cli finalize "${YAML}"

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
