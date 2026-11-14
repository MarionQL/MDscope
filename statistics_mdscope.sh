#!/bin/bash
#SBATCH --job-name=mdscope_statistics
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH -A bevanlab
#SBATCH --output=slurm-out/mdscope/%x-%j.out
#SBATCH --error=slurm-out/mdscope/%x-%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
YAML="${1:-${ROOT}/synd4_full_mechanisms_finalize_fast.yaml}"
cd "${ROOT}"
echo "Statistics setup started: $(date) | config: ${YAML} | requested CPUs: ${MDSCOPE_STATS_WORKERS:-${SLURM_CPUS_PER_TASK:-default}}"
if [[ ! -f "${YAML}" ]]; then
  echo "Missing config: ${YAML}" >&2
  exit 2
fi
source /apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
conda activate mdanalysis
export PYTHONUNBUFFERED=1
python -m pip install --no-deps -e "${ROOT}"
echo "Statistics started: $(date)"
python -u -m mdscope.cli statistics "${YAML}"
echo "Statistics finished: $(date)"
