#!/bin/bash
#SBATCH --job-name=mdscope_plotting
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH -p normal_q
#SBATCH -A bevanlab
#SBATCH --output=slurm-out/mdscope/%x-%j.out
#SBATCH --error=slurm-out/mdscope/%x-%j.err
# ARC defaults; override on the sbatch command line (docs/workflows.md#slurm-jobs).

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
YAML="${1:-${ROOT}/synd4_full_mechanisms_finalize_fast.yaml}"
cd "${ROOT}"
if [[ ! -f "${YAML}" ]]; then
  echo "Missing config: ${YAML}" >&2
  exit 2
fi
export PYTHONUNBUFFERED=1
# shellcheck source=slurm/job_env.sh
source "${ROOT}/slurm/job_env.sh"
mdscope_activate_env
mdscope_preflight
echo "Plotting started: $(date)"
python -u -m mdscope.cli plotting "${YAML}"
echo "Plotting finished: $(date)"
