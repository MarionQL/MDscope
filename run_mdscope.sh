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
# ARC defaults above suit the six-replica synd4_full_mechanisms.yaml; override
# them on the sbatch command line (docs/workflows.md#slurm-jobs).

set -euo pipefail

# SLURM keeps stdout and stderr in separate files. Batch progress writes to
# stderr; combine streams so progress and errors are visible in the .out log.
exec 2>&1

# Run from the directory with YAML
ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
YAML="${1:-${ROOT}/synd4_full_mechanisms.yaml}"

mkdir -p "${ROOT}/slurm-out/mdscope"
cd "${ROOT}"
if [[ ! -f "${YAML}" ]]; then
    echo "ERROR: YAML not found: ${YAML}" >&2
    exit 2
fi

export PYTHONUNBUFFERED=1
export MDSCOPE_PROGRESS_EVERY_PERCENT="${MDSCOPE_PROGRESS_EVERY_PERCENT:-1}"
# Enable lightweight trajectory/analysis timing during SLURM runs by default.
export MDSCOPE_PROFILE="${MDSCOPE_PROFILE:-1}"
# shellcheck source=scripts/slurm_job_env.sh
source "${ROOT}/scripts/slurm_job_env.sh"
mdscope_activate_env
mdscope_preflight

# One task per YAML replicate, indexed 0..N-1 in `mdscope plan` order. Every
# task checks the whole submitted array before reading any trajectory.
if ! n_replicates="$(python -m mdscope.cli plan "${YAML}" --format json \
        | python -c 'import json, sys; print(len(json.load(sys.stdin)["replicates"]))')"; then
    echo "ERROR: could not derive the replicate count from: mdscope plan ${YAML}" >&2
    exit 2
fi
expected="0-$((n_replicates - 1))"
if [[ "${SLURM_ARRAY_TASK_COUNT:-}" != "${n_replicates}" \
      || "${SLURM_ARRAY_TASK_MIN:-}" != 0 \
      || "${SLURM_ARRAY_TASK_MAX:-}" != "$((n_replicates - 1))" ]]; then
    echo "ERROR: Slurm array ${SLURM_ARRAY_TASK_MIN:-?}-${SLURM_ARRAY_TASK_MAX:-?}" \
         "(${SLURM_ARRAY_TASK_COUNT:-?} tasks) does not match the ${n_replicates}" \
         "replicates in ${YAML}; expected --array=${expected}." >&2
    echo "Resubmit: sbatch --array=${expected} run_mdscope.sh ${YAML}" >&2
    exit 2
fi

echo "============================================================"
echo "MDscope replicate array"
echo "Job ID:       ${SLURM_ARRAY_JOB_ID}"
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Config:       ${YAML}"
echo "Host:         $(hostname)"
echo "Started:      $(date)"
echo "============================================================"

python -u -m mdscope.cli replicate \
    "${YAML}" \
    --index "${SLURM_ARRAY_TASK_ID}"

echo "============================================================"
echo "Finished task ${SLURM_ARRAY_TASK_ID}: $(date)"
echo "============================================================"
