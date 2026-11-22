#!/usr/bin/env bash
# Run the real GROMACS test-data checks locally, mirroring the manual GitHub
# workflows, against trajectories kept outside git (default: examples/data/).
#
#   scripts/run_real_data_checks.sh CHECK
#   CHECK: smoke | hbond | solvation | performance | finalization |
#          box-sweep | box-clearance | all
#
# Data: the four test-data-v1 TPR/XTC files, verified against
# tests/data/real_test_data.sha256. For quicker iteration, point
# MDSCOPE_REAL_DATA at a strided copy from scripts/stride_real_data.py (e.g.
# examples/data/stride5); it is verified against its own SHA256SUMS, runs only
# smoke/hbond, and is not CI-equivalent evidence. Outputs: real_test_results/
# (git-ignored).
# Activate an environment with `pip install -e '.[sasa]'` first: the smoke
# check validates the production YAML, which enables SASA (FreeSASA).
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${root}"
check="${1:-smoke}"
export MDSCOPE_REAL_DATA="${MDSCOPE_REAL_DATA:-examples/data}"
export MDSCOPE_REAL_OUTPUT="${MDSCOPE_REAL_OUTPUT:-real_test_results}"
export PYTHONUNBUFFERED=1
export MDSCOPE_PROFILE="${MDSCOPE_PROFILE:-1}"

case "${check}" in
    smoke|hbond|solvation|performance|finalization|box-sweep|box-clearance|all) ;;
    *) echo "usage: $0 [smoke|hbond|solvation|performance|finalization|box-sweep|box-clearance|all]" >&2
       exit 2 ;;
esac

if [[ ! -d "${MDSCOPE_REAL_DATA}" ]]; then
    echo "ERROR: real test data not found in ${MDSCOPE_REAL_DATA}/." >&2
    echo "Copy the four test-data-v1 TPR/XTC files there, or set MDSCOPE_REAL_DATA." >&2
    exit 2
fi
manifest="${root}/tests/data/real_test_data.sha256"
derived=""
if [[ -f "${MDSCOPE_REAL_DATA}/provenance.json" ]]; then
    derived="$(python -c 'import json, sys; print(json.load(open(sys.argv[1]))["stride"])' \
        "${MDSCOPE_REAL_DATA}/provenance.json")"
    manifest="$(cd "${MDSCOPE_REAL_DATA}" && pwd)/SHA256SUMS"
    case "${check}" in
        solvation|performance|finalization|box-sweep|box-clearance|all)
            echo "ERROR: ${check} needs more frames than the stride-${derived} copy;" \
                 "use the full files (MDSCOPE_REAL_DATA=examples/data)." >&2
            exit 2 ;;
    esac
    echo "WARNING: stride-${derived} derived data: quick smoke/parity check only," \
         "not CI-equivalent evidence." >&2
fi
if command -v sha256sum >/dev/null 2>&1; then
    hasher=(sha256sum -c)
else
    hasher=(shasum -a 256 -c)  # macOS
fi
echo "Verifying ${MDSCOPE_REAL_DATA}/ against ${manifest#"${root}/"}"
(cd "${MDSCOPE_REAL_DATA}" && "${hasher[@]}" "${manifest}")
python -c "import mdscope, sys; print(f'MDscope {mdscope.__version__} from {mdscope.__file__} (Python {sys.version.split()[0]})')"

run() {
    echo "==> $*"
    "$@"
}

smoke() {
    run python tests/real_hs_topology_audit.py
    run python tests/real_full_config_validate.py
    run python tests/real_trajectory_smoke.py
    run python tests/real_sdc4_mechanisms.py
    run python tests/real_motif_shared_states.py
    run python tests/real_glycan_hydration.py
    run python tests/real_trajectory_physics.py
}

hbond() {
    MDSCOPE_HBOND_BENCH_FRAMES="${MDSCOPE_HBOND_BENCH_FRAMES:-8}" run python tests/real_hbond_benchmark.py
}

solvation() {
    run python tests/real_release_solvation_bench.py --data-dir "${MDSCOPE_REAL_DATA}" \
        --output-dir "${MDSCOPE_REAL_OUTPUT}/solvation_benchmark" \
        --standard-frames 24 --kinetics-frames 12
}

performance() {
    MDSCOPE_REAL_OUTPUT="${MDSCOPE_REAL_OUTPUT}/benchmark" MDSCOPE_STATS_WORKERS=1 \
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
        run python tests/performance_pipeline_benchmark.py --frames 40 160 --draws 100 500
}

finalization() {
    run python tests/real_finalization_benchmark.py
}

conditions=(carb nocarb)

box-sweep() {
    # Hypothetical periodic-box lattice sensitivity (not a box-size certification).
    for condition in "${conditions[@]}"; do
        run python scripts/box_size_sweep.py             --tpr "${MDSCOPE_REAL_DATA}/${condition}-sdc4-step7_production.tpr"             --xtc "${MDSCOPE_REAL_DATA}/${condition}-sdc4-rep1_test.xtc"             --out "${MDSCOPE_REAL_OUTPUT}/box_size_sweep/${condition}"             --max-frames "${MDSCOPE_BOX_FRAMES:-40}" --lateral 26 28 30 32 --vertical 28
    done
}

box-clearance() {
    # Syndecan-HS periodic-image clearance in the released trajectories.
    for condition in "${conditions[@]}"; do
        run python scripts/release_box_clearance.py             --tpr "${MDSCOPE_REAL_DATA}/${condition}-sdc4-step7_production.tpr"             --xtc "${MDSCOPE_REAL_DATA}/${condition}-sdc4-rep1_test.xtc"             --max-frames "${MDSCOPE_BOX_FRAMES:-200}"             --out "${MDSCOPE_REAL_OUTPUT}/box_clearance/${condition}"
    done
}

if [[ "${check}" == all ]]; then
    smoke; hbond; solvation; performance; finalization; box-sweep; box-clearance
else
    "${check}"
fi
echo "PASS: ${check}${derived:+ (stride-${derived} derived data)} (outputs in ${MDSCOPE_REAL_OUTPUT}/)"
