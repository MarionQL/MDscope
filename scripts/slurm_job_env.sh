# shellcheck shell=bash
# Job-start environment shared by the MDscope Slurm batch scripts; sourced
# from the submission checkout. Jobs never install MDscope: install it once
# into the environment selected here. See docs/workflows.md#slurm-jobs.

# MDSCOPE_ENV_SETUP: unset/empty = ARC example (Miniconda env "mdanalysis"),
# "none" = use the environment inherited from sbatch, otherwise a bash file
# to source (for example, one that loads modules and activates an env).
mdscope_activate_env() {
    local setup="${MDSCOPE_ENV_SETUP:-}"
    if [[ "${setup}" == "none" ]]; then
        echo "Environment:  inherited from submission (MDSCOPE_ENV_SETUP=none)"
    elif [[ -z "${setup}" ]]; then
        local conda_sh=/apps/common/software/Miniconda3/24.7.1-0/etc/profile.d/conda.sh
        if [[ ! -f "${conda_sh}" ]]; then
            echo "ERROR: ARC default environment not found (${conda_sh})." >&2
            echo "Set MDSCOPE_ENV_SETUP to an activation script, or to 'none'" \
                 "to use the submitting environment (see docs/workflows.md#slurm-jobs)." >&2
            return 2
        fi
        # shellcheck source=/dev/null
        source "${conda_sh}"
        conda activate mdanalysis
        echo "Environment:  ARC conda env mdanalysis"
    elif [[ -f "${setup}" ]]; then
        # shellcheck source=/dev/null
        source "${setup}"
        echo "Environment:  sourced ${setup}"
    else
        echo "ERROR: MDSCOPE_ENV_SETUP file not found: ${setup}" >&2
        return 2
    fi
}

# Verify the activated python imports MDscope; report version and source.
mdscope_preflight() {
    local info
    if ! info="$(python -c 'import pathlib, platform, mdscope, mdscope.cli
print(mdscope.__version__, platform.python_version(), pathlib.Path(mdscope.__file__).parent)')"; then
        cat >&2 <<'EOF'
ERROR: MDscope is not importable by `python` in the job environment.
Jobs do not install MDscope. Before submitting, activate the job environment
and install once from the checkout:  python -m pip install -e '.[sasa]'
Select that environment with MDSCOPE_ENV_SETUP (see docs/workflows.md#slurm-jobs).
EOF
        return 3
    fi
    local version python_version source_dir
    read -r version python_version source_dir <<<"${info}"
    echo "MDscope:      ${version} (Python ${python_version})"
    echo "MDscope src:  ${source_dir}"
    if git -C "${source_dir}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        echo "MDscope commit: $(git -C "${source_dir}" rev-parse --short HEAD)"
    fi
}
