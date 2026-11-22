"""Slurm batch scripts: array cardinality, environment hook, no-install startup.

Slurm is simulated: scripts run under bash with SLURM_* variables, and a
``python`` shim records every call, forwarding all but the analysis entry
points (``python -u -m mdscope.cli ...``) to the test interpreter.
"""
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

import mdscope

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = ("slurm/run_mdscope.sh", "slurm/finalize_mdscope.sh", "slurm/plotting_mdscope.sh",
           "slurm/statistics_mdscope.sh")
HELPER = "slurm/job_env.sh"
N_REPLICATES = 4  # Two systems x two replicates.
# A pip install command line (instructions inside messages are allowed).
INSTALL = re.compile(r"^\s*(python3?\s+-m\s+)?pip3?\s+install", re.M)


def _bash():
    if os.name != "nt":
        return shutil.which("bash")
    # PATH's bash.exe on Windows may be WSL; use the Git for Windows bash.
    git = Path(shutil.which("git")).resolve()
    bash = next(p / "bin" / "bash.exe" for p in git.parents
                if (p / "bin" / "bash.exe").is_file())
    return str(bash)


def _shim(path, body):
    path.write_text("#!/bin/bash\n" + body, newline="\n")
    path.chmod(0o755)


@pytest.fixture
def job(tmp_path, project):
    """Return run(script, yaml_path=..., **env) -> (result, python calls)."""
    checkout = tmp_path / "checkout"
    (checkout / "slurm").mkdir(parents=True)
    for name in (*SCRIPTS, HELPER):
        shutil.copyfile(REPO / name, checkout / name)
    config, _ = project("slurm")
    raw = yaml.safe_load(config.read_text())
    tiny = raw["systems"]["tiny"]
    tiny["replicates"] = [{**tiny["replicates"][0], "label": f"rep{i}"} for i in (1, 2)]
    raw["systems"]["tiny_b"] = tiny
    config.write_text(yaml.safe_dump(raw, sort_keys=False))

    stubs = tmp_path / "bin"
    stubs.mkdir()
    log = tmp_path / "python_calls.log"
    real = Path(sys.executable).as_posix()
    _shim(stubs / "python", f'printf "%s\\n" "$*" >> "{log.as_posix()}"\n'
          'if [[ "$*" == "-u -m mdscope.cli "* ]]; then exit 0; fi\n'
          f'exec "{real}" "$@"\n')

    def run(script, yaml_path=config, **env):
        environment = {k: v for k, v in os.environ.items()
                       if not k.startswith(("SLURM_", "MDSCOPE_"))}
        environment.update(PATH=f"{stubs}{os.pathsep}{environment['PATH']}",
                           SLURM_SUBMIT_DIR=checkout.as_posix(),
                           SLURM_JOB_ID="77", **env)
        log.write_text("")
        result = subprocess.run(
            [_bash(), (checkout / script).as_posix(), Path(yaml_path).as_posix()],
            cwd=checkout, env=environment, text=True, capture_output=True,
            timeout=300, check=False)
        return result, log.read_text().splitlines()

    run.config = config
    run.real_python = real
    return run


def _array(task_min, task_max, count, task_id=None):
    return dict(SLURM_ARRAY_JOB_ID="77", SLURM_ARRAY_TASK_MIN=str(task_min),
                SLURM_ARRAY_TASK_MAX=str(task_max), SLURM_ARRAY_TASK_COUNT=str(count),
                SLURM_ARRAY_TASK_ID=str(task_min if task_id is None else task_id))


def _no_install(calls, output):
    assert not any("pip" in call for call in calls), calls
    assert "pip install" not in output


@pytest.mark.parametrize("name", (*SCRIPTS, HELPER))
def test_scripts_parse_and_never_install(name):
    path = REPO / name
    assert b"\r\n" not in path.read_bytes()  # sbatch rejects DOS line endings.
    assert not INSTALL.search(path.read_text())
    result = subprocess.run([_bash(), "-n", path.as_posix()], capture_output=True,
                            text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_plan_replicate_count_is_the_array_cardinality(job):
    result = subprocess.run(
        [sys.executable, "-m", "mdscope.cli", "plan", str(job.config), "--format", "json"],
        capture_output=True, text=True, check=True)
    assert '"replicates"' in result.stdout
    assert len(json.loads(result.stdout)["replicates"]) == N_REPLICATES


@pytest.mark.parametrize("task_min,task_max,count", [
    (0, 2, 3),  # smaller than the YAML replica count
    (0, 5, 6),  # larger (the ARC six-replica default)
    (1, 4, 4),  # right cardinality, wrong offset
    (0, 4, 4),  # sparse, e.g. --array=0,1,2,4
])
def test_mismatched_array_fails_before_analysis(job, task_min, task_max, count):
    result, calls = job("slurm/run_mdscope.sh", MDSCOPE_ENV_SETUP="none",
                        **_array(task_min, task_max, count))
    assert result.returncode == 2, result.stdout + result.stderr
    assert f"expected --array=0-{N_REPLICATES - 1}" in result.stdout
    assert f"sbatch --array=0-{N_REPLICATES - 1} slurm/run_mdscope.sh" in result.stdout
    assert not any("mdscope.cli replicate" in call for call in calls), calls
    _no_install(calls, result.stdout)


def test_non_array_submission_fails_before_analysis(job):
    result, calls = job("slurm/run_mdscope.sh", MDSCOPE_ENV_SETUP="none")
    assert result.returncode == 2
    assert f"expected --array=0-{N_REPLICATES - 1}" in result.stdout
    assert not any("mdscope.cli replicate" in call for call in calls), calls


def test_matching_array_runs_task_without_install(job):
    result, calls = job("slurm/run_mdscope.sh", MDSCOPE_ENV_SETUP="none",
                        **_array(0, N_REPLICATES - 1, N_REPLICATES, task_id=2))
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"MDscope:      {mdscope.__version__}" in result.stdout
    assert f"MDscope src:  {Path(mdscope.__file__).parent}" in result.stdout
    assert calls[-1] == f"-u -m mdscope.cli replicate {job.config.as_posix()} --index 2"
    _no_install(calls, result.stdout)


@pytest.mark.parametrize("script,command", [
    ("slurm/finalize_mdscope.sh", "finalize"),
    ("slurm/statistics_mdscope.sh", "statistics"),
    ("slurm/plotting_mdscope.sh", "plotting"),
])
def test_single_jobs_start_without_install(job, script, command):
    result, calls = job(script, MDSCOPE_ENV_SETUP="none")
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"MDscope:      {mdscope.__version__}" in result.stdout
    assert calls[-1] == f"-u -m mdscope.cli {command} {job.config.as_posix()}"
    _no_install(calls, result.stdout)


def test_env_setup_file_is_sourced(job, tmp_path):
    setup = tmp_path / "site_env.sh"
    setup.write_text('echo "site modules loaded"\n', newline="\n")
    result, calls = job("slurm/finalize_mdscope.sh", MDSCOPE_ENV_SETUP=setup.as_posix())
    assert result.returncode == 0, result.stdout + result.stderr
    assert "site modules loaded" in result.stdout
    assert calls[-1].startswith("-u -m mdscope.cli finalize")


@pytest.mark.parametrize("script", SCRIPTS)
def test_unimportable_mdscope_fails_with_setup_instructions(job, tmp_path, script):
    setup = tmp_path / "broken_env.sh"
    # A job environment whose python lacks MDscope (no site-packages/PYTHONPATH).
    setup.write_text(f'python() {{ env -u PYTHONPATH "{job.real_python}" -S "$@"; }}\n',
                     newline="\n")
    result, calls = job(script, MDSCOPE_ENV_SETUP=setup.as_posix(),
                        **_array(0, N_REPLICATES - 1, N_REPLICATES))
    assert result.returncode == 3, result.stdout + result.stderr
    output = result.stdout + result.stderr
    assert "MDscope is not importable" in output
    assert "python -m pip install -e '.[sasa]'" in output
    assert "MDSCOPE_ENV_SETUP" in output
    assert calls == []


@pytest.mark.parametrize("value,message", [
    (None, "ARC default environment not found"),
    ("missing_env.sh", "MDSCOPE_ENV_SETUP file not found"),
])
def test_unavailable_environment_fails_before_python(job, value, message):
    env = {} if value is None else {"MDSCOPE_ENV_SETUP": value}
    result, calls = job("slurm/statistics_mdscope.sh", **env)
    assert result.returncode == 2
    assert message in result.stdout + result.stderr
    assert calls == []
