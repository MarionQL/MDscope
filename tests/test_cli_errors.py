"""CLI contract: expected input mistakes are one-line errors (exit 2); faults keep tracebacks."""
import subprocess
import sys

import pytest
import yaml

from mdscope_test_support import invoke_cli as invoke


def _user_error(result, *fragments):
    assert result.returncode == 2, result.stderr
    assert result.stderr.startswith("mdscope: error: "), result.stderr
    assert "Traceback" not in result.stderr
    for fragment in fragments:
        assert fragment in result.stderr


def _edit(path, change):
    raw = yaml.safe_load(path.read_text()); change(raw); path.write_text(yaml.safe_dump(raw))


def test_config_validation_error_is_formatted(tmp_path, project):
    path, root = project("bad_option")
    _edit(path, lambda raw: raw["trajectory"].update(stride=0))
    _user_error(invoke(tmp_path, "validate", path, ok=False), "bad_option.yaml", "trajectory.stride must be >= 1")
    assert not root.exists()


def test_malformed_yaml_and_missing_config_are_formatted(tmp_path):
    (tmp_path / "broken.yaml").write_text("systems: [unclosed\n")
    _user_error(invoke(tmp_path, "validate", "broken.yaml", ok=False), "invalid YAML")
    _user_error(invoke(tmp_path, "run", "absent.yaml", ok=False), "Configuration file does not exist: absent.yaml")


@pytest.mark.parametrize("command", [("validate",), ("replicate", "--index", "0")])
def test_missing_trajectory_is_formatted(tmp_path, project, command):
    path, _ = project("missing_traj")
    _edit(path, lambda raw: raw["systems"]["tiny"]["replicates"][0].update(trajectory="gone.dcd"))
    _user_error(invoke(tmp_path, command[0], path, *command[1:], ok=False),
                "System 'tiny' replicate 'rep1' trajectory does not exist", "gone.dcd")


@pytest.mark.parametrize("change,fragment", [
    (lambda raw: raw["selections"].update(features="name CA and and"), "System 'tiny': cannot resolve"),
    (lambda raw: raw["analyses"].update(bilayer_normal={"enabled": True, "species_selections": {"bad": "resname and"}}),
     "Invalid MDAnalysis selection"),
])
def test_selection_syntax_error_is_formatted(tmp_path, project, change, fragment):
    path, _ = project("bad_syntax")
    _edit(path, change)
    _user_error(invoke(tmp_path, "validate", path, ok=False), fragment)


def test_out_of_range_replicate_index_is_formatted(tmp_path, project):
    path, root = project("index")
    _user_error(invoke(tmp_path, "replicate", path, "--index", "3", ok=False), "replicate index 3 is out of range 0..0")
    assert not root.exists()


@pytest.mark.parametrize("debug_args,env", [(("--debug",), None), ((), {"MDSCOPE_DEBUG": "1"})])
def test_debug_shows_traceback_for_expected_errors(tmp_path, project, debug_args, env):
    path, _ = project("debug")
    _edit(path, lambda raw: raw["trajectory"].update(stride=0))
    result = invoke(tmp_path, *debug_args, "validate", path, ok=False, env=env)
    assert result.returncode == 1
    assert "Traceback" in result.stderr and "mdscope.errors.ConfigError" in result.stderr
    assert "trajectory.stride must be >= 1" in result.stderr


@pytest.mark.parametrize("stdin", ["", "example\ntiny.pdb\n"])
def test_noninteractive_init_with_exhausted_stdin_is_formatted(tmp_path, stdin):
    _user_error(invoke(tmp_path, "init", "--output", "new.yaml", stdin=stdin, ok=False),
                "init reads one answer per prompt from stdin")
    assert not (tmp_path / "new.yaml").exists()


def test_unexpected_exception_keeps_traceback(tmp_path):
    code = (
        "import sys, mdscope.cli as cli\n"
        "def fault(args): raise ValueError('internal numerical fault')\n"
        "cli.cmd_validate = fault\n"
        "sys.argv = ['mdscope', 'validate', 'unused.yaml']\n"
        "cli.main()\n"
    )
    result = subprocess.run([sys.executable, "-c", code], cwd=tmp_path, capture_output=True,
                            text=True, timeout=60, check=False)
    assert result.returncode == 1
    assert "Traceback" in result.stderr and "ValueError: internal numerical fault" in result.stderr
    assert "mdscope: error" not in result.stderr
