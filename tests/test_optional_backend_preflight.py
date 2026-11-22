"""Enabled analyses must not pass preflight and then fail for a missing optional backend (#229)."""
import json

import MDAnalysis as mda
import pytest
import yaml

from mdscope.config import load_config
from mdscope.pipeline.preflight import validate_selections
from mdscope.pipeline.runner import run_one
from mdscope_test_support import freesasa_backend, invoke_cli

HINT = "python -m pip install 'mdscope[sasa]'"


def _with_sasa(path, spec):
    raw = yaml.safe_load(path.read_text())
    raw["analyses"]["sasa"] = spec
    path.write_text(yaml.safe_dump(raw))
    return path


@pytest.mark.parametrize("command", [("validate",), ("plan", "--check-topology")])
def test_environment_preflight_reports_missing_freesasa_and_fails(
        tmp_path, project, monkeypatch, command):
    freesasa_backend(monkeypatch, tmp_path, available=False)
    path, root = project("missing_backend")
    _with_sasa(path, {"enabled": True})
    result = invoke_cli(tmp_path, command[0], path, *command[1:], ok=False)
    assert result.returncode == 2 and result.stderr.startswith("mdscope: error: ")
    assert "FreeSASA" in result.stderr and HINT in result.stderr
    assert "Traceback" not in result.stderr
    assert not root.exists()


def test_schema_only_plan_stays_independent_of_installed_backends(tmp_path, project, monkeypatch):
    freesasa_backend(monkeypatch, tmp_path, available=False)
    path, root = project("schema_plan")
    _with_sasa(path, {"enabled": True})
    plan = json.loads(invoke_cli(tmp_path, "plan", path, "--format", "json").stdout)
    assert plan["validation"]["mode"] == "schema_only"
    assert any(item["analysis"] == "sasa" for item in plan["analyses"])
    assert not root.exists()


@pytest.mark.parametrize("spec", [
    {"enabled": False},
    {"enabled": True, "runs": {"paused": {"enabled": False}}},
])
def test_disabled_sasa_does_not_require_freesasa(tmp_path, project, monkeypatch, spec):
    freesasa_backend(monkeypatch, tmp_path, available=False)
    path, _ = project("disabled_sasa")
    _with_sasa(path, spec)
    assert "Valid config" in invoke_cli(tmp_path, "validate", path).stdout


def test_available_backend_passes_environment_preflight(tmp_path, project, monkeypatch):
    freesasa_backend(monkeypatch, tmp_path, available=True)
    path, _ = project("available_backend")
    _with_sasa(path, {"enabled": True})
    assert "Valid config" in invoke_cli(tmp_path, "validate", path).stdout


def test_missing_backend_fails_before_topology_or_outputs(tmp_path, project, monkeypatch):
    freesasa_backend(monkeypatch, tmp_path, available=False)
    path, root = project("before_topology")
    _with_sasa(path, {"enabled": True, "runs": {"environment": {}}})
    cfg = load_config(path)

    def no_universe(*args, **kwargs):
        raise AssertionError("backend preflight opened a topology or trajectory")

    monkeypatch.setattr(mda, "Universe", no_universe)
    with pytest.raises(RuntimeError, match="mdscope\\[sasa\\]"):
        validate_selections(cfg, tmp_path)
    with pytest.raises(RuntimeError, match="mdscope\\[sasa\\]"):
        run_one(path, 0)
    assert not root.exists()


@pytest.mark.parametrize("available", [True, False])
def test_init_enables_sasa_only_with_its_backend(tmp_path, monkeypatch, available):
    freesasa_backend(monkeypatch, tmp_path, available=available)
    output = tmp_path / "initialized.yaml"
    result = invoke_cli(tmp_path, "init", "--output", output,
                        stdin="example\ntiny.pdb\ntiny.dcd\n\n")
    assert yaml.safe_load(output.read_text())["analyses"]["sasa"] == {"enabled": available}
    assert (HINT in result.stdout) is (not available)
