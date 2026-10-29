"""Stable import and packaging contract."""
from importlib import metadata
from pathlib import Path
import subprocess
import sys

def test_version_and_public_entrypoints():
    import mdscope
    from mdscope.analysis.cli import main
    from mdscope.analysis.runner import run_all
    from mdscope.analysis.observables.msd import fit_lateral_diffusion
    assert mdscope.__version__
    assert callable(main) and callable(run_all) and callable(fit_lateral_diffusion)

def test_console_entrypoint_metadata():
    eps = metadata.distribution("mdscope").entry_points
    assert any(ep.group == "console_scripts" and ep.name == "mdscope"
               and ep.value == "mdscope.analysis.cli:main" for ep in eps)

def test_module_help():
    result = subprocess.run([sys.executable, "-m", "mdscope.analysis.cli", "--help"],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "validate" in result.stdout and "replicate" in result.stdout
