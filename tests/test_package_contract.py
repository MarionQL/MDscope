"""Stable import and packaging contract."""
from importlib import metadata
from pathlib import Path
import subprocess
import sys

def test_version_and_public_entrypoints():
    import mdscope
    from mdscope.cli import main
    from mdscope.pipeline.runner import run_all
    from mdscope.observables.msd import fit_lateral_diffusion
    assert mdscope.__version__
    assert callable(main) and callable(run_all) and callable(fit_lateral_diffusion)

def test_console_entrypoint_metadata():
    eps = metadata.distribution("mdscope").entry_points
    assert any(ep.group == "console_scripts" and ep.name == "mdscope"
               and ep.value == "mdscope.cli:main" for ep in eps)

def test_module_help():
    result = subprocess.run([sys.executable, "-m", "mdscope.cli", "--help"],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "validate" in result.stdout and "replicate" in result.stdout


"""Validate the packaged module structure, not a working-directory alias."""
from importlib import import_module
from importlib.resources import files
from pathlib import Path

def test_source_layout_and_subpackage_discovery():
    root = Path(__file__).resolve().parents[1]
    assert (root / "src/mdscope/cli.py").is_file()
    assert not (root / "analysis").exists()
    for sub in ("core", "interactions", "mechanisms", "observables",
                "pipeline", "reporting", "statistics"):
        assert import_module(f"mdscope.{sub}").__file__
    assert (files("mdscope") / "cli.py").is_file()

