"""Validate the packaged module structure, not a working-directory alias."""
from importlib import import_module
from importlib.resources import files
from pathlib import Path

def test_source_layout_and_subpackage_discovery():
    root = Path(__file__).resolve().parents[1]
    assert (root / "src/mdscope/analysis/cli.py").is_file()
    assert not (root / "analysis").exists()
    for sub in ("core", "interactions", "mechanisms", "observables",
                "pipeline", "reporting", "statistics"):
        assert import_module(f"mdscope.analysis.{sub}").__file__
    assert (files("mdscope") / "analysis" / "cli.py").is_file()
