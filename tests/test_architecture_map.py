"""The architecture map is derived from Python syntax, never executed MD code."""
from pathlib import Path
from runpy import run_path
import json
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/architecture_map.py"


def test_import_and_call_maps_on_known_small_package(tmp_path):
    package = tmp_path / "samplemd"
    nested = package / "nested"
    nested.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    (package / "alpha.py").write_text(
        "from .beta import helper\n"
        "def entry():\n    return helper()\n"
        "class Engine:\n"
        "    def prepare(self):\n        return self.finish()\n"
        "    def finish(self):\n        return 1\n"
    )
    (package / "beta.py").write_text(
        "from . import alpha\n"
        "def helper():\n    return 42\n"
    )
    (nested / "__init__.py").write_text("")
    (nested / "worker.py").write_text(
        "from ..beta import helper\n"
        "def execute():\n    return helper()\n"
    )
    graph = run_path(str(SCRIPT))["inspect_package"](package)
    edges = {tuple(item) for item in graph["imports"]}
    assert ("samplemd.alpha", "samplemd.beta") in edges
    assert ("samplemd.beta", "samplemd.alpha") in edges
    assert ("samplemd.nested.worker", "samplemd.beta") in edges
    assert ["samplemd.alpha", "samplemd.beta"] in graph["cycles"]
    calls = {tuple(item) for item in graph["calls"]}
    assert ("samplemd.alpha.entry", "samplemd.beta.helper") in calls
    assert ("samplemd.alpha.Engine.prepare", "samplemd.alpha.Engine.finish") in calls
    assert ("samplemd.nested.worker.execute", "samplemd.beta.helper") in calls


def test_command_writes_metadata_and_both_graphs(tmp_path):
    package = tmp_path / "example"
    package.mkdir()
    (package / "__init__.py").write_text("")
    (package / "module.py").write_text("def compute():\n    return 42\n")
    output = tmp_path / "graphs"
    run = subprocess.run(
        [sys.executable, str(SCRIPT), "--package-root", str(package),
         "--output", str(output)],
        capture_output=True, text=True, check=False,
    )
    assert run.returncode == 0, run.stderr
    assert "Modules: 2" in run.stdout
    assert "digraph module_imports" in (output / "modules.dot").read_text()
    assert "digraph statically_resolved_calls" in (output / "calls.dot").read_text()
    payload = json.loads((output / "architecture.json").read_text())
    assert payload["modules"] == ["example", "example.module"]
    assert payload["environment"]["python"]
    assert payload["limitations"]
    assert payload["calls"] == []


def test_current_stats_module_boundaries_are_not_historic_analysis_tree():
    graph = run_path(str(SCRIPT))["inspect_package"](ROOT / "src/mdscope")
    assert "mdscope.pipeline.runner" in graph["modules"]
    assert "mdscope.statistics.common" in graph["modules"]
    assert "mdscope.statistics.profiles" in graph["modules"]
    assert not any(name.startswith("analysis.") for name in graph["modules"])
    imports = {tuple(item) for item in graph["imports"]}
    assert ("mdscope.statistics.time_series", "mdscope.statistics.common") in imports
    assert ("mdscope.statistics.profiles", "mdscope.observables.membrane_order") in imports
    # The removed imports were never referenced by common and should not
    # keep observers loaded merely to analyze autocorrelation.
    assert ("mdscope.statistics.common", "mdscope.observables.membrane_order") not in imports
    assert ("mdscope.statistics.time_series", "mdscope.observables.membrane_order") not in imports
