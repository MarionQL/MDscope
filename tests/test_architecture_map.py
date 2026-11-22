"""The architecture map is derived from Python syntax, never executed MD code."""
import json
import subprocess
import sys
from pathlib import Path
from runpy import run_path

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
        "def later():\n    from .beta import helper\n    return helper()\n"
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
    assert any(site["source"] == "samplemd.alpha" and
               site["target"] == "samplemd.beta" and
               site["phase"] == "deferred" for site in graph["import_sites"])
    assert any(site["source"] == "samplemd.alpha" and
               site["target"] == "samplemd.beta" and
               site["phase"] == "module_load" for site in graph["import_sites"])


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


def test_solvation_postprocessing_does_not_import_the_observer():
    """The former deferred cycle must not reappear in the static map."""
    graph = run_path(str(SCRIPT))["inspect_package"](ROOT / "src/mdscope")
    edges = {tuple(edge) for edge in graph["imports"]}
    solvation = "mdscope.observables.solvation"
    offline = "mdscope.observables.solvation_postprocess"
    shared = "mdscope.observables.solvation_common"
    assert (solvation, offline) in edges  # Finalize still uses offline derivatives.
    assert (solvation, shared) in edges
    assert (offline, shared) in edges
    assert (offline, solvation) not in edges
    assert not any({solvation, offline} <= set(group) for group in graph["cycles"])


def test_interactions_package_exports_are_explicit_without_changing_named_api():
    """All intended functions remain importable; third-party names do not leak."""
    import mdscope.interactions as package
    from mdscope.interactions.interaction_tables import (
        event_lifetimes,
        interaction_key_columns,
        interaction_tables,
        sample_dt,
    )

    assert package.event_lifetimes is event_lifetimes
    assert package.interaction_key_columns is interaction_key_columns
    assert package.sample_dt is sample_dt
    assert package.interaction_tables is interaction_tables
    assert set(package.__all__) <= set(vars(package))
    assert all(callable(getattr(package, name)) for name in package.__all__)
    assert not any(name in package.__all__ for name in ("np", "pd", "Path", "atomic_parquet"))


def test_pipeline_graphs_reuse_live_architecture_inventory(tmp_path):
    """Diagram builder must work from a checkout, without trajectories."""
    script = ROOT / "scripts/pipeline_graphs.py"
    output = tmp_path / "pipeline-graphs"
    result = subprocess.run(
        [sys.executable, str(script), "--output", str(output)],
        capture_output=True, text=True, check=False, timeout=60,
    )
    assert result.returncode == 0, result.stderr
    inventory = json.loads((output / "architecture.json").read_text())
    assert "mdscope.pipeline.streaming" in inventory["modules"]
    execution = (output / "execution.dot").read_text()
    recovery = (output / "recovery.dot").read_text()
    assert "digraph recovery" in recovery
    assert '"invalid" -> "abort"' in recovery
    assert '"pending" -> "shared"' in recovery
    assert '"checkpoint" -> "replay_stream"' in recovery
    assert "streaming: shared primary frame pass" in execution
    assert '"fusedonly" [' in execution
    assert '"choice" -> "audit"' in execution
    assert '"audit" -> "fusedonly"' in execution
    assert "FusedMembrane" in execution and "FusedExtensions" in execution
    focused = (output / "imports_focused.dot").read_text()
    assert '"mdscope.pipeline.runner"' in focused
    assert '"mdscope.pipeline.fused_all"' in focused
    assert '"mdscope.observables.solvation"' in focused
    assert '"mdscope.observables.solvation_common"' in focused
    full = (output / "imports_full.dot").read_text()
    assert '"mdscope.pipeline.streaming"' in full
    assert '"mdscope.observables.lipid_enrichment"' in full
