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


def test_package_diagram_is_transitively_reduced_and_lists_exceptions(tmp_path):
    """Upward exceptions are listed, mutual imports merge, implied edges drop."""
    package = tmp_path / "mdscope"
    for name in ("base", "mid", "top"):
        (package / name).mkdir(parents=True)
        (package / name / "__init__.py").write_text("")
    (package / "__init__.py").write_text("")
    (package / "base" / "kernel.py").write_text("from mdscope.mid import helper\n")
    (package / "mid" / "helper.py").write_text("from mdscope.base import kernel\n")
    (package / "top" / "run.py").write_text(
        "from mdscope.mid import helper\nfrom mdscope import tool\n"
        "def later():\n    from mdscope.base import kernel\n")
    (package / "tool.py").write_text("from mdscope.base.kernel import x\n")
    namespace = run_path(str(SCRIPT))
    graph = namespace["inspect_package"](package)
    merged = namespace["package_dependency_mermaid"](graph, set())
    assert "base_mid[\"base/ + mid/\"]" in merged
    assert "top --> tool" in merged and "tool --> base_mid" in merged
    assert "top --> base_mid" not in merged and "top -.-> base_mid" not in merged
    excepted = namespace["package_dependency_mermaid"](
        graph, {("mdscope.base.kernel", "mdscope.mid.helper")})
    assert "base_mid" not in excepted and "mid --> base" in excepted
    assert "not drawn: `base.kernel` → `mid.helper`" in excepted


def test_architecture_doc_package_diagram_is_current():
    namespace = run_path(str(SCRIPT))
    graph = namespace["inspect_package"](ROOT / "src/mdscope")
    diagram = namespace["package_dependency_mermaid"](graph, namespace["layering_exceptions"]())
    assert "    contracts -->" not in diagram  # contracts imports nothing from mdscope
    assert not namespace["update_docs"](graph, check=True), (
        "Package diagram in docs/developer/architecture.md is stale; "
        "run: python scripts/architecture_map.py --write-docs")
