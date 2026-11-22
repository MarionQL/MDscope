"""Execution scheduling is separate from the persisted-output contract registry."""

import subprocess
import sys
from pathlib import Path

import pandas as pd

from mdscope.config import load_config
from mdscope.output_paths import FUSED_OUTPUTS as CANONICAL_OUTPUT_PATHS
from mdscope.pipeline.replica_tasks import additional_time_series, specialized_tasks
from mdscope.registry import FUSED_OUTPUTS, analysis_output_plan


def test_registry_import_does_not_import_execution_implementations():
    """An output-contract import must not eagerly load trajectory task callables."""
    code = (
        "import sys\n"
        "from mdscope.registry import FUSED_OUTPUTS, analysis_output_plan\n"
        "assert 'mdscope.pipeline.replica_tasks' not in sys.modules\n"
        "assert 'mdscope.observables.lipid_enrichment' not in sys.modules\n"
        "assert FUSED_OUTPUTS['lipid_enrichment'] == "
        "'extensions/lipid_footprint/timeseries.parquet'\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, check=False, timeout=60,
    )
    assert result.returncode == 0, result.stderr


def test_scheduler_preserves_task_order_keys_and_scalar_contracts(project):
    path, _ = project("task_contracts")
    cfg = load_config(path)
    analyses = cfg["analyses"]
    analyses["residue_membrane_distance"].update({
        "enabled": True,
        "runs": {
            "phosphate": {"membrane_selection": "name P"},
            "lipids": {"membrane_selection": "resname POPC"},
        },
    })
    for name in (
        "membrane_order", "sasa", "hydrogen_bonds", "local_environment",
        "potassium_bridge_audit", "lipid_enrichment", "bilayer_normal",
        "lipid_hydrogen_bonds",
    ):
        analyses.setdefault(name, {})["enabled"] = True

    tasks = list(specialized_tasks(cfg))
    # All enabled analyses above use fused observer paths except the K audit.
    assert [task.key for task in tasks] == ["potassium_bridge_audit"]
    assert len({task.key for task in tasks}) == len(tasks)
    assert all(callable(task.run) for task in tasks)
    assert not any(hasattr(task, "scalar_output") for task in tasks)
    # Planning must use the very same output contract as the task scheduler.
    paths = {entry.path.relative_to(Path("__replica__")).as_posix()
             for entry in analysis_output_plan(cfg, Path("__replica__"))}
    assert FUSED_OUTPUTS["lipid_enrichment"] in paths
    assert FUSED_OUTPUTS is CANONICAL_OUTPUT_PATHS


def test_specialized_profiles_only_when_dssp_requires_workers(project):
    """Non-serial DSSP is the only profile scheduling exception."""
    path, _ = project("profile_dispatch")
    cfg = load_config(path)
    for name in ("dssp", "sasa"):
        cfg["analyses"][name]["enabled"] = True

    from mdscope.config.runs import extensions_enabled_names

    assert "dssp" in extensions_enabled_names(cfg)
    assert [task.key for task in specialized_tasks(cfg)] == []

    cfg["analyses"]["dssp"]["backend"] = "parallel"
    assert "dssp" not in extensions_enabled_names(cfg)
    assert "sasa" not in extensions_enabled_names(cfg)
    assert [task.key for task in specialized_tasks(cfg)] == ["profiles"]

    cfg["analyses"]["sasa"]["enabled"] = False
    assert [task.key for task in specialized_tasks(cfg)] == ["profiles"]

    cfg["analyses"]["dssp"]["enabled"] = False
    cfg["analyses"]["sasa"]["enabled"] = True
    assert [task.key for task in specialized_tasks(cfg)] == []
    assert "sasa" in extensions_enabled_names(cfg)


def test_fused_only_analysis_families_have_no_specialized_task(project):
    path, _ = project("fused_only_dispatch")
    cfg = load_config(path)
    for name in (
        "residue_membrane_distance", "membrane_order", "bilayer_normal",
        "hydrogen_bonds", "lipid_hydrogen_bonds", "sasa", "local_environment",
    ):
        cfg["analyses"].setdefault(name, {})["enabled"] = True
    cfg["analyses"]["dssp"]["enabled"] = True
    cfg["analyses"]["dssp"]["backend"] = "serial"
    assert list(specialized_tasks(cfg)) == []


def test_lipid_footprint_directory_is_derived_from_output_registry():
    from mdscope.observables.lipid_enrichment import FOOTPRINT_OUTPUT
    assert FOOTPRINT_OUTPUT == "extensions/lipid_footprint"
    assert (Path(FOOTPRINT_OUTPUT) / "timeseries.parquet").as_posix() == (
        FUSED_OUTPUTS["lipid_enrichment"]
    )


def test_optional_scalar_collection_keeps_missing_files_missing(tmp_path, project):
    path, _ = project("optional_series")
    cfg = load_config(path)
    cfg["analyses"]["lipid_enrichment"]["enabled"] = True
    cfg["analyses"]["sdc4_mechanisms"]["enabled"] = True
    cfg["analyses"].setdefault("lipid_contacts", {})["enabled"] = True
    repout = tmp_path / "replica"
    assert additional_time_series(repout, cfg) == []

    expected = [
        (FUSED_OUTPUTS["lipid_enrichment"], "lipid_footprint_value"),
        ("extensions/sdc4_mechanisms/timeseries.parquet", "sdc4_value"),
        ("extensions/lipid_contacts/timeseries.parquet", "lipid_contacts_value"),
    ]
    for relative, metric in expected:
        output = repout / relative
        output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"time_ps": [0.0, 10.0], metric: [1.0, 2.0]}).to_parquet(
            output, index=False,
        )
    tables = additional_time_series(repout, cfg)
    assert [sorted(set(df.columns) - {"time_ps"}) for df in tables] == [
        [metric] for _, metric in expected
    ]
    assert all(df.time_ps.tolist() == [0.0, 10.0] for df in tables)


def test_output_plan_uses_lightweight_production_predicates(project):
    """Registry and execution share the config.runs scheduling predicates."""
    from mdscope.systems.sdc4.mechanisms import enabled as sdc4_producer
    from mdscope.config.runs import (
        distance_enabled_names,
        extensions_enabled_names,
        membrane_enabled_names,
    )
    from mdscope.config.runs import sdc4_enabled

    assert sdc4_enabled is sdc4_producer

    path, _ = project("lightweight_predicates")
    cfg = load_config(path)
    specs = cfg["analyses"]
    specs["membrane_order"]["enabled"] = True
    specs["sasa"]["enabled"] = True
    specs["dssp"].update({"enabled": True, "backend": "parallel"})
    specs["sdc4_mechanisms"].update({
        "enabled": True, "hs_sugar_water_cutoff_A": 3.5,
    })
    specs["residue_membrane_distance"].update({
        "enabled": True,
        "runs": {"near": {"membrane_selection": "name P"}},
    })
    assert membrane_enabled_names(cfg) == ("membrane_order",)
    assert "sasa" not in extensions_enabled_names(cfg)
    assert "dssp" not in extensions_enabled_names(cfg)
    assert "hydration" in extensions_enabled_names(cfg)
    assert distance_enabled_names(cfg) == ("near",)
    assert sdc4_enabled(cfg)
    entries = analysis_output_plan(cfg, Path("__replica__"))
    assert {entry.name for entry in entries} >= {
        "membrane_order", "sasa", "dssp", "hydration",
        "residue_membrane_distance", "sdc4_mechanisms",
    }


def test_read_only_output_plan_does_not_load_producers():
    code = (
        "import sys\n"
        "from pathlib import Path\n"
        "from mdscope.registry import analysis_output_plan\n"
        "cfg = {'analyses': {name: {'enabled': False} for name in (\n"
        "    'membrane_order', 'sasa', 'dssp', 'residue_membrane_distance',\n"
        "    'sdc4_mechanisms', 'potassium_bridge_audit')}}\n"
        "analysis_output_plan(cfg, Path('__replica__'))\n"
        "for name in ('mdscope.pipeline.fused_membrane',\n"
        "             'mdscope.pipeline.fused_extensions',\n"
        "             'mdscope.pipeline.fused_distance',\n"
        "             'mdscope.systems.sdc4.mechanisms'):\n"
        "    assert name not in sys.modules, name\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, check=False, timeout=60,
    )
    assert result.returncode == 0, result.stderr


def test_scalar_collection_does_not_construct_executable_tasks(
        tmp_path, project, monkeypatch):
    """Fused lipid series remains readable without the obsolete run wrapper."""
    import mdscope.pipeline.replica_tasks as tasks

    path, _ = project("scalar_without_tasks")
    cfg = load_config(path)
    for name in ("local_environment", "lipid_enrichment", "bilayer_normal"):
        cfg["analyses"].setdefault(name, {})["enabled"] = True
    repout = tmp_path / "replica"
    for i, name in enumerate(("local_environment", "lipid_enrichment", "bilayer_normal")):
        output = repout / FUSED_OUTPUTS[name]
        output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"time_ps": [0.0, 10.0],
                      f"value_{i}": [1.0, 2.0]}).to_parquet(output, index=False)

    def never_schedule(*args, **kwargs):
        raise AssertionError("Scalar discovery must not build executable tasks")

    monkeypatch.setattr(tasks, "specialized_tasks", never_schedule)
    tables = tasks.additional_time_series(repout, cfg)
    assert [sorted(set(table) - {"time_ps"}) for table in tables] == [
        ["value_0"], ["value_1"], ["value_2"],
    ]
    assert all(table.time_ps.tolist() == [0.0, 10.0] for table in tables)
    cfg["analyses"]["lipid_enrichment"]["enabled"] = False
    assert [sorted(set(table) - {"time_ps"}) for table in
            tasks.additional_time_series(repout, cfg)] == [
        ["value_0"], ["value_2"],
    ]


def test_lipid_observable_has_no_back_edge_to_fused_executor():
    """The scientific kernel remains independent of the trajectory scheduler."""
    code = (
        "import sys\n"
        "import mdscope.observables.lipid_enrichment as lipid\n"
        "assert not hasattr(lipid, 'run_lipid_enrichment')\n"
        "assert 'mdscope.pipeline.fused_membrane' not in sys.modules\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, check=False, timeout=60,
    )
    assert result.returncode == 0, result.stderr


def test_foundational_imports_do_not_load_configuration_or_pipeline():
    """Planning is lightweight and the geometry core has no config back-edge."""
    code = (
        "import sys\n"
        "from mdscope.core.alignment import _fit_transform\n"
        "assert 'mdscope.config' not in sys.modules\n"
        "from mdscope.registry import analysis_output_plan\n"
        "assert not any(name == 'mdscope.pipeline' or name.startswith('mdscope.pipeline.') "
        "for name in sys.modules), sorted(name for name in sys.modules "
        "if name.startswith('mdscope.pipeline'))\n"
        "from mdscope.contracts.outputs import SAVED_SCALAR_SOURCE_FILES as canonical\n"
        "from mdscope.output_paths import SAVED_SCALAR_SOURCE_FILES as facade\n"
        "assert canonical is facade\n"
    )
    result = subprocess.run([sys.executable, "-c", code],
                            capture_output=True, text=True, check=False, timeout=60)
    assert result.returncode == 0, result.stderr


def test_onboarding_scalar_and_named_profile_use_existing_contracts(project):
    """Two onboarding probes: summary scalar and two configured named profile runs."""
    from mdscope.config.runs import distance_enabled_names
    from mdscope.registry import saved_output_plan

    path, _ = project("onboarding_existing_analysis_probes")
    cfg = load_config(path)
    cfg["analyses"]["radius_of_gyration"]["enabled"] = True
    cfg["analyses"]["residue_membrane_distance"].update({
        "enabled": True,
        "runs": {
            "near": {"enabled": True, "membrane_selection": "name P"},
            "skipped": {"enabled": False, "membrane_selection": "name O"},
        },
    })
    assert distance_enabled_names(cfg) == ("near",)
    output = saved_output_plan(cfg, Path("__replica__"))
    scalar = [entry for entry in output if entry.name == "radius_of_gyration"]
    profile = [entry for entry in output if entry.name == "residue_membrane_distance"]
    assert len(scalar) == 1
    assert scalar[0].path == Path("__replica__/summary/timeseries.parquet")
    assert [(entry.run_id, entry.path.as_posix()) for entry in profile] == [
        ("near", "__replica__/profiles/residue_membrane_distance/near/distance_matrix_A.npy"),
    ]
