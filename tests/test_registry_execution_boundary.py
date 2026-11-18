"""Execution scheduling is separate from the persisted-output contract registry."""

import subprocess
import sys
from pathlib import Path

import pandas as pd

from mdscope.config import load_config
from mdscope.pipeline.replica_tasks import additional_time_series, specialized_tasks
from mdscope.registry import FUSED_OUTPUTS, analysis_output_plan
from mdscope.output_paths import FUSED_OUTPUTS as CANONICAL_OUTPUT_PATHS


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
    assert [task.key for task in tasks] == [
        "residue-membrane-phosphate", "residue-membrane-lipids",
        "membrane-order", "profiles", "hydrogen-bonds", "local_environment",
        "potassium_bridge_audit", "lipid_enrichment", "bilayer_normal",
        "lipid_hydrogen_bonds",
    ]
    assert len({task.key for task in tasks}) == len(tasks)
    assert all(callable(task.run) for task in tasks)
    scalars = {task.key: task.scalar_output for task in tasks
               if task.scalar_output is not None}
    assert scalars == {
        "local_environment": "extensions/local_environment/timeseries.parquet",
        "lipid_enrichment": FUSED_OUTPUTS["lipid_enrichment"],
        "bilayer_normal": FUSED_OUTPUTS["bilayer_normal"],
    }
    # Planning must use the very same output contract as the task scheduler.
    paths = {entry.path.relative_to(Path("__replica__")).as_posix()
             for entry in analysis_output_plan(cfg, Path("__replica__"))}
    assert FUSED_OUTPUTS["lipid_enrichment"] in paths
    assert FUSED_OUTPUTS is CANONICAL_OUTPUT_PATHS


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
    """Producer exports and registry must see identical scheduling decisions."""
    from mdscope.pipeline.analysis_selection import (
        membrane_enabled_names, extensions_enabled_names,
        distance_enabled_names, sdc4_enabled,
    )
    from mdscope.pipeline.fused_membrane import enabled_names as membrane_producer
    from mdscope.pipeline.fused_extensions import enabled_names as extensions_producer
    from mdscope.pipeline.fused_distance import enabled_names as distance_producer
    from mdscope.mechanisms.sdc4_mechanisms import enabled as sdc4_producer

    assert membrane_enabled_names is membrane_producer
    assert extensions_enabled_names is extensions_producer
    assert distance_enabled_names is distance_producer
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
        "             'mdscope.mechanisms.sdc4_mechanisms'):\n"
        "    assert name not in sys.modules, name\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, check=False, timeout=60,
    )
    assert result.returncode == 0, result.stderr
