"""Execution scheduling is separate from the persisted-output contract registry."""

import subprocess
import sys
from pathlib import Path

import pandas as pd

from mdscope.config import load_config
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
