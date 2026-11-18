"""The read-only plan must match production registry contracts without an XTC pass."""

import json
from pathlib import Path

import pytest
import yaml

from mdscope.config import load_config
from mdscope.planning import build_plan
from mdscope.registry import analysis_output_plan, required_artifacts
from mdscope_test_support import invoke_cli


def test_schema_only_plan_ignores_missing_inputs_but_checks_yaml(tmp_path, project):
    path, root = project("plan_schema_only")
    raw = yaml.safe_load(path.read_text())
    raw["systems"]["tiny"]["topology"] = "missing.tpr"
    raw["systems"]["tiny"]["replicates"][0]["trajectory"] = "missing.xtc"
    path.write_text(yaml.safe_dump(raw))
    cfg = load_config(path, check_paths=False)
    plan = build_plan(cfg, path)
    assert plan["validation"]["mode"] == "schema_only"
    assert not plan["validation"]["topology_checked"]
    assert not plan["validation"]["trajectory_frames_checked"]
    assert plan["systems"][0]["replicates"] == 1
    assert not root.exists()
    with pytest.raises(FileNotFoundError):
        load_config(path)
    raw["analyses"]["rmsff"] = True
    path.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match="rmsff"):
        load_config(path, check_paths=False)
    assert not root.exists()


def test_cli_plan_json_does_not_write_or_read_trajectory(tmp_path, project):
    path, root = project("plan_cli")
    result = invoke_cli(tmp_path, "plan", path, "--format", "json")
    plan = json.loads(result.stdout)
    assert plan["validation"]["mode"] == "schema_only"
    assert plan["statistics_families"]
    assert plan["replicates"][0]["id"] == "rep1"
    assert not root.exists()

    inspected = invoke_cli(tmp_path, "plan", path, "--format", "json",
                           "--check-topology")
    topology = json.loads(inspected.stdout)
    assert topology["validation"]["topology_checked"]
    assert topology["topology_counts"]["tiny"]["protein"] > 0
    assert not root.exists()


def test_named_run_registered_products_reuse_actual_registry(tmp_path, project):
    path, root = project("plan_runs")
    cfg = yaml.safe_load(path.read_text())
    cfg["analyses"]["sasa"] = {
        "enabled": True,
        "runs": {
            "isolated": {
                "target_selection": "protein",
                "occluder_selection": "protein",
            },
            "environment": {
                "target_selection": "protein",
                "occluder_selection": "all",
                "reference_run": "isolated",
                "atom_sasa_selection": "protein and name CA",
            },
        },
    }
    cfg["analyses"]["residue_membrane_distance"] = {
        "enabled": True,
        "runs": {
            "phosphate": {"membrane_selection": "name P"},
            "all_membrane": {"membrane_selection": "resname POPC"},
        },
    }
    path.write_text(yaml.safe_dump(cfg))
    resolved = load_config(path)
    plan = build_plan(resolved, path)

    expected = sorted({
        (entry.name, entry.run_id or "default",
         item.relative_to(Path("__replica__")).as_posix())
        for entry in analysis_output_plan(resolved, Path("__replica__"))
        for item in required_artifacts(entry, resolved)
    })
    actual = sorted({
        (item["analysis"], item["run_id"], item["path_within_replica"])
        for item in plan["registered_replica_outputs"]
    })
    assert actual == expected
    assert all(item["columns_status"] == "not_registered"
               for item in plan["registered_replica_outputs"])
    assert {item["run_id"] for item in plan["analyses"]
            if item["analysis"] == "sasa"} == {"isolated", "environment"}
    assert any(item["path_within_replica"].endswith("residue_burial.parquet")
               for item in plan["registered_replica_outputs"])
    assert not root.exists()


def test_family_controls_and_partial_catalog_are_explicit(tmp_path, project):
    path, _ = project("plan_switches")
    cfg = yaml.safe_load(path.read_text())
    cfg["plots"] = {"enabled": True, "families": {
        "profiles": False, "solvation": False,
    }}
    cfg["statistics"] = {"enabled": True, "bootstrap_draws": 100,
                         "families": {"profiles": False}}
    path.write_text(yaml.safe_dump(cfg))
    plan = build_plan(load_config(path), path)
    families = {x["family"]: x for x in plan["statistics_families"]}
    figures = {x["family"]: x for x in plan["figure_families"]}
    assert families["profiles"]["status"] == "disabled"
    assert families["time_series"]["status"] == "enabled_conditional"
    assert figures["profiles"]["status"] == "disabled"
    assert figures["solvation"]["status"] == "disabled"
    assert figures["time_series"]["status"] == "enabled_conditional"
    assert plan["coverage"]["all_generated_files_cataloged"] is False
    assert plan["coverage"]["individual_figure_ids_cataloged"] is False


def test_plan_lists_extra_passes_and_does_not_mutate_config(tmp_path, project):
    from copy import deepcopy

    path, _ = project("plan_passes")
    cfg = load_config(path)
    cfg["analyses"]["membrane"] = {
        "enabled": True, "msd": {"enabled": True, "multi_origin": True}
    }
    cfg["analyses"]["potassium_bridge_audit"]["enabled"] = True
    before = deepcopy(cfg)
    plan = build_plan(cfg, path)
    assert cfg == before
    assert set(plan["additional_trajectory_passes"]) == {
        "multi_origin_msd", "potassium_bridge_audit"
    }
