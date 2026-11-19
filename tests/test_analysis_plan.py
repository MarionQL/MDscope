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
    sasa = [item for item in plan["registered_replica_outputs"]
            if item["analysis"] == "sasa" and item["path_within_replica"].endswith("residue_sasa.parquet")]
    assert len(sasa) == 2
    assert all(item["columns"] == ["time_ps", "topology_residue_index",
                                    "residue", "sasa_nm2"] for item in sasa)
    distance = next(item for item in plan["registered_replica_outputs"]
                    if item["path_within_replica"].endswith("distance_matrix_A.npy"))
    assert distance["shape_axes"] == ["sampled_frames", "selected_residues"]
    assert distance["dtype"] == "float32"
    assert distance["units"] == "angstrom"
    assert plan["coverage"]["partial_verified_column_schemas"] is True
    assert all(item["id"].endswith(":" + item["run_id"])
               for item in plan["analyses"])
    assert any(item["purpose"] != "Not cataloged" for item in plan["analyses"])
    assert any(item["topology_requirements_status"] == "partial_contract"
               for item in plan["analyses"] if item["analysis"] == "sasa")
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
    assert "project_summary/statistics/profiles/system_profiles.parquet" in (
        families["profiles"]["files"])
    assert "project_summary/statistics/time_series/replicate_statistics.parquet" in (
        families["time_series"]["files"])
    assert families["time_series"]["files_status"].startswith("conditional_")


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


def test_plan_summary_output_matches_real_six_frame_run(tmp_path, project):
    """End-to-end check: a promised conditional summary is really produced."""
    from mdscope.pipeline.runner import run_one

    path, root = project("plan_real_output")
    cfg = load_config(path)
    plan = build_plan(cfg, path)
    summary = {entry["path_within_replica"]: entry
               for entry in plan["conditional_replica_outputs"]}
    assert summary["summary/timeseries.parquet"]["status"] == "written_if_nonempty"

    run_one(path, 0)
    observed = root / "tiny" / "rep1" / "summary/timeseries.parquet"
    assert observed.is_file()
    import pandas as pd

    data = pd.read_parquet(observed)
    assert "time_ps" in data and len(data) > 0



def test_plan_provenance_follows_alignment_toggle(tmp_path, project):
    """RMSD and RMSF are source-frame measurements when fitting is disabled."""
    path, _ = project("plan_alignment_source")
    cfg = load_config(path)
    cfg["analyses"]["alignment"]["enabled"] = False
    plan = build_plan(cfg, path)
    descriptors = {item["analysis"]: item for item in plan["analyses"]}
    for metric in ("rmsd", "rmsf"):
        assert descriptors[metric]["coordinate_source"] == "source"
    cfg["analyses"]["alignment"]["enabled"] = True
    fitted = build_plan(cfg, path)
    descriptors = {item["analysis"]: item for item in fitted["analyses"]}
    for metric in ("rmsd", "rmsf"):
        assert descriptors[metric]["coordinate_source"] == "fitted"


def test_implicit_hydration_uses_complete_analysis_descriptor(tmp_path, project):
    """GUI consumers can iterate every selected analysis without missing keys."""
    path, _ = project("plan_implicit_hydration")
    cfg = load_config(path)
    cfg["analyses"]["sdc4_mechanisms"]["enabled"] = True
    cfg["analyses"]["sdc4_mechanisms"]["hs_sugar_water_cutoff_A"] = 3.5
    plan = build_plan(cfg, path)
    hydration = next(item for item in plan["analyses"]
                     if item["analysis"] == "hydration")
    assert hydration["artifact_catalog_status"] == "registered_implicit"
    assert hydration["purpose"]
    assert hydration["required_topology"]
    assert hydration["statistics_eligibility"]
    assert hydration["plot_family_eligibility"]
    assert set(hydration) == set(next(
        item for item in plan["analyses"] if item["analysis"] == "rmsd"))


def test_schema_only_statistics_metric_contracts_are_discoverable(tmp_path, project):
    path, root = project("plan_metric_selection")
    cfg = yaml.safe_load(path.read_text())
    cfg["statistics"]["metric_selection"] = {
        "time_series": [{"group": "PROA", "metric": "rmsd_A"}],
        "profiles": [],
        "interactions": [{"kind": "salt_bridges", "metric": "occupancy"}],
        "solvation_rdfs": [],
        "diffusion": [],
        "coupling": [],
    }
    path.write_text(yaml.safe_dump(cfg))
    plan = build_plan(load_config(path, check_paths=False), path)
    families = {row["family"]: row for row in plan["statistics_families"]}
    assert families["profiles"]["status"] == "disabled"
    assert families["time_series"]["metric_selection"] == [
        {"group": "PROA", "metric": "rmsd_A"}]
    assert families["interactions"]["metric_contract"]["identity_fields"] == [
        "kind", "interaction", "metric"
    ]
    assert "bh_fdr" in families["interactions"]["metric_contract"]["supported_estimators"]
    assert "shell_replicate_contrast" in families["solvation_rdfs"]["metric_contract"]["supported_estimators"]
    assert all(
        row["metric_identity_status"].startswith("saved_data_dependent")
        for family, row in families.items()
        if family in {"time_series", "interactions", "profiles",
                      "solvation_rdfs", "diffusion", "coupling"}
    )
    assert not root.exists()


def test_schema_only_metric_contract_includes_sources_and_inferential_figures(tmp_path, project):
    path, root = project("metric_catalog")
    cfg = yaml.safe_load(path.read_text())
    cfg["statistics"]["metric_selection"] = {
        "time_series": [{"group": "PROA", "metric": "rmsd_A", "inference": False}]
    }
    cfg["plots"]["include"] = ["time_series:rmsd_A:comparison_table"]
    path.write_text(yaml.safe_dump(cfg))
    plan = build_plan(load_config(path, check_paths=False), path)
    by_family = {row["family"]: row for row in plan["statistics_families"]}
    chosen = by_family["time_series"]
    assert chosen["metric_contract"]["inferential_selection"] == (
        "per_entry_inference_boolean_default_true"
    )
    assert chosen["metric_contract"]["configured_statistical_figure_ids"] == [
        "time_series:rmsd_A:comparison_table"
    ]
    assert any(item["analysis"] == "rmsd" and item["run_id"] == "default"
               for item in chosen["candidate_source_runs"])
    assert "replicate_independence" in chosen["metric_contract"]["eligibility_requirements"]
    assert not root.exists()


def test_metric_estimator_flags_validate_through_real_yaml_schema(tmp_path, project):
    path, _ = project("metric_estimator_flags")
    cfg = yaml.safe_load(path.read_text())
    cfg["statistics"]["metric_selection"] = {
        "time_series": [{"group": "PROA", "metric": "rmsd_A",
                         "ci": False, "hypothesis_tests": True}],
        "solvation_rdfs": [{"solvent": "water", "target_class": "protein",
                            "site_class": "acidic", "ci": False}],
    }
    path.write_text(yaml.safe_dump(cfg))
    resolved = load_config(path, check_paths=False)
    plan = build_plan(resolved, path)
    by_family = {row["family"]: row for row in plan["statistics_families"]}
    assert set(by_family["time_series"]["metric_contract"]["independent_estimator_toggles"]) == {
        "ci", "hypothesis_tests"
    }
    assert by_family["solvation_rdfs"]["metric_contract"]["independent_estimator_toggles"] == ["ci"]
    cfg["statistics"]["metric_selection"]["solvation_rdfs"][0]["hypothesis_tests"] = True
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="metric_selection"):
        load_config(path, check_paths=False)


def test_saved_statistics_catalog_reads_real_parquet_schema_and_bounds_identities(tmp_path):
    import pandas as pd
    from mdscope.planning import inspect_saved_statistics
    folder = tmp_path / "project_summary/statistics"
    folder.mkdir(parents=True)
    (folder / "manifest.json").write_text(json.dumps({"status": "complete"}))
    source = folder / "time_series/replicate_statistics.parquet"
    source.parent.mkdir()
    pd.DataFrame({
        "group": ["PROA", "PROA", "PROB", "PROB"],
        "metric": ["rmsd_A", "rmsd_A", "rmsd_A", "radius_of_gyration_A"],
        "mean_ci_low": [1., 1.5, 2., 3.],
        "bootstrap_draws": [20, 20, 20, 20],
    }).to_parquet(source)
    result = inspect_saved_statistics(tmp_path, max_identities=2)
    item = result["families"]["time_series"]
    assert item["status"] == "saved"
    assert len(item["identities"]) == 2
    assert item["identities_truncated"]
    assert item["identities"][0]["id"] == 'statistics:time_series:["PROA","rmsd_A"]'
    assert item["identities"][0]["group"] == "PROA"
    assert len({row["id"] for row in item["identities"]}) == 2
    assert item["identity_columns"] == ["group", "metric"]
    assert "mean_ci_low" in item["files"][source.name]["columns"]
    assert item["files"][source.name]["rows"] == 4
    assert item["files"][source.name]["column_types"]["mean_ci_low"] == "double"
    assert item["files"][source.name]["column_types"]["bootstrap_draws"] == "int64"
    assert item["files"][source.name]["row_grain"] == "per-replica observable"
    assert item["files"][source.name]["row_grain_status"] == "family_identity_source"
    assert result["families"]["profiles"]["status"] == "not_materialized"
    with pytest.raises(ValueError, match="max_identities"):
        inspect_saved_statistics(tmp_path, max_identities=0)
    (folder / "manifest.json").write_text(json.dumps({"status": "running"}))
    with pytest.raises(ValueError, match="incomplete"):
        inspect_saved_statistics(tmp_path)


def test_plan_uses_production_registry_for_configured_plot_ids(tmp_path, project):
    from mdscope.reporting.plot_ids import describe_plot
    path, root = project("plan_figure_ids")
    cfg = yaml.safe_load(path.read_text())
    cfg["plots"]["include"] = [
        "time_series:rmsd_A:timeseries",
        "time_series:rmsd_A:comparison_table",
        "solvation:shell_rdf",
    ]
    cfg["plots"].setdefault("families", {})["solvation"] = False
    path.write_text(yaml.safe_dump(cfg))
    plan = build_plan(load_config(path, check_paths=False), path)
    families = {row["family"]: row for row in plan["figure_families"]}
    by_id = {row["id"]: row for row in families["time_series"]["specific_figures"]}
    assert set(by_id) == {
        "time_series:rmsd_A:timeseries",
        "time_series:rmsd_A:comparison_table",
    }
    for identifier, row in by_id.items():
        assert row["destination"] == describe_plot(identifier)["destination"]
        assert row["source_candidates"] == describe_plot(identifier)["source_candidates"]
        assert row["availability"] == "conditional_on_saved_data"
    assert by_id["time_series:rmsd_A:timeseries"]["requires_statistics"] is False
    assert by_id["time_series:rmsd_A:comparison_table"]["requires_statistics"] is True
    assert families["solvation"]["specific_figures"][0]["availability"] == (
        "disabled_by_configuration"
    )
    assert families["time_series"]["specific_figures_status"].startswith(
        "configured_ids_only"
    )
    assert not root.exists()


def test_plan_lists_only_known_static_candidates_without_an_include_list(tmp_path, project):
    from mdscope.reporting.plot_ids import STATIC_IDS
    path, root = project("plan_static_plot_ids")
    cfg = load_config(path, check_paths=False)
    plan = build_plan(cfg, path)
    families = {row["family"]: row for row in plan["figure_families"]}
    ids = {
        item["id"]
        for family in families.values()
        for item in family["specific_figures"]
    }
    assert ids == STATIC_IDS
    assert families["time_series"]["specific_figures"] == []
    assert "dynamic" in families["time_series"]["specific_figures_status"]
    assert not root.exists()
