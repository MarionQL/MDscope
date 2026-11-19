"""Authoring contracts through the same resolver/planner used by the CLI."""
from copy import deepcopy

import pytest
import yaml

from mdscope.config import load_config, resolve_config
from mdscope.studio.model import change_group, import_document, new_document, review, topology_preflight


def rmsd_document():
    doc = new_document()
    doc["analyses"]["rmsd"]["enabled"] = True
    return doc


def test_new_document_rmsd_without_files_and_opt_in_figures(tmp_path):
    doc = rmsd_document()
    text, plan = review(doc, tmp_path)
    raw = yaml.safe_load(text)
    assert raw["plots"] == {"enabled": False, "include": []}
    assert not raw["statistics"]["enabled"]
    variable, = plan["observable_contracts"]
    assert variable["columns"] == ["time_ps", "rmsd_A"]
    assert variable["units"]["rmsd_A"] == "angstrom"
    assert variable["path_within_replica"] == "summary/timeseries.parquet"
    assert len([fig for fig in variable["figures"] if not fig["requires_statistics"]]) == 2
    choices = {"rmsd": ["time_series:rmsd_A:timeseries"]}
    text, _ = review(doc, tmp_path, choices)
    raw = yaml.safe_load(text)
    assert raw["plots"]["include"] == choices["rmsd"]
    assert raw["plots"]["enabled"]
    doc["analyses"]["rmsd"]["enabled"] = False
    text, _ = review(doc, tmp_path, choices)
    assert yaml.safe_load(text)["plots"]["include"] == []
    doc["analyses"]["rmsd"]["enabled"] = True
    assert yaml.safe_load(review(doc, tmp_path, choices)[0])["plots"]["include"] == choices["rmsd"]
    assert list(tmp_path.iterdir()) == []


def test_groups_preserve_expressions_order_roundtrip_and_supported_keys(tmp_path):
    doc = rmsd_document()
    groups = {"patch": {"label": "Patch", "selection": "  resid 10:60 and backbone  ",
                         "alignment_selection": "segid PROA and name CA"},
              "lipids": {"selection": "resname POPC and name P"}}
    doc["selections"]["groups"] = deepcopy(groups)
    doc["analyses"]["rmsd"]["groups"] = ["patch", "lipids"]
    doc["analyses"]["alignment"]["enabled"] = True
    text, plan = review(doc, tmp_path)
    imported = import_document(text)
    assert imported["selections"]["groups"] == groups
    assert list(imported["selections"]["groups"]) == list(groups)
    assert plan["observable_contracts"][1]["alignment_selection"] == groups["patch"]["alignment_selection"]
    assert plan["observable_contracts"][1]["group_selection"] == groups["patch"]["selection"]
    assert not yaml.safe_load(text)["plots"]["include"]
    change_group(doc, "patch", "patch_copy", duplicate=True)
    assert list(doc["selections"]["groups"]) == ["patch", "patch_copy", "lipids"]
    change_group(doc, "patch", "renamed")
    assert doc["analyses"]["rmsd"]["groups"] == ["renamed", "lipids"]
    change_group(doc, "lipids")
    assert doc["analyses"]["rmsd"]["groups"] == ["renamed"]
    doc["analyses"]["rmsf"]["groups"] = ["renamed"]
    with pytest.raises(ValueError, match="groups"):
        resolve_config(doc)


def test_resolver_identical_to_yaml_entrypoint(tmp_path):
    doc = rmsd_document()
    path = tmp_path / "analysis.yaml"
    path.write_text(yaml.safe_dump(doc))
    assert resolve_config(doc, base=tmp_path) == load_config(path, check_paths=False)
    with pytest.raises(ValueError, match="mapping"):
        import_document("- unexpected list")


def test_topology_preflight_real_counts_no_trajectory_read(project):
    path, root = project("studio_topology")
    doc = yaml.safe_load(path.read_text())
    doc["systems"]["tiny"]["replicates"][0]["trajectory"] = "missing.xtc"
    doc["selections"]["groups"] = {"ca": {"selection": "protein and name CA"}}
    counts = topology_preflight(doc, path.parent)
    details = counts["tiny"]["groups"]["ca"]["measurement"]
    assert details["atoms"] > 0
    assert details["residues"] == details["atoms"]
    assert details["first_atoms"][0]["name"] == "CA"
    assert details["residue_labels"]
    assert not root.exists()
    for selection in ("resname NOTHING", "this is invalid syntax"):
        doc["selections"]["groups"]["ca"]["selection"] = selection
        with pytest.raises(ValueError):
            topology_preflight(doc, path.parent)


def test_app_rmsd_slice_and_statistics(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    at = AppTest.from_file(app.__file__).run()
    assert not at.exception
    at.checkbox(key="enable:rmsd").check().run()
    assert not at.exception
    assert any("ILLUSTRATIVE" in item.value for item in at.caption)
    assert not at.checkbox(key="figure:rmsd:time_series:rmsd_A:timeseries").value
    at.checkbox(key="figure:rmsd:time_series:rmsd_A:timeseries").check().run()
    text = next(item.value for item in at.code if "analyses:" in item.value)
    assert yaml.safe_load(text)["plots"]["include"] == ["time_series:rmsd_A:timeseries"]
    at.checkbox(key="stats:master").check().run()
    at.checkbox(key="stats:rmsd_A:None").check().run()
    assert not at.exception
    text = next(item.value for item in at.code if "analyses:" in item.value)
    raw = yaml.safe_load(text)
    assert raw["statistics"]["metric_selection"]["time_series"] == [
        {"metric": "rmsd_A", "group": None, "ci": True, "hypothesis_tests": True}]
    assert all(value == [] for key, value in raw["statistics"]["metric_selection"].items() if key != "time_series")
    at.checkbox(key="enable:rmsd").uncheck().run()
    assert not at.exception
    assert list(tmp_path.iterdir()) == []


def test_exported_distribution_runs_through_production(project):
    from mdscope.pipeline.runner import run_one, finalize_project
    path, root = project("studio_export")
    doc = rmsd_document()
    fixture = yaml.safe_load(path.read_text())
    doc["systems"] = fixture["systems"]
    doc["project"] = fixture["project"]
    doc["time"] = fixture.get("time", {})
    doc["cache"] = {"aligned_coordinates": False}
    doc["analyses"]["rmsf"]["enabled"] = True
    text, plan = review(doc, path.parent, {"rmsd": ["time_series:rmsd_A:distribution"]})
    path.write_text(text)
    cfg = load_config(path)
    assert cfg["plots"]["distributions"]
    run_one(path, 0)
    import pandas as pd
    for item in plan["structural_profile_outputs"]:
        for replica in plan["replicates"]:
            saved = root / replica["system"] / replica["id"] / item["path_within_replica"]
            assert list(pd.read_parquet(saved).columns) == item["columns"]
    finalize_project(path, phase="plotting")
    assert [file.name for file in root.rglob("*.png")] == ["rmsd_A__distribution.png"]


def test_stale_group_figures_and_unrelated_statistics_are_removed(tmp_path):
    doc = rmsd_document()
    doc["statistics"]["enabled"] = True
    choices = {"rmsd": ["time_series:rmsd_A__group_deleted:timeseries"],
               "statistics": ["time_series:rmsd_A:comparison_table"]}
    raw = yaml.safe_load(review(doc, tmp_path, choices)[0])
    assert raw["plots"]["include"] == []


def test_app_group_editor_migrates_selected_group_ids(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    at = AppTest.from_file(app.__file__).run()
    at.text_input(key="edit:new_group").set_value("patch").run()
    next(button for button in at.button if button.label == "Add group").click().run()
    at.text_input(key="edit:group:patch:selection").set_value("resid 10:60 and backbone").run()
    at.text_input(key="edit:group:patch:alignment_selection").set_value("name CA").run()
    at.checkbox(key="enable:rmsd").check().run()
    at.multiselect(key="edit:groups:rmsd").set_value(["patch"]).run()
    at.text_input(key="edit:group:patch:id").set_value("renamed").run()
    at.button(key="edit:group:patch:Rename").click().run()
    assert not at.exception
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert list(raw["selections"]["groups"]) == ["renamed"]
    assert raw["selections"]["groups"]["renamed"]["selection"] == "resid 10:60 and backbone"
    assert raw["selections"]["groups"]["renamed"]["alignment_selection"] == "name CA"
    assert raw["analyses"]["rmsd"]["groups"] == ["renamed"]
    assert not raw["plots"]["include"]


def test_disabled_imported_figures_do_not_get_enabled_by_authoring(tmp_path):
    doc = rmsd_document()
    doc["plots"] = {"enabled": False, "include": ["solvation:water_and_ions"]}
    raw = yaml.safe_load(review(doc, tmp_path, {"rmsd": []})[0])
    assert not raw["plots"]["enabled"]


def test_app_preserves_legacy_inference_and_plot_switches(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    document = {"systems": new_document()["systems"],
                "plots": {"enabled": True}, "statistics": {"bootstrap_draws": 100, "metric_selection": None}}
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(document)
    at.run()
    assert not at.exception
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert raw["plots"] == document["plots"]
    assert raw["statistics"]["bootstrap_draws"] == 100
    assert raw["statistics"]["enabled"]
    assert raw["statistics"]["metric_selection"] is None
    assert raw["analyses"]["alignment"]["enabled"]
    assert raw["analyses"]["rmsf"]["enabled"]
    at.checkbox(key="stats:exact").check().run()
    assert not at.exception
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert all(value == [] for value in raw["statistics"]["metric_selection"].values())


def test_topology_required_export_invalidates_after_selection_edit(project, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    path, root = project("studio_guard")
    monkeypatch.chdir(path.parent)
    doc = rmsd_document()
    doc["systems"] = yaml.safe_load(path.read_text())["systems"]
    doc["selections"]["groups"] = {"ca": {"selection": "name CA"}}
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = doc
    at.run()
    next(item for item in at.checkbox if item.label == "Require topology preflight before export").check().run()
    assert at.get("download_button")[0].proto.disabled
    next(item for item in at.button if item.label == "Validate topology selections (no XTC)").click().run()
    assert not at.exception
    assert not at.get("download_button")[0].proto.disabled
    at.text_input(key="edit:group:ca:selection").set_value("resname NOTHING").run()
    assert at.get("download_button")[0].proto.disabled
    next(item for item in at.button if item.label == "Validate topology selections (no XTC)").click().run()
    assert any("zero" in item.value for item in at.error)
    assert not at.get("download_button")
    assert not root.exists()


def test_guided_family_fields_reuse_backend_estimator_capabilities():
    from mdscope.studio.model import statistics_selector_fields
    assert statistics_selector_fields("profiles") == ["profile", "inference", "ci", "hypothesis_tests"]
    assert statistics_selector_fields("solvation_rdfs") == ["solvent", "target_class", "site_class", "inference", "ci"]
    assert "hypothesis_tests" not in statistics_selector_fields("diffusion")


def test_dssp_card_contract_matches_writer_fields():
    from mdscope.registry import FusedOutput, artifact_schema
    from pathlib import Path
    schema = artifact_schema(FusedOutput("dssp", Path("assignments.parquet")), Path("assignments.parquet"))
    assert schema["columns"] == ["time_ps", "topology_residue_index", "residue", "dssp"]
    assert schema["row_grain"] == "one residue per analyzed frame per replicate"
    assert schema["units"]["time_ps"] == "ps"


def test_interaction_card_schema_matches_empty_producer():
    import numpy as np
    from pathlib import Path
    from mdscope.registry import FusedOutput, artifact_schema
    from mdscope.interactions.lipid_hydrogen_bonds import tabulate_events
    expected = tabulate_events(np.empty((0, 6)), None, [], [])
    schema = artifact_schema(FusedOutput("lipid_hydrogen_bonds", Path("events.parquet")), Path("events.parquet"))
    assert schema["columns"] == list(expected.columns)


def test_descriptive_interaction_figures_follow_lifetime_control():
    from mdscope.planning import analysis_figure_candidates
    cfg = resolve_config(new_document())
    assert [item["id"] for item in analysis_figure_candidates(cfg, "hydrogen_bonds")] == [
        "interactions:hydrogen_bonds:occupancy_heatmap"]
    cfg["analyses"]["lifetimes"]["enabled"] = True
    candidates = analysis_figure_candidates(cfg, "hydrogen_bonds")
    assert candidates[-1]["id"] == "interactions:hydrogen_bonds:lifetime_distribution"
    assert not any(item["requires_statistics"] for item in candidates)


def test_nonstructural_statistical_figures_survive_export(tmp_path):
    from mdscope.planning import statistical_figure_candidates
    doc = new_document()
    doc["statistics"]["enabled"] = True
    doc["statistics"]["metric_selection"]["profiles"] = [{"profile": "sasa_mean"}]
    doc["statistics"]["metric_selection"]["interactions"] = [{"kind": "hydrogen_bonds", "metric": "occupancy"}]
    figures = statistical_figure_candidates(resolve_config(doc))
    assert all(item["requires_statistics"] for item in figures)
    chosen = "interactions:hydrogen_bonds:occupancy:comparison_table"
    raw = yaml.safe_load(review(doc, tmp_path, {"statistics": [chosen]})[0])
    assert raw["plots"]["include"] == [chosen]


def test_rmsf_profile_contracts_are_not_frame_series(tmp_path):
    from mdscope.registry import structural_profile_contracts
    doc = new_document()
    doc["analyses"]["rmsf"]["enabled"] = True
    doc["selections"]["groups"] = {"patch": {"selection": "name CA"}}
    contracts = structural_profile_contracts(resolve_config(doc))
    assert [item["path_within_replica"] for item in contracts] == ["raw/rmsf.parquet", "raw/group_rmsf.parquet"]
    assert "equivalent_bfactor_A2" in contracts[0]["columns"]
    assert "equivalent_bfactor_A2" not in contracts[1]["columns"]
    assert "time_ps" not in contracts[0]["columns"]
    assert all(item["coordinate_source"] == "source" for item in contracts)
    assert all("cache" in item["statistics_requirement"] for item in contracts)


def test_membrane_and_dssp_scalar_controls_share_exact_contracts(tmp_path):
    doc = new_document()
    doc["analyses"]["membrane"].update(enabled=True, headgroup_selection="name P")
    doc["analyses"]["dssp"]["enabled"] = True
    doc["statistics"]["enabled"] = True
    doc["statistics"]["metric_selection"]["time_series"] = [
        {"metric": "membrane_thickness_A", "group": None}]
    choices = {"membrane": ["time_series:area_per_lipid_A2:timeseries"],
               "statistics": ["time_series:membrane_thickness_A:comparison_table"]}
    text, plan = review(doc, tmp_path, choices)
    metrics = {item["metric"]: item for item in plan["observable_contracts"]}
    assert len(metrics) == 7
    assert metrics["area_per_lipid_A2"]["units"]["area_per_lipid_A2"] == "angstrom^2"
    assert metrics["helix_fraction"]["units"]["helix_fraction"] == "fraction"
    assert yaml.safe_load(text)["plots"]["include"] == [*choices["membrane"], *choices["statistics"]]
    doc["analyses"]["membrane"]["thickness"] = {"enabled": False}
    _, plan = review(doc, tmp_path)
    assert "membrane_thickness_A" not in {item["metric"] for item in plan["observable_contracts"]}


def test_named_sasa_and_dssp_profile_identity_catalog():
    from mdscope.registry import configured_profile_variables
    doc = new_document()
    doc["analyses"]["sasa"] = {"enabled": True, "runs": {
        "isolated": {"target_selection": "protein", "occluder_selection": "protein"},
        "environment": {"target_selection": "protein", "occluder_selection": "all", "reference_run": "isolated"}}}
    doc["analyses"]["dssp"]["enabled"] = True
    items = configured_profile_variables(resolve_config(doc))
    assert [item["profile"] for item in items] == ["sasa_mean__isolated", "sasa_mean__environment",
        "sasa_buried_nm2__environment", "dssp_helix_occupancy", "dssp_sheet_occupancy", "dssp_coil_occupancy"]
    assert all(item["statistics_selector"] == {"profile": item["profile"]} for item in items)


def test_app_guided_rmsf_statistics_are_opt_in(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    at = AppTest.from_file(app.__file__).run()
    at.checkbox(key="enable:rmsf").check().run()
    at.checkbox(key="stats:master").check().run()
    assert not at.checkbox(key="stats:profile:rmsf").value
    at.checkbox(key="stats:profile:rmsf").check().run()
    assert not at.exception
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert raw["statistics"]["metric_selection"]["profiles"] == [{"profile": "rmsf", "ci": True, "hypothesis_tests": True}]
    assert not raw["plots"]["include"]
    assert raw["statistics"]["metric_selection"]["time_series"] == []


def test_disabled_analysis_prunes_export_but_preserves_editable_choices(tmp_path):
    doc = new_document()
    doc["statistics"].update(enabled=True)
    doc["statistics"]["metric_selection"]["time_series"] = [{"metric": "membrane_thickness_A"}]
    doc["statistics"]["metric_selection"]["profiles"] = [{"profile": "rmsf"}]
    doc["plots"] = {"enabled": True, "include": ["time_series:membrane_thickness_A:timeseries", "profiles:rmsf:profile"]}
    original = deepcopy(doc)
    exported = yaml.safe_load(review(doc, tmp_path, {})[0])
    assert exported["statistics"]["metric_selection"]["time_series"] == []
    assert exported["statistics"]["metric_selection"]["profiles"] == []
    assert exported["plots"]["include"] == []
    assert doc == original
    doc["analyses"]["membrane"].update(enabled=True, headgroup_selection="name P")
    doc["analyses"]["rmsf"]["enabled"] = True
    exported = yaml.safe_load(review(doc, tmp_path, {})[0])
    assert exported["statistics"]["metric_selection"]["profiles"] == [{"profile": "rmsf"}]
    assert exported["statistics"]["metric_selection"]["time_series"] == [{"metric": "membrane_thickness_A"}]


def test_group_rename_and_delete_migrate_inference_identity():
    doc = new_document()
    doc["selections"]["groups"] = {"patch": {"selection": "name CA"}}
    selectors = doc["statistics"]["metric_selection"]
    selectors["time_series"] = [{"metric": "rmsd_A", "group": "patch"}]
    selectors["profiles"] = [{"profile": "rmsf__group_patch"}]
    change_group(doc, "patch", "renamed")
    assert selectors["time_series"][0]["group"] == "renamed"
    assert selectors["profiles"][0]["profile"] == "rmsf__group_renamed"
    change_group(doc, "renamed")
    assert selectors["time_series"] == selectors["profiles"] == []


def test_app_named_sasa_runs_update_output_and_statistic_identity(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    at = AppTest.from_file(app.__file__).run()
    at.checkbox(key="enable:sasa").check().run()
    at.checkbox(key="edit:runs:sasa:enabled").check().run()
    at.text_input(key="edit:runs:sasa:new").set_value("isolated").run()
    at.button(key="edit:runs:sasa:add").click().run()
    assert not at.exception
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert raw["analyses"]["sasa"]["runs"] == {"isolated": {}}
    assert any("sasa/isolated/residue_sasa.parquet" in item.value for item in at.code)
    at.checkbox(key="stats:master").check().run()
    assert not at.checkbox(key="stats:profile:sasa_mean__isolated").value
    at.checkbox(key="stats:profile:sasa_mean__isolated").check().run()
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert raw["statistics"]["metric_selection"]["profiles"][0]["profile"] == "sasa_mean__isolated"
    assert not raw["plots"]["include"]


def test_run_deletion_migrates_selectors_and_preserves_referenced_runs():
    from mdscope.studio.model import delete_named_run
    doc = new_document()
    doc["analyses"]["sasa"] = {"enabled": True, "runs": {
        "isolated": {}, "environment": {"reference_run": "isolated"}}}
    selectors = doc["statistics"]["metric_selection"]
    selectors["profiles"] = [{"profile": "sasa_mean__environment"}, {"profile": "sasa_buried_nm2__environment"}]
    doc["plots"] = {"enabled": True, "include": ["profiles:sasa_mean__environment:profile"]}
    before = deepcopy(doc)
    with pytest.raises(ValueError, match="referenced"):
        delete_named_run(doc, "sasa", "isolated")
    assert doc == before
    delete_named_run(doc, "sasa", "environment")
    assert selectors["profiles"] == []
    assert doc["plots"]["include"] == []
    with pytest.raises(ValueError, match="Last run"):
        delete_named_run(doc, "sasa", "isolated")


def test_disabled_run_statistics_are_filtered_without_losing_editable_choice(tmp_path):
    doc = new_document()
    doc["analyses"]["sasa"] = {"enabled": True, "runs": {"active": {}, "paused": {"enabled": False}}}
    doc["statistics"]["metric_selection"]["profiles"] = [{"profile": "sasa_mean__paused"}]
    raw = yaml.safe_load(review(doc, tmp_path, {})[0])
    assert raw["statistics"]["metric_selection"]["profiles"] == []
    assert doc["statistics"]["metric_selection"]["profiles"] == [{"profile": "sasa_mean__paused"}]


def test_density_contract_parity_with_real_accumulator(tmp_path):
    import pandas as pd
    from mdscope.registry import density_output_contracts, configured_profile_variables
    from mdscope.config import analysis_runs
    from mdscope.observables.density_1d import DensityProfileAccumulator, temporal_density_profiles
    from mdscope_test_support import membrane_universe
    doc = new_document()
    doc["analyses"]["density_1d"] = {"enabled": True, "runs": {"lipids": {
        "selections": {"heads": "resname POPC"}}}}
    cfg = resolve_config(doc)
    run, = analysis_runs(cfg, "density_1d")
    u = membrane_universe()
    observer = DensityProfileAccumulator(u, cfg, run.spec, run.run_id)
    for ts in u.trajectory:
        observer.update(u, float(ts.frame))
    observer.flush_frames(tmp_path, 0)
    observer.finalize(tmp_path)
    contract, = density_output_contracts(cfg)
    saved = pd.read_parquet(tmp_path / contract["path_within_replica"])
    assert list(saved.columns) == contract["columns"]
    assert set(saved.density_unit) == {contract["units"]["density"]}
    actual = temporal_density_profiles(tmp_path, cfg)
    assert set(actual) == {item["profile"] for item in configured_profile_variables(cfg)}


def test_distance_run_plan_and_heatmap_lifecycle(tmp_path):
    from mdscope.studio.model import delete_named_run
    doc = new_document()
    doc["analyses"]["residue_membrane_distance"] = {"enabled": True, "runs": {
        "heads": {"membrane_selection": "name P", "contact_cutoff_A": None},
        "all_lipids": {"membrane_selection": "resname POPC"}}}
    _, plan = review(doc, tmp_path)
    profiles = {item["profile"] for item in plan["profile_variables"]}
    assert "residue_membrane_distance__heads__mean_min_distance_A" in profiles
    assert "residue_membrane_distance__heads__fraction_frames_within_cutoff" not in profiles
    identifier = "profiles:residue_membrane_distance__heads:condition1:distance_heatmap"
    raw = yaml.safe_load(review(doc, tmp_path, {"residue_membrane_distance": [identifier]})[0])
    assert raw["plots"]["include"] == [identifier]
    doc["plots"] = {"enabled": True, "include": [identifier]}
    delete_named_run(doc, "residue_membrane_distance", "heads")
    assert doc["plots"]["include"] == []
    raw = yaml.safe_load(review(doc, tmp_path, {"residue_membrane_distance": [identifier]})[0])
    assert raw["plots"]["include"] == []


def test_distance_group_control_is_backend_boolean(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    at = AppTest.from_file(app.__file__).run()
    at.checkbox(key="enable:residue_membrane_distance").check().run()
    assert not at.exception
    at.checkbox(key="edit:setting:residue_membrane_distance:groups").check().run()
    assert not at.exception
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert raw["analyses"]["residue_membrane_distance"]["groups"] is True


def test_group_lifecycle_migrates_density_distance_and_figure_choices():
    doc = new_document()
    doc["selections"]["groups"] = {"patch": {"selection": "name CA"}}
    doc["analyses"]["density_1d"] = {"enabled": True, "selections": {"target": {"selection": "all", "groups": True}}}
    doc["analyses"]["residue_membrane_distance"] = {"enabled": True, "groups": True}
    profiles = ["density_1d__default__target__patch", "residue_membrane_distance__default__group_patch__mean_min_distance_A"]
    doc["statistics"]["metric_selection"]["profiles"] = [{"profile": item} for item in profiles]
    ids = [f"profiles:{item}:profile" for item in profiles]
    doc["plots"] = {"enabled": True, "include": ids}
    choices = {"statistics": list(ids)}
    change_group(doc, "patch", "new", choices=choices)
    expected = [item.replace("patch", "new") for item in profiles]
    assert [item["profile"] for item in doc["statistics"]["metric_selection"]["profiles"]] == expected
    assert doc["plots"]["include"] == choices["statistics"] == [f"profiles:{item}:profile" for item in expected]
    change_group(doc, "new", choices=choices)
    assert doc["plots"]["include"] == choices["statistics"] == []
    assert doc["statistics"]["metric_selection"]["profiles"] == []


def test_imported_inactive_interaction_figures_cannot_leak_through_exact_mode(tmp_path):
    doc = new_document()
    doc["plots"] = {"enabled": True, "include": ["interactions:hydrogen_bonds:occupancy_heatmap"]}
    assert yaml.safe_load(review(doc, tmp_path, {})[0])["plots"]["include"] == []
    assert yaml.safe_load(review(doc, tmp_path)[0])["plots"]["include"] == doc["plots"]["include"]


def test_commit_draft_owns_copy_and_accepts_incomplete_editor_state():
    from mdscope.studio.model import commit_draft
    document = new_document()
    draft = deepcopy(document)
    draft["systems"] = {}
    commit_draft(document, draft)
    draft["project"]["name"] = "later edit"
    assert document["systems"] == {}
    assert document["project"]["name"] != "later edit"


def test_exact_mode_blocks_unowned_imported_ids_without_rewriting_legacy(tmp_path):
    doc = new_document()
    doc["plots"] = {"enabled": True, "include": ["time_series:custom_observable:timeseries"]}
    with pytest.raises(ValueError, match="ownership contracts"):
        review(doc, tmp_path, {})
    assert yaml.safe_load(review(doc, tmp_path)[0])["plots"] == doc["plots"]


def test_shared_catalog_supplies_all_guided_profile_figures():
    from mdscope.planning import authoring_catalog, statistical_figure_candidates
    doc = new_document()
    doc["analyses"]["rmsf"]["enabled"] = True
    doc["statistics"]["metric_selection"]["profiles"] = [{"profile": "rmsf"}]
    cfg = resolve_config(doc)
    assert {item["id"] for item in authoring_catalog(cfg)["figures"]["statistics"]} == {
        item["id"] for item in statistical_figure_candidates(cfg)}


def test_invalid_run_creation_does_not_mutate_draft():
    from mdscope.studio.model import add_named_run
    doc = new_document()
    before = deepcopy(doc)
    with pytest.raises(ValueError, match="reserved"):
        add_named_run(doc, "sasa", "default")
    assert doc == before
