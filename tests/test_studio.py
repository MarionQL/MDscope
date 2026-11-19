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
    assert not at.checkbox(key="figure:tree:time_series:rmsd_A:timeseries").value
    at.checkbox(key="figure:tree:time_series:rmsd_A:timeseries").check().run()
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



def test_exact_figures_override_conflicting_imported_plot_switches(tmp_path):
    from mdscope.config import figure_enabled, plot_enabled

    doc = rmsd_document()
    doc["analyses"]["hydrogen_bonds"]["enabled"] = True
    doc["statistics"]["enabled"] = True
    doc["statistics"]["metric_selection"]["time_series"] = [
        {"metric": "rmsd_A", "group": None}]
    doc["plots"] = {
        "enabled": True, "include": None, "time_series": False,
        "interactions": False, "distributions": False,
        "autocorrelation": False, "bootstrap_tables": False,
        "families": {"time_series": False, "interactions": False, "profiles": False},
    }
    before = deepcopy(doc)
    choices = {
        "rmsd": ["time_series:rmsd_A:distribution"],
        "hydrogen_bonds": ["interactions:hydrogen_bonds:occupancy_heatmap"],
        "statistics": ["time_series:rmsd_A:autocorrelation",
                       "time_series:rmsd_A:bootstrap_table"],
    }
    text, _ = review(doc, tmp_path, choices)
    raw = yaml.safe_load(text)
    plots = raw["plots"]
    selected = [identifier for ids in choices.values() for identifier in ids]
    assert plots["enabled"] is True
    assert plots["include"] == selected
    assert all(plots[key] is True for key in (
        "time_series", "interactions", "distributions",
        "autocorrelation", "bootstrap_tables"))
    assert plots["families"] == {
        "time_series": True, "interactions": True, "profiles": False}
    exported = tmp_path / "analysis.yaml"
    exported.write_text(text)
    cfg = load_config(exported, check_paths=False)
    assert all(figure_enabled(cfg, identifier) for identifier in selected)
    assert not plot_enabled(cfg, "profiles")
    assert doc == before  # Exact export must not rewrite imported legacy mode.
    assert yaml.safe_load(review(doc, tmp_path)[0])["plots"] == before["plots"]


def test_exact_figure_export_reaches_production_renderer(tmp_path):
    import pandas as pd
    from mdscope.reporting.plotting import plot_time_series

    doc = rmsd_document()
    doc["plots"] = {
        "enabled": True, "include": None, "time_series": False,
        "distributions": False, "families": {"time_series": False},
    }
    selected = "time_series:rmsd_A:distribution"
    text, _ = review(doc, tmp_path, {"rmsd": [selected]})
    cfg = resolve_config(yaml.safe_load(text))
    root = tmp_path / "results"
    rep = root / "condition1" / "rep1"
    (rep / "summary").mkdir(parents=True)
    pd.DataFrame({"time_ps": [0., 1., 2., 3., 4.],
                  "rmsd_A": [1., 1.2, 1.5, 1.3, 1.4]}).to_parquet(
                      rep / "summary" / "timeseries.parquet")
    plot_time_series(cfg, root, {"condition1": [rep]})
    assert (root / "project_summary/figures/time_series/rmsd_A__distribution.png").is_file()
    assert not (root / "project_summary/figures/time_series/rmsd_A__timeseries.png").exists()


def test_app_exact_figures_override_imported_disabled_family(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["plots"] = {
        "enabled": True, "include": None, "time_series": False,
        "families": {"time_series": False},
    }
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    assert not at.exception
    exported = lambda: yaml.safe_load(
        next(item.value for item in at.code if "analyses:" in item.value))
    assert exported()["plots"] == doc["plots"]
    at.checkbox(key="figure:tree:time_series:rmsd_A:timeseries").check().run()
    assert not at.exception
    selected = "time_series:rmsd_A:timeseries"
    assert exported()["plots"]["include"] == [selected]
    assert exported()["plots"]["families"]["time_series"] is True
    assert exported()["plots"]["time_series"] is True
    from mdscope.config import figure_enabled
    assert figure_enabled(resolve_config(exported()), selected)


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
    at.checkbox(key="stats:tree:family:time_series").uncheck().run()
    assert not at.exception
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert raw["statistics"]["metric_selection"]["time_series"] == []
    assert "profiles" not in raw["statistics"]["metric_selection"]


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
    assert doc["statistics"]["metric_selection"]["time_series"][0]["group"] == "renamed"
    assert doc["statistics"]["metric_selection"]["profiles"][0]["profile"] == "rmsf__group_renamed"
    change_group(doc, "renamed")
    assert doc["statistics"]["metric_selection"]["time_series"] == doc["statistics"]["metric_selection"]["profiles"] == []


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
    assert doc["statistics"]["metric_selection"]["profiles"] == []
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


def test_inactive_catalog_normalizes_boolean_imports_without_mutation():
    from mdscope.config import DEFAULTS, merge
    from mdscope.planning import authoring_catalog
    doc = new_document()
    doc["analyses"]["sdc4_mechanisms"] = False
    doc["analyses"]["rmsd"] = True
    before = deepcopy(doc)
    catalog = authoring_catalog(merge(DEFAULTS, doc), include_inactive=True)
    assert any(item["metric"] == "rmsd_A" for item in catalog["scalars"])
    assert doc == before


def test_delete_disabled_run_removes_optional_contract_choices():
    from mdscope.studio.model import delete_named_run
    doc = new_document()
    doc["analyses"]["sdc4_mechanisms"] = False
    doc["analyses"]["residue_membrane_distance"] = {"enabled": False, "runs": {
        "old": {"enabled": False}, "keep": {"enabled": False}}}
    profile = "residue_membrane_distance__old__fraction_frames_within_cutoff"
    doc["statistics"]["metric_selection"]["profiles"] = [{"profile": profile}]
    doc["plots"] = {"enabled": True, "include": [f"profiles:{profile}:profile"]}
    choices = {"statistics": [f"profiles:{profile}:profile"], "residue_membrane_distance": ["profiles:residue_membrane_distance__old:condition1:distance_heatmap"]}
    delete_named_run(doc, "residue_membrane_distance", "old", choices=choices)
    assert all(value == [] for value in choices.values())
    assert doc["statistics"]["metric_selection"]["profiles"] == []
    assert doc["plots"]["include"] == []
    assert list(doc["analyses"]["residue_membrane_distance"]["runs"]) == ["keep"]


def test_disabled_membrane_suboption_filters_export_and_preserves_draft(tmp_path):
    doc = new_document()
    doc["analyses"]["membrane"] = {"enabled": True, "headgroup_selection": "name P", "area_per_lipid": {"enabled": False}}
    doc["statistics"]["metric_selection"]["time_series"] = [{"metric": "area_per_lipid_A2", "group": None}]
    before = deepcopy(doc)
    raw = yaml.safe_load(review(doc, tmp_path, {})[0])
    assert raw["statistics"]["metric_selection"]["time_series"] == []
    assert doc == before


def test_return_to_default_run_prunes_choices_and_protects_external_references(tmp_path):
    from mdscope.studio.model import clear_named_runs
    doc = new_document()
    doc["analyses"]["sasa"] = {"enabled": True, "runs": {"custom": {}}}
    profile = "sasa_mean__custom"
    doc["statistics"]["metric_selection"]["profiles"] = [{"profile": profile}]
    choices = {"statistics": [f"profiles:{profile}:profile"]}
    doc["plots"] = {"enabled": True, "include": choices["statistics"][:]}
    doc["analyses"]["sdc4_mechanisms"] = {"enabled": False, "sasa_run": "custom"}
    before = deepcopy(doc)
    with pytest.raises(ValueError, match="referenced"):
        clear_named_runs(doc, "sasa", choices=choices)
    assert doc == before
    doc["analyses"]["sdc4_mechanisms"] = False
    clear_named_runs(doc, "sasa", choices=choices)
    assert "runs" not in doc["analyses"]["sasa"]
    assert choices["statistics"] == []
    assert doc["plots"]["include"] == []
    assert yaml.safe_load(review(doc, tmp_path, choices)[0])["statistics"]["metric_selection"]["profiles"] == []


def test_selections_tab_topology_preview_uses_real_parser(project, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    path, root = project("studio_selection_preview")
    monkeypatch.chdir(path.parent)
    doc = new_document()
    doc["systems"] = yaml.safe_load(path.read_text())["systems"]
    doc["selections"]["groups"] = {"ca": {"selection": "protein and name CA"}}
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = doc
    at.run()
    assert not at.exception
    at.button(key="edit:groups:preflight").click().run()
    assert not at.exception
    assert any("Topology selections parsed" in item.value for item in at.success)
    assert any("ca:" in item.value and "atoms" in item.value for item in at.caption)
    assert not root.exists()

    at.text_input(key="edit:group:ca:selection").set_value("resname NOTHING").run()
    at.button(key="edit:groups:preflight").click().run()
    assert any("zero" in item.value for item in at.error)
    assert not root.exists()


def test_app_global_selection_yaml_deletion_preserves_named_groups(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["selections"]["features"] = "protein"
    doc["selections"]["groups"] = {"patch": {"selection": "resid 10:20"}}
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    assert not at.exception
    editor = at.text_area(key="edit:selections")
    edited = yaml.safe_load(editor.value)
    edited.pop("features")
    editor.set_value(yaml.safe_dump(edited, sort_keys=False)).run()
    assert not at.exception
    exported = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert "features" not in exported["selections"]
    assert exported["selections"]["groups"] == doc["selections"]["groups"]


def test_app_statistics_tree_preserves_dynamic_all_and_exact_choices(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["statistics"]["enabled"] = True
    doc["statistics"]["metric_selection"] = None
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    assert not at.exception
    exported = lambda: yaml.safe_load(
        next(item.value for item in at.code if "analyses:" in item.value))
    assert exported()["statistics"]["metric_selection"] is None
    at.checkbox(key="stats:tree:family:time_series").uncheck().run()
    assert not at.exception
    assert exported()["statistics"]["metric_selection"]["time_series"] == []
    at.checkbox(key="stats:rmsd_A:None").check().run()
    assert not at.exception
    assert exported()["statistics"]["metric_selection"]["time_series"]
    at.checkbox(key="stats:tree:family:time_series").check().run()
    assert not at.exception
    assert "time_series" not in exported()["statistics"]["metric_selection"]


@pytest.mark.parametrize("operation", ["group", "run", "condition"])
def test_lifecycle_failure_preserves_document_and_figure_choices(monkeypatch, operation):
    """A late catalog failure must not leave partially pruned selectors or IDs."""
    from mdscope.studio import model
    from mdscope import planning

    doc = rmsd_document()
    doc["selections"]["groups"] = {"patch": {"selection": "protein"}}
    doc["analyses"]["rmsd"]["groups"] = ["patch"]
    doc["systems"]["condition2"] = deepcopy(doc["systems"]["condition1"])
    doc["analyses"]["sasa"] = {
        "enabled": True,
        "runs": {"first": {}, "second": {}},
    }
    choices = {"rmsd": ["time_series:rmsd_A:timeseries"],
               "sasa": ["profiles:sasa__first:profile"]}
    before, before_choices = deepcopy(doc), deepcopy(choices)
    original = planning.authoring_catalog
    calls = 0

    def fail_after_first(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls >= (2 if operation in {"group", "condition"} else 1):
            raise RuntimeError("injected catalog failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(planning, "authoring_catalog", fail_after_first)
    with pytest.raises(RuntimeError, match="injected catalog failure"):
        if operation == "group":
            model.change_group(doc, "patch", "renamed", choices=choices)
        elif operation == "run":
            model.delete_named_run(doc, "sasa", "first", choices=choices)
        else:
            model.delete_system(doc, "condition2", choices=choices)
    assert doc == before
    assert choices == before_choices


def test_replica_lifecycle_rejects_invalid_removal_without_mutation():
    from mdscope.studio.model import add_replica, remove_replica

    doc = new_document()
    before = deepcopy(doc)
    with pytest.raises(ValueError, match="at least one replica"):
        remove_replica(doc, "condition1", 0)
    assert doc == before
    with pytest.raises(ValueError, match="Unknown replica index"):
        doc["systems"]["condition1"]["replicates"].append(
            {"label": "rep2", "trajectory": "rep2.xtc"})
        before = deepcopy(doc)
        remove_replica(doc, "condition1", 9)
    assert doc == before
    add_replica(doc, "condition1")
    labels = [rep["label"] for rep in doc["systems"]["condition1"]["replicates"]]
    assert labels == ["rep1", "rep2", "rep3"]
    remove_replica(doc, "condition1", 1)
    assert [rep["label"] for rep in doc["systems"]["condition1"]["replicates"]] == ["rep1", "rep3"]


def test_app_rejected_group_rename_preserves_committed_draft_and_choices(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["selections"]["groups"] = {
        "patch": {"selection": "protein"},
        "existing": {"selection": "resid 1:5"},
    }
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.session_state["choices"] = {"rmsd": ["time_series:rmsd_A:timeseries"]}
    at.run()
    assert not at.exception
    before_doc = deepcopy(at.session_state["document"])
    before_choices = deepcopy(at.session_state["choices"])
    at.text_input(key="edit:group:patch:id").set_value("existing").run()
    at.button(key="edit:group:patch:Rename").click().run()
    assert not at.exception
    assert any("already exists" in error.value for error in at.error)
    assert at.session_state["document"] == before_doc
    assert at.session_state["choices"] == before_choices
    assert at.text_input(key="edit:group:patch:id").value == "existing"
    at.text_input(key="edit:group:patch:id").set_value("renamed").run()
    at.button(key="edit:group:patch:Rename").click().run()
    assert not at.exception
    assert "renamed" in at.session_state["document"]["selections"]["groups"]
    assert "patch" not in at.session_state["document"]["selections"]["groups"]


def test_app_rejected_condition_removal_preserves_committed_state(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app, model

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["systems"]["condition2"] = deepcopy(doc["systems"]["condition1"])
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.session_state["choices"] = {"rmsd": ["time_series:rmsd_A:timeseries"]}
    at.run()
    before_doc = deepcopy(at.session_state["document"])
    before_choices = deepcopy(at.session_state["choices"])

    def reject(*args, **kwargs):
        raise ValueError("injected condition failure")

    monkeypatch.setattr(model, "delete_system", reject)
    # The project editor imports the operation into its own module.
    monkeypatch.setattr("mdscope.studio.project_editor.delete_system", reject)
    at.button(key="edit:system:condition2:remove").click().run()
    assert not at.exception
    assert any("injected condition failure" in error.value for error in at.error)
    assert at.session_state["document"] == before_doc
    assert at.session_state["choices"] == before_choices


def test_app_rejected_global_selection_yaml_preserves_export(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["selections"]["groups"] = {"patch": {"selection": "resid 10:20"}}
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    before = deepcopy(at.session_state["document"])
    at.text_area(key="edit:selections").set_value("groups:\n  patch: anything").run()
    assert not at.exception
    assert any("named group editor" in item.value for item in at.error)
    assert at.session_state["document"] == before
    assert at.session_state["document"]["selections"]["groups"] == before["selections"]["groups"]
    at.text_area(key="edit:selections").set_value(
        yaml.safe_dump({key: value for key, value in before["selections"].items() if key != "groups"})
    ).run()
    assert not at.exception
    assert at.session_state["document"]["selections"] == before["selections"]


def test_app_legacy_import_preserves_all_metric_selection_on_noop(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["statistics"]["enabled"] = True
    doc["statistics"]["metric_selection"] = None
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    assert not at.exception
    assert at.checkbox(key="stats:tree:family:time_series").value
    exported = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert exported["statistics"]["metric_selection"] is None


def test_imported_category_switches_remain_backend_compatible():
    from mdscope.config import figure_enabled
    doc = rmsd_document()
    doc["plots"] = {"enabled": True, "include": None, "time_series": True}
    identifier = "time_series:rmsd_A:timeseries"
    assert figure_enabled(resolve_config(doc), identifier)
    doc["plots"]["time_series"] = False
    assert not figure_enabled(resolve_config(doc), identifier)
    doc["plots"]["time_series"] = True
    doc["plots"]["include"] = []
    assert not figure_enabled(resolve_config(doc), identifier)


def test_app_imported_figure_categories_preserved_until_tree_edit(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["plots"] = {"enabled": True, "include": None, "time_series": True,
                    "profiles": False, "interactions": False}
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    assert not at.exception
    exported = lambda: yaml.safe_load(
        next(item.value for item in at.code if "analyses:" in item.value))
    assert exported()["plots"] == doc["plots"]
    assert at.checkbox(key="figure:tree:time_series:rmsd_A:timeseries").value
    at.checkbox(key="figure:tree:time_series:rmsd_A:timeseries").uncheck().run()
    assert not at.exception
    assert exported()["plots"]["include"] == ["time_series:rmsd_A:distribution"]
    assert exported()["plots"]["enabled"] is True


def test_app_hierarchical_figures_support_category_and_individual_override(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["plots"] = {"enabled": False, "include": []}
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    assert not at.exception
    category = at.checkbox(key="figure:tree:family:time_series")
    category.check().run()
    assert not at.exception
    leaf = at.checkbox(key="figure:tree:time_series:rmsd_A:timeseries")
    assert leaf.value
    leaf.uncheck().run()
    assert not at.exception
    assert not at.checkbox(key="figure:tree:family:time_series").value
    at.checkbox(key="figure:tree:all").check().run()
    assert not at.exception
    assert at.checkbox(key="figure:tree:time_series:rmsd_A:timeseries").value


@pytest.mark.parametrize("family", ["time_series", "profiles", "interactions",
                                    "solvation_rdfs", "diffusion", "coupling"])
def test_statistics_family_contract_preserves_unrestricted_dynamic_metrics(family):
    from mdscope.studio.model import (
        statistics_family_mode, set_statistics_family_mode,
        set_statistics_exact_family,
    )
    doc = {"statistics": {"enabled": True, "metric_selection": None}}
    assert statistics_family_mode(doc, family) == "all"
    set_statistics_family_mode(doc, family, "all")
    assert doc["statistics"]["metric_selection"] is None
    before = deepcopy(doc)
    with pytest.raises(ValueError, match="complete dynamic identity"):
        set_statistics_exact_family(doc, family, [])
    assert doc == before
    set_statistics_family_mode(doc, family, "none")
    assert statistics_family_mode(doc, family) == "none"
    assert doc["statistics"]["metric_selection"][family] == []
    set_statistics_family_mode(doc, family, "all")
    assert statistics_family_mode(doc, family) == "all"
    assert family not in doc["statistics"]["metric_selection"]


def test_statistics_exact_conversion_validates_atomically_and_preserves_other_families():
    from mdscope.studio.model import (
        statistics_family_mode, set_statistics_exact_family,
    )
    doc = {"statistics": {"metric_selection": {"profiles": [{"profile": "rmsf"}]}}}
    before = deepcopy(doc)
    with pytest.raises(ValueError, match="complete dynamic identity"):
        set_statistics_exact_family(doc, "time_series", [{"metric": "rmsd_A"}])
    assert doc == before
    with pytest.raises(ValueError, match="require"):
        set_statistics_exact_family(doc, "time_series", [{"unknown": "rmsd_A"}],
                                    catalog_complete=True)
    assert doc == before
    set_statistics_exact_family(doc, "time_series", [{"metric": "rmsd_A"}],
                                catalog_complete=True)
    assert statistics_family_mode(doc, "time_series") == "exact"
    assert doc["statistics"]["metric_selection"]["profiles"] == [{"profile": "rmsf"}]


def test_results_inventory_requires_registered_companions_and_is_read_only(tmp_path):
    from mdscope.studio.model import inspect_results_directory
    doc = rmsd_document()
    doc["analyses"]["sasa"] = {"enabled": True}
    root = tmp_path / "results"
    root.mkdir()
    before = deepcopy(doc)
    inventory = inspect_results_directory(doc, tmp_path, root)
    sasa = [item for item in inventory if item["analysis"] == "sasa"]
    assert sasa
    assert all(item["status"] == "missing" for item in sasa)
    assert doc == before
    artifact = root / "condition1" / "rep1" / "profiles" / "sasa" / "residue_sasa.parquet"
    artifact.parent.mkdir(parents=True)
    artifact.write_bytes(b"placeholder")
    inventory = inspect_results_directory(doc, tmp_path, root)
    sasa = [item for item in inventory if item["analysis"] == "sasa"]
    assert any(item["status"] == "incompatible" for item in sasa)
    assert not (root / "project_summary").exists()


def test_results_inventory_rejects_missing_directory(tmp_path):
    from mdscope.studio.model import inspect_results_directory
    with pytest.raises(ValueError, match="does not exist"):
        inspect_results_directory(rmsd_document(), tmp_path, tmp_path / "missing")


@pytest.mark.parametrize("family", ["time_series", "profiles", "interactions",
                                    "solvation_rdfs", "diffusion", "coupling"])
def test_statistics_imported_family_semantics_are_preserved_without_conversion(family):
    from mdscope.studio.model import statistics_family_mode
    for selection, expected in ((None, "all"), ({}, "all"),
                                ({family: []}, "none")):
        doc = {"statistics": {"metric_selection": deepcopy(selection)}}
        before = deepcopy(doc)
        assert statistics_family_mode(doc, family) == expected
        assert doc == before
    entry = {
        "time_series": {"metric": "rmsd_A"},
        "profiles": {"profile": "rmsf"},
        "interactions": {"kind": "hydrogen_bonds"},
        "solvation_rdfs": {"solvent": "water", "target_class": "protein", "site_class": "backbone"},
        "diffusion": {"metric": "lipid_diffusion"},
        "coupling": {"mechanism": "protein_lipid"},
    }[family]
    doc = {"statistics": {"metric_selection": {family: [entry]}}}
    before = deepcopy(doc)
    assert statistics_family_mode(doc, family) == "exact"
    assert doc == before


def test_app_interaction_occupancy_cutoff_round_trip(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc["statistics"]["enabled"] = True
    doc["statistics"]["interaction_min_occupancy"] = 0.20
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    assert not at.exception
    assert at.slider(key="stats:interaction_min_occupancy").value == 0.20
    at.slider(key="stats:interaction_min_occupancy").set_value(0.35).run()
    assert not at.exception
    exported = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert exported["statistics"]["interaction_min_occupancy"] == pytest.approx(0.35)


def test_statistics_source_status_requires_scheduled_producer(tmp_path):
    from mdscope.studio.model import statistics_source_status
    doc = rmsd_document()
    doc["statistics"]["enabled"] = True
    _, plan = review(doc, tmp_path)
    status = statistics_source_status(plan)
    assert status["time_series"]["status"] == "scheduled"
    assert "rmsd:default" in status["time_series"]["scheduled_sources"]
    assert status["interactions"]["status"] == "unavailable"
    assert status["coupling"]["status"] == "unavailable"
    doc["analyses"]["rmsd"]["enabled"] = False
    _, plan = review(doc, tmp_path)
    assert statistics_source_status(plan)["time_series"]["status"] == "unavailable"


def test_exact_statistics_without_producer_blocks_export_but_unrestricted_does_not(tmp_path):
    from mdscope.studio.model import unavailable_exact_statistics
    doc = rmsd_document()
    doc["statistics"]["enabled"] = True
    doc["statistics"]["metric_selection"] = None
    _, plan = review(doc, tmp_path)
    assert unavailable_exact_statistics(doc, plan) == []
    doc["statistics"]["metric_selection"] = {
        "interactions": [{"kind": "hydrogen_bonds"}]}
    assert unavailable_exact_statistics(doc, plan) == ["interactions"]
    doc["analyses"]["hydrogen_bonds"]["enabled"] = True
    _, plan = review(doc, tmp_path)
    assert unavailable_exact_statistics(doc, plan) == []


def test_results_coverage_does_not_collapse_named_runs_or_replicas():
    from mdscope.studio.model import results_coverage
    rows = [
        {"analysis": "sasa", "run_id": "A", "condition": "control",
         "replica": "rep1", "status": "verified_existing"},
        {"analysis": "sasa", "run_id": "A", "condition": "control",
         "replica": "rep2", "status": "incompatible"},
        {"analysis": "sasa", "run_id": "B", "condition": "control",
         "replica": "rep1", "status": "verified_existing"},
        {"analysis": "sasa", "run_id": "A", "condition": "treated",
         "replica": "rep1", "status": "verified_existing"},
    ]
    coverage = results_coverage(rows)
    assert len(coverage) == 3
    a = next(row for row in coverage
             if row["run_id"] == "A" and row["condition"] == "control")
    assert a["status"] == "partial"
    assert a["verified"] == 1 and a["total"] == 2
    assert a["replicas"] == {"rep1": "verified_existing", "rep2": "incompatible"}


def test_results_inventory_requires_marker_matching_config_and_valid_schema(tmp_path):
    import json
    import pandas as pd
    from mdscope.config import resolve_config, iter_replicates
    from mdscope.pipeline.runner import _replicate_state_prefix
    from mdscope.registry import analysis_output_plan
    from mdscope.studio.model import inspect_results_directory

    doc = rmsd_document()
    doc["analyses"]["sasa"] = {"enabled": True}
    root = tmp_path / "results"
    root.mkdir()
    cfg = resolve_config(doc, base=tmp_path)
    rep = next(iter(iter_replicates(cfg, tmp_path / "analysis.yaml")))
    entry = next(item for item in analysis_output_plan(cfg, root / rep.system / rep.label)
                 if item.name == "sasa")
    entry.path.parent.mkdir(parents=True)
    pd.DataFrame({"time_ps": [0.], "topology_residue_index": [0],
                  "sasa_nm2": [1.]}).to_parquet(entry.path)
    summary = root / rep.system / rep.label / "summary" / "timeseries.parquet"
    summary.parent.mkdir(parents=True)
    pd.DataFrame({"time_ps": [0.], "rmsd_A": [1.]}).to_parquet(summary)
    (root / "resolved_config.yaml").write_text(yaml.safe_dump(cfg))
    def status():
        return next(item for item in inspect_results_directory(doc, tmp_path, root)
                    if item["analysis"] == "sasa")["status"]
    assert status() == "incompatible"
    marker = root / ".state" / (
        _replicate_state_prefix(rep) + "-" + entry.checkpoint_key + ".json")
    marker.parent.mkdir()
    marker.write_text(json.dumps({"status": "complete"}))
    assert status() == "verified_existing"
    unscheduled = deepcopy(doc)
    unscheduled["analyses"]["sasa"]["enabled"] = False
    saved_only = next(item for item in inspect_results_directory(
        unscheduled, tmp_path, root) if item["analysis"] == "sasa")
    assert saved_only["status"] == "verified_existing"
    assert saved_only["scheduled"] is False
    entry.path.write_bytes(b"not a parquet file")
    assert status() == "incompatible"
    pd.DataFrame({"time_ps": [0.], "topology_residue_index": [0],
                  "sasa_nm2": [1.]}).to_parquet(entry.path)
    changed = deepcopy(doc)
    changed["analyses"]["sasa"]["probe_radius_A"] = 1.8
    changed_status = next(item for item in inspect_results_directory(
        changed, tmp_path, root) if item["analysis"] == "sasa")["status"]
    assert changed_status == "incompatible"


def test_statistics_source_status_reports_verified_coverage_without_inventing_inputs(tmp_path):
    from mdscope.studio.model import statistics_source_status
    doc = rmsd_document()
    _, plan = review(doc, tmp_path)
    rows = [
        {"analysis": "sasa", "run_id": "A", "condition": "control",
         "replica": "rep1", "status": "verified_existing"},
        {"analysis": "sasa", "run_id": "A", "condition": "control",
         "replica": "rep2", "status": "missing"},
    ]
    status = statistics_source_status(plan, rows)
    assert status["profiles"]["status"] == "partial_existing"
    coverage = status["profiles"]["existing_coverage"]
    assert len(coverage) == 1
    assert coverage[0]["run_id"] == "A"
    assert coverage[0]["verified"] == 1
    assert coverage[0]["status"] == "partial"


def test_explicit_results_root_rejects_missing_without_initialization(tmp_path):
    from mdscope.pipeline.runner import finalize_project
    doc = rmsd_document()
    with pytest.raises(FileNotFoundError, match="Existing results directory"):
        finalize_project(tmp_path / "analysis.yaml", cfg=resolve_config(doc, base=tmp_path),
                         phase="statistics", results_dir=tmp_path / "absent")
    assert not (tmp_path / "absent").exists()


def test_explicit_results_root_rejects_unverified_before_mutating_saved_provenance(tmp_path):
    from mdscope.pipeline.runner import finalize_project
    doc = rmsd_document()
    root = tmp_path / "saved"
    root.mkdir()
    original = "persisted scientific provenance"
    (root / "resolved_config.yaml").write_text(original)
    with pytest.raises(RuntimeError, match="not verified"):
        finalize_project(tmp_path / "analysis.yaml", cfg=resolve_config(doc, base=tmp_path),
                         phase="statistics", results_dir=root)
    assert (root / "resolved_config.yaml").read_text() == original
    assert not (root / "project_summary").exists()


def test_explicit_results_root_refuses_conflicting_roots(tmp_path):
    from mdscope.pipeline.runner import finalize_project
    with pytest.raises(ValueError, match="either root or results_dir"):
        finalize_project(tmp_path / "analysis.yaml",
                         cfg=resolve_config(rmsd_document(), base=tmp_path),
                         root=tmp_path, results_dir=tmp_path)


def test_explicit_results_root_uses_saved_only_inventory_before_statistics(tmp_path, monkeypatch):
    from mdscope.pipeline import runner
    from mdscope.studio import model
    root = tmp_path / "saved"
    root.mkdir()
    doc = rmsd_document()
    cfg = resolve_config(doc, base=tmp_path)
    seen = {}
    def inventory(document, base, directory):
        seen["root"] = directory
        return [{"analysis": "sasa", "run_id": "A", "condition": "system1",
                 "replica": "rep1", "status": "partial"}]
    monkeypatch.setattr(model, "inspect_results_directory", inventory)
    with pytest.raises(RuntimeError, match="sasa/A/system1/rep1/partial"):
        runner.finalize_project(tmp_path / "analysis.yaml", cfg=cfg,
                                phase="statistics", results_dir=root)
    assert seen["root"] == root
    assert not (root / "resolved_config.yaml").exists()
    assert not (root / "project_summary").exists()


def test_app_global_selection_key_deletion_preserves_named_groups(tmp_path, monkeypatch):
    pytest.importorskip('streamlit')
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc['selections']['protein'] = 'protein'
    doc['selections']['alignment'] = 'backbone'
    doc['selections']['groups'] = {'patch': {'label': 'Patch', 'selection': 'resid 1:5'}}
    at = AppTest.from_file(app.__file__)
    at.session_state['document'] = deepcopy(doc)
    at.run()
    assert not at.exception
    at.text_area(key='edit:selections').set_value('protein: protein\n').run()
    assert not at.exception
    exported = yaml.safe_load(next(item.value for item in at.code if 'analyses:' in item.value))
    assert 'alignment' not in exported['selections']
    assert exported['selections']['groups']['patch']['selection'] == 'resid 1:5'
    assert exported['selections']['groups']['patch']['label'] == 'Patch'


def test_global_selections_replacement_is_atomic_and_exact():
    from mdscope.studio.model import replace_global_selections
    doc = rmsd_document()
    doc['selections'] = {'protein': 'protein', 'alignment': 'backbone',
                         'groups': {'patch': {'selection': 'resid 1:10'}}}
    original = deepcopy(doc)
    for invalid in (None, [], {'groups': {}}):
        with pytest.raises(ValueError):
            replace_global_selections(doc, invalid)
        assert doc == original
    replace_global_selections(doc, {'protein': 'name CA'})
    assert doc['selections'] == {'protein': 'name CA',
        'groups': {'patch': {'selection': 'resid 1:10'}}}
    text, _ = review(doc, '.')
    raw = yaml.safe_load(text)
    assert 'alignment' not in raw['selections']
    assert raw['selections']['groups'] == original['selections']['groups']


def test_exact_figures_override_legacy_false_in_production_resolver(tmp_path):
    from mdscope.config import figure_enabled
    doc = rmsd_document()
    doc['plots'] = {'enabled': True, 'include': None,
                    'time_series': False, 'families': {'time_series': False}}
    figure = 'time_series:rmsd_A:timeseries'
    assert not figure_enabled(resolve_config(doc), figure)
    exported, _ = review(doc, tmp_path, {'rmsd': [figure]})
    cfg = resolve_config(yaml.safe_load(exported))
    assert figure_enabled(cfg, figure)
    assert cfg['plots']['include'] == [figure]
    assert cfg['plots']['time_series'] is True
    assert cfg['plots']['families']['time_series'] is True
    assert doc['plots']['time_series'] is False


def test_rejected_group_and_named_run_mutations_preserve_all_choices():
    from mdscope.studio.model import change_group, add_named_run
    doc = rmsd_document()
    doc['selections']['groups'] = {'patch': {'selection': 'resid 1:10'}}
    doc['statistics']['metric_selection'] = {'time_series': [{'metric': 'rmsd_A', 'group': 'patch'}]}
    choices = {'rmsd': ['time_series:rmsd_A:timeseries']}
    before = deepcopy(doc)
    chosen = deepcopy(choices)
    with pytest.raises(ValueError, match='Unknown group'):
        change_group(doc, 'missing', 'renamed', choices=choices)
    assert doc == before and choices == chosen
    with pytest.raises(ValueError):
        add_named_run(doc, 'rmsd', 'run1')
    assert doc == before and choices == chosen


def test_app_invalid_global_selections_yaml_does_not_commit_partial_draft(tmp_path, monkeypatch):
    pytest.importorskip('streamlit')
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    doc = rmsd_document()
    doc['selections']['groups'] = {'patch': {'selection': 'resid 1:10'}}
    doc['statistics']['metric_selection'] = {'time_series': [{'metric': 'rmsd_A'}]}
    at = AppTest.from_file(app.__file__)
    at.session_state['document'] = deepcopy(doc)
    at.run()
    assert not at.exception
    before = deepcopy(at.session_state['document'])
    at.text_area(key='edit:selections').set_value('groups:\n  patch: invalid\n').run()
    assert not at.exception
    assert at.session_state['document'] == before
    at.text_area(key='edit:selections').set_value('protein: protein\n').run()
    assert not at.exception
    assert at.session_state['document']['selections']['groups'] == before['selections']['groups']
