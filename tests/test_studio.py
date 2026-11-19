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
    text, _ = review(doc, path.parent, {"rmsd": ["time_series:rmsd_A:distribution"]})
    path.write_text(text)
    cfg = load_config(path)
    assert cfg["plots"]["distributions"]
    run_one(path, 0)
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
                "plots": {"enabled": True}, "statistics": {"bootstrap_draws": 100}}
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(document)
    at.run()
    assert not at.exception
    raw = yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    assert raw["plots"] == document["plots"]
    assert raw["statistics"]["bootstrap_draws"] == 100
    assert raw["statistics"]["enabled"]
    assert "metric_selection" not in raw["statistics"]
    assert raw["analyses"]["alignment"]["enabled"]
    assert raw["analyses"]["rmsf"]["enabled"]


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
