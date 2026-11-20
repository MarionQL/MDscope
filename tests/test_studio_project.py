"""Guided project setup stays equivalent to CLI configuration input."""
from copy import deepcopy

import pytest
import yaml

from mdscope.studio.model import add_system, delete_system, new_document


def test_condition_mutations_preserve_shared_choices_and_are_failure_atomic():
    doc = new_document()
    before = deepcopy(doc)
    for name in ("", "../escape", "condition1"):
        with pytest.raises(ValueError):
            add_system(doc, name)
        assert doc == before
    with pytest.raises(ValueError, match="at least one"):
        delete_system(doc, "condition1")
    add_system(doc, "second")
    removed = "profiles:residue_membrane_distance__default:second:distance_heatmap"
    keep = "time_series:rmsd_A:timeseries"
    choices = {"residue_membrane_distance": [removed], "rmsd": [keep]}
    doc["plots"] = {"enabled": True, "include": [removed, keep]}
    delete_system(doc, "second", choices=choices)
    assert choices == {"residue_membrane_distance": [], "rmsd": [keep]}
    assert doc["plots"]["include"] == [keep]


def test_guided_setup_roundtrip_and_mode_switch(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app
    monkeypatch.chdir(tmp_path)
    at = AppTest.from_file(app.__file__).run()
    at.text_input(key="edit:project:name").set_value("guided").run()
    at.text_input(key="edit:rep:condition1:0:msd_trajectory").set_value("source.xtc").run()
    at.text_input(key="edit:trajectory:start").set_value("2").run()
    at.text_input(key="edit:trajectory:stop").set_value("20").run()
    at.text_input(key="edit:time:dt_ps").set_value("5").run()
    at.button(key="edit:system:condition1:add_rep").click().run()
    at.button(key="edit:system:add").click().run()
    assert not at.exception
    def exported():
        return yaml.safe_load(next(item.value for item in at.code if "analyses:" in item.value))
    doc = exported()
    assert doc["project"]["name"] == "guided"
    assert len(doc["systems"]["condition1"]["replicates"]) == 2
    assert doc["systems"]["condition1"]["replicates"][0]["msd_trajectory"] == "source.xtc"
    assert doc["trajectory"]["start"] == 2 and doc["trajectory"]["stop"] == 20
    assert doc["time"]["dt_ps"] == 5
    assert not doc["plots"]["enabled"] and not doc["statistics"]["enabled"]
    at.checkbox(key="edit:project:advanced").check().run()
    at.text_area(key="edit:project").set_value("name: advanced\noutput: elsewhere").run()
    at.checkbox(key="edit:project:advanced").uncheck().run()
    assert at.text_input(key="edit:project:name").value == "advanced"
    at.button(key="edit:rep:condition1:0:remove").click().run()
    assert len(exported()["systems"]["condition1"]["replicates"]) == 1
    assert exported()["systems"]["condition1"]["replicates"][0]["label"] == "rep2"
    at.button(key="edit:system:condition2:remove").click().run()
    assert list(exported()["systems"]) == ["condition1"]
    assert not at.exception
    assert list(tmp_path.iterdir()) == []



def test_group_replica_condition_rejections_are_atomic_across_selectors():
    from mdscope.studio.model import (
        add_group, add_replica, add_system, change_group, delete_system,
        remove_replica,
    )

    doc = new_document()
    doc["selections"]["groups"]["patch"] = {"selection": "resid 1:10"}
    doc["statistics"]["metric_selection"]["time_series"] = [
        {"metric": "rmsd_A", "group": "patch", "ci": True}
    ]
    figure = "time_series:rmsd_A:timeseries"
    doc["plots"] = {"enabled": True, "include": [figure]}
    choices = {"rmsd": [figure]}
    baseline = deepcopy(doc)
    initial_choices = deepcopy(choices)

    for operation in (
        lambda: add_group(doc, "patch"),
        lambda: add_group(doc, "../invalid"),
        lambda: change_group(doc, "missing", "new", choices=choices),
        lambda: add_system(doc, "condition1"),
        lambda: add_system(doc, "../invalid"),
        lambda: add_replica(doc, "unknown"),
        lambda: remove_replica(doc, "unknown", 0),
        lambda: remove_replica(doc, "condition1", 0),
        lambda: remove_replica(doc, "condition1", 4),
        lambda: delete_system(doc, "unknown", choices=choices),
        lambda: delete_system(doc, "condition1", choices=choices),
    ):
        with pytest.raises((ValueError, KeyError)):
            operation()
        assert doc == baseline
        assert choices == initial_choices

    add_replica(doc, "condition1")
    assert [rep["label"] for rep in doc["systems"]["condition1"]["replicates"]] == [
        "rep1", "rep2"
    ]
    remove_replica(doc, "condition1", 1)
    assert doc == baseline
    assert choices == initial_choices


def test_guided_editor_rejected_condition_add_keeps_hierarchy_and_recovers(tmp_path, monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    from mdscope.studio import app

    monkeypatch.chdir(tmp_path)
    doc = new_document()
    doc["analyses"]["rmsd"]["enabled"] = True
    doc["statistics"]["metric_selection"]["time_series"] = [
        {"metric": "rmsd_A", "ci": True}
    ]
    at = AppTest.from_file(app.__file__)
    at.session_state["document"] = deepcopy(doc)
    at.run()
    assert not at.exception
    baseline = deepcopy(at.session_state["document"])
    choices = deepcopy(at.session_state["choices"])

    def reject(*args, **kwargs):
        raise ValueError("injected condition failure")

    with monkeypatch.context() as patch:
        patch.setattr("mdscope.studio.project_editor.add_system", reject)
        at.button(key="edit:system:add").click().run()
        assert not at.exception
        assert any("injected condition failure" in item.value for item in at.error)
        assert at.session_state["document"] == baseline
        assert at.session_state["choices"] == choices

    at.button(key="edit:system:add").click().run()
    assert not at.exception
    assert "condition2" in at.session_state["document"]["systems"]
    assert at.session_state["document"]["statistics"]["metric_selection"] == baseline[
        "statistics"
    ]["metric_selection"]



def test_topology_identity_hash_streams_large_selection_without_changing_digest(record_property):
    import hashlib
    import json
    import tracemalloc
    from time import perf_counter
    from types import SimpleNamespace
    from mdscope.core.alignment import selection_details

    class Atoms(list):
        @property
        def residues(self):
            return range(100)

    atoms = Atoms(
        SimpleNamespace(index=i, segid="PROα", resid=i % 100 + 1,
                        resname="ALA", name=f"C{i % 4}")
        for i in range(30000)
    )
    expected = hashlib.sha256(json.dumps([
        (a.segid, a.resid, a.resname, a.name) for a in atoms
    ], ensure_ascii=False).encode()).hexdigest()

    tracemalloc.start()
    started = perf_counter()
    try:
        details = selection_details(atoms)
        elapsed_seconds = perf_counter() - started
        _, peak_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    record_property("selection_details_30k_seconds", elapsed_seconds)
    record_property("selection_details_30k_peak_bytes", peak_bytes)
    assert details["ordered_identity_sha256"] == expected
    assert details["atoms"] == 30000
    assert details["residues"] == 100
    assert details["segments"] == ["PROα"]
    assert len(details["residue_labels"]) == 100
    assert len(details["first_atoms"]) == 5
    assert peak_bytes < 3_000_000, f"Topology preview peak working memory: {peak_bytes:,} bytes"
