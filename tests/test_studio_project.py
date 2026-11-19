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
    at.text_input(key="edit:system:new").set_value("other").run()
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
    at.button(key="edit:system:other:remove").click().run()
    assert list(exported()["systems"]) == ["condition1"]
    assert not at.exception
    assert list(tmp_path.iterdir()) == []
