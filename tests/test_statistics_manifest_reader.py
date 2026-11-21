"""Shared lifecycle facts and explicit legacy/strict admission (#216)."""
import json

import pytest

from mdscope.statistics.manifest_reader import read_statistics_manifest
from mdscope.reporting.figure_inventory import statistics_status
from mdscope.reporting.report import generate_report
from mdscope.planning import inspect_saved_statistics


@pytest.mark.parametrize("payload,state", [
    (None, "missing"),
    ("{bad json", "invalid"),
    ([], "invalid"),
    ({"families": []}, "invalid"),
    ({"methods": {}}, "legacy"),
    ({"status": "running", "methods": {}}, "incomplete"),
    ({"status": "complete", "methods": {}}, "complete"),
    ({"status": "complete", "families": {"time_series": "disabled"}}, "complete"),
])
def test_shared_manifest_lifecycle_matrix(tmp_path, payload, state):
    path = tmp_path / "manifest.json"
    if isinstance(payload, str):
        path.write_text(payload)
    elif payload is not None:
        path.write_text(json.dumps(payload))
    facts = read_statistics_manifest(path)
    assert facts.state == state
    assert facts.complete == (state == "complete")


@pytest.mark.parametrize("payload,expected", [
    (None, "legacy_unverified"),
    ({"methods": {}}, "legacy_unverified"),
    ({"status": "running", "methods": {}}, "incomplete"),
    ({"status": "complete", "methods": {}}, "valid"),
    ({"status": "complete"}, "legacy_unverified"),
    ({"status": "complete", "methods": {},
      "families": {"time_series": "disabled"}}, "unavailable"),
    ([], "invalid_manifest"),
])
def test_plotting_legacy_policy_is_explicit(tmp_path, payload, expected):
    path = tmp_path / "project_summary/statistics/manifest.json"
    if payload is not None:
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(payload))
    cfg = {"statistics": {"enabled": True}}
    assert statistics_status(cfg, tmp_path, "time_series") == expected


@pytest.mark.parametrize("payload", [
    {"methods": {}},
    {"status": "running", "methods": {}},
    [],
    {"status": "complete", "families": []},
])
def test_planner_and_findings_require_explicit_complete(tmp_path, payload):
    root = tmp_path / "project_summary/statistics"
    root.mkdir(parents=True)
    (root / "manifest.json").write_text(json.dumps(payload))
    with pytest.raises(ValueError):
        inspect_saved_statistics(tmp_path)
    with pytest.raises(ValueError):
        generate_report(root)


def test_completed_manifest_keeps_legacy_render_read_without_methods(tmp_path):
    path = tmp_path / "project_summary/statistics/manifest.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"status": "complete"}))
    facts = read_statistics_manifest(path)
    assert facts.complete
    assert not facts.methods_present
    assert statistics_status({"statistics": {}}, tmp_path, "time_series") == "legacy_unverified"
