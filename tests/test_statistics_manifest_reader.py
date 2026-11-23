"""Shared lifecycle facts and explicit legacy/strict admission (#216)."""
import json

import pytest

from mdscope.planning import inspect_saved_statistics
from mdscope.reporting.figure_inventory import statistics_status
from mdscope.reporting.report import generate_report
from mdscope.statistics.manifest import read_statistics_manifest


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


def test_manifest_import_does_not_load_estimator_package():
    import subprocess
    import sys
    subprocess.run([sys.executable,'-c',
                    'import sys; import mdscope.statistics.manifest; assert "mdscope.statistics.time_series" not in sys.modules; assert "mdscope.statistics.dispatch" not in sys.modules'],check=True)


def test_disabled_families_hide_stale_catalog_and_findings(tmp_path):
    import pandas as pd

    from mdscope.statistics.coverage import COLUMNS
    root=tmp_path/'project_summary/statistics'
    (root/'time_series').mkdir(parents=True); (root/'coverage').mkdir()
    (root/'manifest.json').write_text(json.dumps({'status':'complete','methods':{},'families':{'time_series':'disabled'},'replicate_paths':[]}))
    pd.DataFrame({'group':[None],'metric':['stale_metric']}).to_parquet(root/'time_series/replicate_statistics.parquet')
    result=inspect_saved_statistics(tmp_path)['families']['time_series']
    assert result['status']=='disabled' and result['identities']==[]
    pd.DataFrame([dict(family='time_series',observable='{"metric":"stale_metric"}',status='above_zero',
                       difference_system2_minus_system1=1.,ci_low=.5,ci_high=1.5,
                       source_table='time_series/system_comparisons.parquet')],columns=COLUMNS).to_parquet(root/'coverage/contrast_coverage.parquet')
    generate_report(root)
    assert pd.read_parquet(root/'findings.parquet').empty


def test_saved_statistics_catalog_covers_every_family(tmp_path):
    """lipid_contacts statistics were written but absent from the catalog."""
    import pandas as pd

    from mdscope.contracts.statistics import STATISTICS_FAMILIES
    root = tmp_path/'project_summary/statistics'
    (root/'lipid_contacts').mkdir(parents=True)
    (root/'manifest.json').write_text(json.dumps({'status': 'complete', 'methods': {},
                                                  'families': {}, 'replicate_paths': []}))
    pd.DataFrame({'system': ['a'], 'replicate': ['rep1'], 'level': ['class'],
                  'leaflet': ['upper'], 'species_a': ['PC'], 'species_b': ['CHOL'],
                  'metric': ['enrichment']}).to_parquet(root/'lipid_contacts/replicate_statistics.parquet')
    families = inspect_saved_statistics(tmp_path)['families']
    assert set(families) == set(STATISTICS_FAMILIES)
    identity, = families['lipid_contacts']['identities']
    assert {key: identity[key] for key in ('level', 'leaflet', 'species_a', 'species_b', 'metric')} == {
        'level': 'class', 'leaflet': 'upper', 'species_a': 'PC',
        'species_b': 'CHOL', 'metric': 'enrichment'}


def test_unfinished_finalization_does_not_advertise_previous_generation_as_current(tmp_path):
    from mdscope.core.persistence import atomic_json

    path = tmp_path / "project_summary/statistics/manifest.json"
    path.parent.mkdir(parents=True)
    previous = '{"status":"complete","methods":{}}'
    path.write_text(previous)
    attempt = tmp_path / ".state/finalization.json"
    atomic_json(attempt, {"status": "running", "config_hash": "new request"})
    assert not read_statistics_manifest(path).complete
    assert statistics_status({"statistics": {}}, tmp_path, "time_series") == "incomplete"
    with pytest.raises(ValueError, match="complete"):
        inspect_saved_statistics(tmp_path)
    assert path.read_text() == previous
    atomic_json(attempt, {"status": "complete", "config_hash": "new request"})
    assert read_statistics_manifest(path).complete
