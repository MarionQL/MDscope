"""Saved-result browsing and actual finalization must agree on readiness facts."""
from copy import deepcopy
from pathlib import Path
import json

import numpy as np
import pandas as pd
import pytest
import yaml

from mdscope.config import resolve_config, iter_replicates
from mdscope.registry import analysis_output_plan
from mdscope.studio.model import new_document, inspect_results_directory
from mdscope.pipeline.runner import finalize_project, _replicate_state_prefix


def saved_project(tmp_path):
    doc = new_document()
    doc["analyses"]["sasa"] = {"enabled": True}
    cfg = resolve_config(doc, base=tmp_path)
    root = tmp_path / "saved"
    rep = next(iter(iter_replicates(cfg, tmp_path / "analysis.yaml")))
    repout = root / rep.system / rep.label
    entry = next(iter(analysis_output_plan(cfg, repout)))
    entry.path.parent.mkdir(parents=True)
    pd.DataFrame({"time_ps": [0., 1., 2.], "topology_residue_index": [0, 0, 0],
                  "sasa_nm2": [1., 2., 3.]}).to_parquet(entry.path)
    summary = repout / "summary/timeseries.parquet"
    summary.parent.mkdir()
    pd.DataFrame({"time_ps": [0., 1., 2.]}).to_parquet(summary)
    (root / "resolved_config.yaml").write_text(yaml.safe_dump(cfg))
    marker = root / ".state" / (_replicate_state_prefix(rep) + "-" + entry.checkpoint_key + ".json")
    marker.parent.mkdir()
    marker.write_text(json.dumps({"status": "complete"}))
    return doc, cfg, root, rep, summary


def _row(doc, base, root):
    return next(r for r in inspect_results_directory(doc, base, root) if r["analysis"] == "sasa")


@pytest.mark.parametrize("times", [[0., 0., 2.], [2., 1., 0.], [0., np.nan, 2.], []])
def test_inventory_and_finalize_reject_invalid_summary_without_mutation(tmp_path, times):
    doc, cfg, root, rep, summary = saved_project(tmp_path)
    pd.DataFrame({"time_ps": times}).to_parquet(summary)
    snapshot = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    row = _row(doc, tmp_path, root)
    assert row["status"] == "incompatible"
    assert "time_ps" in row["reason"]
    with pytest.raises(RuntimeError, match="persisted|verified"):
        finalize_project(tmp_path / "analysis.yaml", cfg=cfg, results_dir=root)
    assert all(p.read_bytes() == data for p, data in snapshot.items())
    assert not (root / "project_summary").exists()


@pytest.mark.parametrize("identity,processed", [(["wrong", "rep1"], 3),
                                              (["condition1", "rep1"], 2)])
def test_inventory_detects_checkpoint_identity_and_frame_mismatch(tmp_path, identity, processed):
    doc, cfg, root, rep, summary = saved_project(tmp_path)
    cache = summary.parent.parent / "cache"
    cache.mkdir()
    np.savez(cache / "stream_checkpoint.npz", replicate_identity=identity, processed=processed)
    row = _row(doc, tmp_path, root)
    assert row["status"] == "incompatible"
    assert "checkpoint" in row["reason"]


@pytest.mark.parametrize("section,key,value", [("trajectory", "stride", 2),
                                              ("time", "origin_ps", 100.)])
def test_actual_sampling_settings_invalidate_inventory(tmp_path, section, key, value):
    doc, _, root, _, _ = saved_project(tmp_path)
    changed = deepcopy(doc)
    changed.setdefault(section, {})[key] = value
    row = _row(changed, tmp_path, root)
    assert row["status"] == "incompatible"
    assert "configuration differs" in row["reason"]


def test_reporting_changes_and_saved_only_products_remain_reusable(tmp_path):
    doc, cfg, root, _, _ = saved_project(tmp_path)
    doc["analyses"]["sasa"]["enabled"] = False
    doc["statistics"]["bootstrap_draws"] = 300
    doc["plots"]["dpi"] = 300
    row = _row(doc, tmp_path, root)
    assert row["status"] == "verified_existing"
    assert row["scheduled"] is False


def test_missing_saved_replica_is_reported_instead_of_breaking_inventory(tmp_path):
    doc, cfg, root, _, _ = saved_project(tmp_path)
    doc["systems"]["other"] = deepcopy(doc["systems"]["condition1"])
    rows = inspect_results_directory(doc, tmp_path, root)
    assert any(r["condition"] == "other" and r["status"] == "missing" for r in rows)


def test_saved_only_source_can_feed_offline_statistics_without_changing_yaml(tmp_path, monkeypatch):
    from mdscope.pipeline import runner
    doc, cfg, root, rep, summary = saved_project(tmp_path)
    doc["analyses"]["sasa"]["enabled"] = False
    doc["statistics"] = {"enabled": True, "families": {"profiles": True}}
    current = resolve_config(doc, base=tmp_path)
    before = deepcopy(current)
    seen = {}
    def consume(effective, *args):
        seen["sasa"] = effective["analyses"]["sasa"]["enabled"]
    monkeypatch.setattr(runner, "analyze_project", consume)
    monkeypatch.setattr(runner, "fit_shared_states", lambda *a: None)
    monkeypatch.setattr("mdscope.reporting.report.generate_report", lambda path: path / "Findings.md")
    runner.finalize_project(tmp_path / "analysis.yaml", cfg=current, results_dir=root)
    assert seen["sasa"] is True
    assert current == before
    assert yaml.safe_load((root / "resolved_config.yaml").read_text()) == cfg


def test_verified_sibling_run_survives_unrelated_run_settings_change(tmp_path):
    doc, cfg, root, rep, summary = saved_project(tmp_path)
    # Named runs use their own paths/markers and must be compared independently.
    doc["analyses"]["sasa"]["runs"] = {"A": {}, "B": {"probe_radius_A": 1.8}}
    saved = resolve_config(doc, base=tmp_path)
    repout = summary.parent.parent
    for entry in analysis_output_plan(saved, repout):
        entry.path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"time_ps": [0., 1., 2.], "topology_residue_index": [0, 0, 0],
                      "sasa_nm2": [1., 2., 3.]}).to_parquet(entry.path)
        (root / ".state" / (_replicate_state_prefix(rep) + "-" + entry.checkpoint_key + ".json")).write_text('{"status":"complete"}')
    (root / "resolved_config.yaml").write_text(yaml.safe_dump(saved))
    doc["analyses"]["sasa"]["runs"]["B"]["probe_radius_A"] = 2.0
    by_run = {r["run_id"]: r for r in inspect_results_directory(doc, tmp_path, root)}
    assert by_run["A"]["status"] == "verified_existing"
    assert by_run["B"]["status"] == "incompatible"


def test_coverage_includes_absent_expected_replica_rows():
    from mdscope.studio.model import results_coverage
    coverage = results_coverage([{
        "analysis": "sasa", "run_id": "A", "condition": "control",
        "replica": "rep1", "status": "verified_existing",
        "expected_replicates": ["rep1", "rep2"],
    }])
    assert coverage[0]["status"] == "partial"
    assert coverage[0]["replicas"]["rep2"] == "missing"


def test_saved_catalog_and_exact_export_keep_source_disabled(tmp_path):
    from mdscope.studio.model import saved_source_catalog, review
    doc, cfg, root, rep, summary = saved_project(tmp_path)
    doc["analyses"]["sasa"]["enabled"] = False
    doc["statistics"]["enabled"] = True
    doc["statistics"]["metric_selection"] = {"profiles": [{"profile": "sasa_mean"}]}
    rows = inspect_results_directory(doc, tmp_path, root)
    catalog = saved_source_catalog(doc, tmp_path, rows)
    ids = {f["id"] for f in catalog["figures"]["statistics"]}
    assert "profiles:sasa_mean:profile" in ids
    text, _ = review(doc, tmp_path, {"statistics": ["profiles:sasa_mean:profile"]}, inventory=rows)
    exported = yaml.safe_load(text)
    assert exported["plots"]["include"] == ["profiles:sasa_mean:profile"]
    assert exported["analyses"]["sasa"]["enabled"] is False


def test_default_finalization_failure_does_not_rewrite_saved_configuration(tmp_path):
    doc, cfg, root, rep, summary = saved_project(tmp_path)
    cfg["project"]["output"] = str(root)
    pd.DataFrame({"time_ps": [0., 0., 1.]}).to_parquet(summary)
    saved = (root / "resolved_config.yaml").read_bytes()
    with pytest.raises(RuntimeError, match="time_ps"):
        finalize_project(tmp_path / "analysis.yaml", cfg=cfg)
    assert (root / "resolved_config.yaml").read_bytes() == saved
    assert not (root / "project_summary").exists()
