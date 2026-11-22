"""Saved-result browsing and actual finalization must agree on readiness facts."""
import json
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
import yaml

from mdscope.config import iter_replicates, resolve_config
from mdscope.pipeline.runner import _replicate_state_prefix, finalize_project
from mdscope.registry import analysis_output_plan
from mdscope.studio.model import inspect_results_directory, new_document


def saved_project(tmp_path):
    doc = new_document()
    doc["analyses"]["sasa"] = {"enabled": True}
    cfg = resolve_config(doc, base=tmp_path)
    root = tmp_path / "saved"
    rep = next(iter(iter_replicates(cfg, tmp_path / "analysis.yaml")))
    repout = root / rep.system / rep.label
    entry = next(iter(analysis_output_plan(cfg, repout)))
    entry.path.parent.mkdir(parents=True)
    pd.DataFrame({"time_ps": [0., 1., 2.], "topology_residue_index": [0, 0, 0], "residue": ["A:ALA1"] * 3,
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


def test_complete_products_without_saved_configuration_have_unknown_provenance(tmp_path):
    doc, _, root, _, _ = saved_project(tmp_path)
    (root / "resolved_config.yaml").unlink()
    row = _row(doc, tmp_path, root)
    assert row["status"] == "unknown"
    assert not row["statistics_ready"]
    assert "configuration" in row["reason"]


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
        pd.DataFrame({"time_ps": [0., 1., 2.], "topology_residue_index": [0, 0, 0], "residue": ["A:ALA1"] * 3,
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
    from mdscope.studio.model import review, saved_source_catalog
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


def test_streamed_scalar_reuse_tracks_alignment_provenance(tmp_path, project):
    from mdscope.pipeline.runner import run_one
    path, root = project("saved_scalar")
    run_one(path, 0)
    doc = yaml.safe_load(path.read_text())
    rows = inspect_results_directory(doc, path.parent, root)
    assert next(r for r in rows if r["analysis"] == "rmsd")["status"] == "verified_existing"
    changed = deepcopy(doc)
    changed["analyses"]["alignment"] = {"enabled": False}
    row = next(r for r in inspect_results_directory(changed, path.parent, root)
               if r["analysis"] == "rmsd")
    assert row["status"] == "incompatible"
    assert "alignment" in row["reason"]
    doc["analyses"]["rmsd"] = {"enabled": False}
    assert next(r for r in inspect_results_directory(doc, path.parent, root)
                if r["analysis"] == "rmsd")["status"] == "verified_existing"


@pytest.mark.parametrize("names", [("A", "Z"), ("Z", "A")])
@pytest.mark.parametrize("enabled", [True, False])
def test_saved_only_sibling_retains_its_effective_settings_regardless_of_order(tmp_path, names, enabled):
    from mdscope.config import analysis_runs
    from mdscope.pipeline.persisted_inputs import configuration_with_saved_sources
    saved_only, current_name = names
    doc, cfg, root, rep, summary = saved_project(tmp_path)
    doc["analyses"]["sasa"] = {"enabled": True, "probe_radius_A": 1.4,
                               "runs": {saved_only: {}, current_name: {"probe_radius_A": 2.0}}}
    saved = resolve_config(doc, base=tmp_path)
    for entry in analysis_output_plan(saved, summary.parent.parent):
        entry.path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"time_ps": [0., 1., 2.], "topology_residue_index": [0, 0, 0], "residue": ["A:ALA1"] * 3,
                      "sasa_nm2": [1., 2., 3.]}).to_parquet(entry.path)
        (root / ".state" / (_replicate_state_prefix(rep) + "-" + entry.checkpoint_key + ".json")).write_text('{"status":"complete"}')
    (root / "resolved_config.yaml").write_text(yaml.safe_dump(saved))
    doc["analyses"]["sasa"] = {"enabled": enabled, "probe_radius_A": 2.0,
                               "runs": {current_name: {}}}
    current = resolve_config(doc, base=tmp_path)
    before = deepcopy(current)
    rows = inspect_results_directory(doc, tmp_path, root)
    effective = configuration_with_saved_sources(current, rows)
    runs = {r.run_id: r.spec for r in analysis_runs(effective, "sasa")}
    assert set(runs) == set(names)
    assert runs[saved_only]["probe_radius_A"] == 1.4
    assert runs[current_name]["probe_radius_A"] == 2.0
    assert current == before


def test_streamed_scalar_inventory_rejects_missing_metric_without_mutation(tmp_path, project):
    from mdscope.config import load_config
    from mdscope.pipeline.runner import run_one
    path, root = project("missing_scalar")
    run_one(path, 0)
    cfg = load_config(path)
    rep = next(iter(iter_replicates(cfg, path)))
    summary = root / rep.system / rep.label / "summary/timeseries.parquet"
    table = pd.read_parquet(summary).drop(columns=["rmsd_A"])
    table.to_parquet(summary)
    before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
    row = next(r for r in inspect_results_directory(cfg, path.parent, root) if r["analysis"] == "rmsd")
    assert row["status"] == "incompatible"
    from mdscope.pipeline.persisted_inputs import validate_replica_products
    # Direct finalization must inspect core-stream contracts as well as fused outputs.
    with pytest.raises(RuntimeError, match="rmsd"):
        validate_replica_products(cfg, root, [rep])
    with pytest.raises(RuntimeError, match="rmsd"):
        finalize_project(path, cfg=cfg, results_dir=root)
    assert all(p.read_bytes() == value for p, value in before.items())


def test_saved_rmsf_static_readiness_is_separate_from_temporal_statistics(tmp_path, monkeypatch):
    from mdscope.pipeline.persisted_inputs import unavailable_saved_statistics
    from mdscope.studio.model import (
        review,
        saved_source_catalog,
        unavailable_exact_statistics,
    )
    doc = new_document()
    doc["analyses"]["rmsf"]["enabled"] = True
    cfg = resolve_config(doc, base=tmp_path)
    rep = next(iter(iter_replicates(cfg, tmp_path / "analysis.yaml")))
    root = tmp_path / "saved"
    repout = root / rep.system / rep.label
    (repout / "summary").mkdir(parents=True)
    (repout / "raw").mkdir()
    pd.DataFrame({"time_ps": [0., 1., 2.]}).to_parquet(repout / "summary/timeseries.parquet")
    atoms = pd.DataFrame({"atom_index": [0], "topology_residue_index": [0],
                          "residue": ["A:ALA1"], "atom_name": ["CA"], "rmsf_A": [1.],
                          "equivalent_bfactor_A2": [1.]})
    atoms.to_parquet(repout / "raw/rmsf.parquet")
    (root / "resolved_config.yaml").write_text(yaml.safe_dump(cfg))
    (root / ".state").mkdir()
    (root / ".state" / (_replicate_state_prefix(rep) + "-summary.json")).write_text('{"status":"complete"}')
    doc["analyses"]["rmsf"]["enabled"] = False
    doc["statistics"].update(enabled=True, metric_selection={"profiles": [{"profile": "rmsf"}]})
    rows = inspect_results_directory(doc, tmp_path, root)
    row = next(r for r in rows if r["analysis"] == "rmsf")
    assert row["status"] == "verified_existing"
    assert not row["statistics_ready"]
    assert "cache" in row["statistics_reason"]
    _, plan = review(doc, tmp_path)
    assert unavailable_exact_statistics(doc, plan, rows) == ["profiles"]
    assert "rmsf" not in {p["profile"] for p in saved_source_catalog(doc, tmp_path, rows)["profiles"]}
    before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
    with pytest.raises(RuntimeError, match="temporal RMSF"):
        finalize_project(tmp_path / "analysis.yaml", cfg=resolve_config(doc), results_dir=root)
    assert {p: p.read_bytes() for p in root.rglob("*") if p.is_file()} == before
    default_root = resolve_config(doc)
    default_root["project"]["output"] = str(root)
    with pytest.raises(RuntimeError, match="Temporal RMSF"):
        finalize_project(tmp_path / "analysis.yaml", cfg=default_root)
    assert {p: p.read_bytes() for p in root.rglob("*") if p.is_file()} == before
    disabled = resolve_config(doc)
    disabled["statistics"].setdefault("families", {})["profiles"] = False
    assert unavailable_saved_statistics(disabled, rows, allow_scheduled=False) == []
    cache = repout / "cache"
    cache.mkdir()
    atoms.to_parquet(cache / "feature_atoms.parquet")
    np.savez(cache / "aligned-part-00000.npz", time_ps=[0., 1., 2.], xyz=np.zeros((3, 1, 3)))
    with np.load(cache / "aligned-part-00000.npz") as part:
        archive_type = type(part)
    original_getitem = archive_type.__getitem__
    def bounded_read(archive, key):
        assert key != "xyz", "Browsing must inspect coordinate headers, not load coordinates"
        return original_getitem(archive, key)
    with monkeypatch.context() as bounded:
        bounded.setattr(archive_type, "__getitem__", bounded_read)
        ready = inspect_results_directory(doc, tmp_path, root)
    assert next(r for r in ready if r["analysis"] == "rmsf")["statistics_ready"]
    assert unavailable_exact_statistics(doc, plan, ready) == []
    assert "rmsf" in {p["profile"] for p in saved_source_catalog(doc, tmp_path, ready)["profiles"]}
    np.savez(cache / "aligned-part-00000.npz", time_ps=[0., 1., 9.], xyz=np.zeros((3, 1, 3)))
    stale = inspect_results_directory(doc, tmp_path, root)
    assert next(r for r in stale if r["analysis"] == "rmsf")["status"] == "verified_existing"
    assert unavailable_saved_statistics(resolve_config(doc), stale, allow_scheduled=False) == ["profiles"]


def test_temporal_rmsf_rejects_corruption_beyond_header_buffer(tmp_path):
    from zipfile import ZipFile

    from mdscope.observables.rmsf import temporal_rmsf_readiness
    cache = tmp_path / "cache"
    cache.mkdir()
    (tmp_path / "raw").mkdir()
    atoms = pd.DataFrame({"atom_index": np.arange(100), "residue": ["A:ALA1"] * 100})
    atoms.to_parquet(cache / "feature_atoms.parquet")
    atoms.to_parquet(tmp_path / "raw/rmsf.parquet")
    times = np.arange(1000.)
    path = cache / "aligned-part-00000.npz"
    np.savez(path, time_ps=times, xyz=np.zeros((1000, 100, 3)))
    assert temporal_rmsf_readiness(tmp_path, times)[0]
    with ZipFile(path) as archive:
        info = archive.getinfo("xyz.npy")
        offset = info.header_offset + 30 + len(info.filename) + len(info.extra) + info.compress_size - 1
    with path.open("r+b") as handle:
        handle.seek(offset)
        value = handle.read(1)
        handle.seek(offset)
        handle.write(bytes([value[0] ^ 1]))
    ready, reason = temporal_rmsf_readiness(tmp_path, times)
    assert not ready
    assert "CRC" in reason


@pytest.mark.parametrize("missing_column", ["tail", "carbon_position", "n_lipids"])
def test_registered_companion_schema_rejects_missing_fields_without_mutation(tmp_path, missing_column):
    """Readable companion Parquet files still need their registered columns."""
    from mdscope.pipeline.persisted_inputs import validate_analysis_output
    from mdscope.registry import FusedOutput

    entry = FusedOutput("membrane_order", tmp_path / "timeseries.parquet")
    cfg = {"analyses": {"membrane_order": {"enabled": True}}}
    pd.DataFrame({"time_ps": [0.0, 1.0]}).to_parquet(entry.path)
    mapping = tmp_path / "site_mapping.parquet"
    expected = pd.DataFrame({
        "resname": ["POPC"], "tail": ["sn1"], "carbon_position": [2],
        "atom_name": ["C32"], "n_lipids": [1],
    })
    expected.to_parquet(mapping)
    pd.DataFrame({"profile": pd.Series(dtype=str),
                  "file": pd.Series(dtype=str)}).to_parquet(tmp_path / "profile_index.parquet")

    validate_analysis_output(entry, cfg)
    expected.drop(columns=[missing_column]).to_parquet(mapping)
    before = {p: p.read_bytes() for p in tmp_path.iterdir() if p.is_file()}
    with pytest.raises(RuntimeError, match=f"site_mapping.parquet.*missing columns.*{missing_column}"):
        validate_analysis_output(entry, cfg)
    assert before == {p: p.read_bytes() for p in tmp_path.iterdir() if p.is_file()}


def test_registered_companion_schema_accepts_additional_dynamic_columns(tmp_path):
    """Registry contracts prescribe required columns, not a closed table schema."""
    from mdscope.pipeline.persisted_inputs import validate_analysis_output
    from mdscope.registry import FusedOutput

    entry = FusedOutput("membrane_order", tmp_path / "timeseries.parquet")
    cfg = {"analyses": {"membrane_order": {"enabled": True}}}
    pd.DataFrame({"time_ps": [0.0], "order_mean": [0.4]}).to_parquet(entry.path)
    pd.DataFrame({"resname": ["POPC"], "tail": ["sn1"],
                  "carbon_position": [2], "atom_name": ["C32"], "n_lipids": [1],
                  "new_optional_column": [0]}).to_parquet(tmp_path / "site_mapping.parquet")
    pd.DataFrame({"profile": pd.Series(dtype=str),
                  "file": pd.Series(dtype=str),
                  "other_metadata": pd.Series(dtype=str)}).to_parquet(tmp_path / "profile_index.parquet")
    validate_analysis_output(entry, cfg)


def test_saved_inventory_and_finalizer_share_companion_schema_verdict(tmp_path):
    """A complete marker cannot promote a schema-incomplete companion to ready."""
    doc, _, root, rep, _ = saved_project(tmp_path)
    doc["analyses"]["sasa"]["atom_sasa_selection"] = "name CA"
    cfg = resolve_config(doc, base=tmp_path)
    (root / "resolved_config.yaml").write_text(yaml.safe_dump(cfg))
    atom_file = root / rep.system / rep.label / "profiles/sasa/atom_sasa.parquet"
    atoms = pd.DataFrame({
        "time_ps": [0.0, 1.0, 2.0],
        "atom_index": [0, 0, 0],
        "topology_residue_index": [0, 0, 0],
        "atom_name": ["CA"] * 3,
        "sasa_nm2": [1.0, 2.0, 3.0],
    })
    atoms.to_parquet(atom_file)
    assert _row(doc, tmp_path, root)["status"] == "verified_existing"

    atoms.drop(columns=["atom_name"]).to_parquet(atom_file)
    snapshot = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
    row = _row(doc, tmp_path, root)
    assert row["status"] == "incompatible"
    assert "atom_sasa.parquet" in row["reason"]
    assert "atom_name" in row["reason"]
    with pytest.raises(RuntimeError, match="atom_sasa.parquet.*missing columns"):
        finalize_project(tmp_path / "analysis.yaml", cfg=cfg, results_dir=root)
    assert snapshot == {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
