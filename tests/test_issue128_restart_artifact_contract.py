"""Regression tests for restart marker and persisted summary integrity (#128)."""

import json

import pandas as pd

from mdscope.pipeline.runner import RunState


def test_runstate_accepts_completed_legacy_marker_regardless_of_hash(tmp_path):
    state = RunState(tmp_path)
    state.marker("stage").write_text(json.dumps({
        "config_hash": "expected", "status": "failed"
    }))
    assert not state.complete("stage")
    state.marker("stage").write_text(json.dumps({
        "config_hash": "other", "status": "complete"
    }))
    assert state.complete("stage")
    state.marker("stage").write_text("broken JSON")
    assert not state.complete("stage")
    state.mark("stage", {"status": "failed"})
    assert state.complete("stage")


def test_summary_marker_requires_inventory_and_valid_parquet(tmp_path):
    state = RunState(tmp_path)
    path = tmp_path / "replica" / "summary" / "timeseries.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame({"time_ps": [0.0, 10.0], "value": [1., 2.]}).to_parquet(path)
    state.mark("summary")
    # Historical markers can reuse readable expected products.
    assert state.complete("summary", required_artifacts=(path,))
    state.mark("summary", {"artifacts": {str(path): path.stat().st_size}})
    assert state.complete("summary", required_artifacts=(path,))

    path.unlink()
    assert not state.complete("summary", required_artifacts=())
    pd.DataFrame({"time_ps": [0.0]}).to_parquet(path)
    assert not state.complete("summary", required_artifacts=(path,))
    state.mark("summary", {"artifacts": {str(path): path.stat().st_size}})
    assert state.complete("summary", required_artifacts=(path,))

    # Preserve the file length to check footer integrity independently of size.
    path.write_bytes(b"X" * path.stat().st_size)
    assert not state.complete("summary", required_artifacts=(path,))


def test_summary_marker_rejects_missing_secondary_artifact(tmp_path):
    state = RunState(tmp_path)
    summary = tmp_path / "summary"
    summary.mkdir()
    whole = summary / "timeseries.parquet"
    grouped = summary / "group_timeseries.parquet"
    pd.DataFrame({"time_ps": [0.0]}).to_parquet(whole)
    pd.DataFrame({"group": ["A"], "time_ps": [0.0]}).to_parquet(grouped)
    state.mark("summary", {"artifacts": {
        str(p): p.stat().st_size for p in (whole, grouped)
    }})
    assert state.complete("summary", required_artifacts=(whole, grouped))
    grouped.unlink()
    assert not state.complete("summary", required_artifacts=(whole,))
