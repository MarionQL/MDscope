"""Exact per-metric inference selection; legacy YAML remains unchanged."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from mdscope.config import DEFAULTS, statistics_family_config, validate_config
from mdscope.statistics.time_series import analyze_time_series
from mdscope.statistics.selection import validate_metric_selection


def _data(root):
    paths = {"a": [], "b": []}
    for system in paths:
        for i in range(2):
            rep = root / system / f"rep{i}"
            summary = rep / "summary"
            summary.mkdir(parents=True)
            pd.DataFrame({
                "time_ps": np.arange(12, dtype=float) * 5,
                "rmsd_A": np.linspace(1 + i, 3 + i, 12),
                "membrane_thickness_A": np.linspace(30, 34, 12),
            }).to_parquet(summary / "timeseries.parquet")
            pd.DataFrame({
                "time_ps": np.arange(12, dtype=float) * 5,
                "group": ["PROA"] * 12,
                "rmsd_A": np.linspace(1, 2, 12) + i,
            }).to_parquet(summary / "group_timeseries.parquet")
            paths[system].append(rep)
    return paths


def test_selected_group_rmsd_does_not_run_unselected_inference(tmp_path, monkeypatch):
    paths = _data(tmp_path / "inputs")
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["bootstrap_draws"] = 100
    cfg["statistics"]["metric_selection"] = {
        "time_series": [{"group": "PROA", "metric": "rmsd_A"}]
    }
    import mdscope.statistics.time_series as mod
    actual = mod._block_bootstrap
    count = []
    def counted(*args, **kwargs):
        count.append(1)
        return actual(*args, **kwargs)
    monkeypatch.setattr(mod, "_block_bootstrap", counted)
    analyze_time_series(cfg, tmp_path, paths)
    rows = pd.read_parquet(
        tmp_path / "project_summary/statistics/time_series/replicate_statistics.parquet"
    )
    assert set(zip(rows.group, rows.metric)) == {("PROA", "rmsd_A")}
    assert count
    assert len(rows) == 4
    assert "membrane_thickness_A" not in set(rows.metric)


def test_empty_selection_avoids_all_time_series_reads(tmp_path, monkeypatch):
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {"time_series": []}
    assert statistics_family_config(cfg, "time_series") is None
    def forbidden(*args, **kwargs):
        raise AssertionError("Unselected inputs must not be read")
    monkeypatch.setattr(pd, "read_parquet", forbidden)
    analyze_time_series(cfg, tmp_path, {"x": [tmp_path / "x"]})
    assert not (tmp_path / "project_summary/statistics/time_series").exists()


def test_legacy_selector_absent_preserves_all_metrics(tmp_path):
    paths = _data(tmp_path / "inputs")
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["bootstrap_draws"] = 100
    analyze_time_series(cfg, tmp_path, paths)
    rows = pd.read_parquet(
        tmp_path / "project_summary/statistics/time_series/replicate_statistics.parquet"
    )
    assert {"rmsd_A", "membrane_thickness_A"} <= set(rows.metric)
    assert (rows.group == "PROA").any()


@pytest.mark.parametrize("selection", [
    {"time_series": [{"metric": ""}]},
    {"time_series": [{"metric": "rmsd_A", "group": ""}]},
    {"time_series": [{"metric": "rmsd_A", "group": "PROA"},
                     {"metric": "rmsd_A", "group": "PROA"}]},
    {"profiles": [{"metric": "rmsf_A"}]},
    {"time_series": "rmsd_A"},
])
def test_invalid_selectors_fail_validation(selection):
    with pytest.raises(ValueError, match="metric_selection"):
        validate_metric_selection({"metric_selection": selection})
