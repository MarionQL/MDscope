"""Exact per-metric inference selection; legacy YAML remains unchanged."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from mdscope.config import DEFAULTS, statistics_family_config, validate_config
from mdscope.statistics.time_series import analyze_time_series
from mdscope.statistics_selection import validate_metric_selection


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


def test_unknown_saved_identity_fails_before_statistics_outputs(tmp_path):
    paths = _data(tmp_path / "inputs")
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {
        "time_series": [{"group": "NOT_A_GROUP", "metric": "rmsd_A"}]
    }
    with pytest.raises(ValueError, match="Unknown saved time-series metric"):
        analyze_time_series(cfg, tmp_path, paths)
    assert not (tmp_path / "project_summary/statistics/time_series").exists()


def test_conditional_availability_can_be_selected_without_raw_metric(tmp_path):
    from mdscope.statistics.time_series import collect_time_series
    rep = tmp_path / "rep"
    (rep / "summary").mkdir(parents=True)
    pd.DataFrame({
        "time_ps": np.arange(8, dtype=float),
        "sdc4_local_upper_height_A": [1, 2, np.nan, 4, 5, 6, 7, 8],
    }).to_parquet(rep / "summary/timeseries.parquet")
    key = (None, "sdc4_local_upper_height_A_available_fraction")
    series = collect_time_series(rep, selected={key})
    assert set(series) == {key}
    np.testing.assert_array_equal(series[key]["value"].to_numpy(), [1, 1, 0, 1, 1, 1, 1, 1])


def test_interaction_identity_metric_selector_skips_other_inference(tmp_path, monkeypatch):
    from test_interaction_statistics_reference import fixture
    from mdscope.statistics.interactions import analyze_interactions
    cfg, paths = fixture(tmp_path / "input")
    cfg["statistics"]["metric_selection"] = {"interactions": [
        {"kind": "salt_bridges", "interaction": "ASP0 | LYS7 | interchain",
         "metric": "occupancy"}
    ]}
    import mdscope.statistics.interactions as module
    original = module._interaction_bootstrap
    called = []
    def counted(*args, **kwargs):
        called.append(1)
        return original(*args, **kwargs)
    monkeypatch.setattr(module, "_interaction_bootstrap", counted)
    analyze_interactions(cfg, tmp_path, paths)
    dest = tmp_path / "project_summary/statistics/interactions"
    replica = pd.read_parquet(dest / "replicate_interaction_metrics.parquet")
    assert set(replica.metric) == {"occupancy"}
    assert set(replica.interaction) == {"ASP0 | LYS7 | interchain"}
    assert len(called) == 6  # Exactly one chosen identity per replica
    assert len(replica) == 6
    assert replica.loc[(replica.system == "a") & (replica.replicate == "rep1"),
                       "value"].item() == pytest.approx(0.25)
    system = pd.read_parquet(dest / "system_interaction_metrics.parquet")
    contrast = pd.read_parquet(dest / "system_interaction_comparisons.parquet")
    assert set(system.metric) == set(contrast.metric) == {"occupancy"}


def test_empty_interaction_selection_skips_input_loading(tmp_path, monkeypatch):
    from mdscope.statistics.interactions import analyze_interactions
    import mdscope.statistics.interactions as module
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {"interactions": []}
    assert statistics_family_config(cfg, "interactions") is None
    monkeypatch.setattr(module, "interaction_tables",
                        lambda *args, **kwargs: (_ for _ in ()).throw(
                            AssertionError("No interaction files may be read")))
    analyze_interactions(cfg, tmp_path, {"x": [tmp_path / "x"]})
    assert not (tmp_path / "project_summary/statistics/interactions").exists()


def test_unknown_interaction_identity_fails_before_output_creation(tmp_path):
    from test_interaction_statistics_reference import fixture
    from mdscope.statistics.interactions import analyze_interactions
    cfg, paths = fixture(tmp_path / "input")
    cfg["statistics"]["metric_selection"] = {"interactions": [
        {"kind": "salt_bridges", "interaction": "UNOBSERVED | LYS7 | interchain"}
    ]}
    with pytest.raises(ValueError, match="Unknown saved interaction"):
        analyze_interactions(cfg, tmp_path, paths)
    assert not (tmp_path / "project_summary/statistics/interactions").exists()


@pytest.mark.parametrize("invalid", [
    [{"kind": "salt_bridges", "metric": "invented_stat"}],
    [{"kind": ""}],
    [{"kind": "salt_bridges", "interaction": ""}],
    [{"kind": "salt_bridges"}, {"kind": "salt_bridges"}],
    [{"kind": "salt_bridges", "unsupported": True}],
])
def test_invalid_interaction_selection_fails_config_validation(invalid):
    from mdscope.statistics_selection import validate_metric_selection
    with pytest.raises(ValueError, match="metric_selection"):
        validate_metric_selection({"metric_selection": {"interactions": invalid}})


def test_legacy_interaction_selection_is_unrestricted(tmp_path):
    from test_interaction_statistics_reference import fixture
    from mdscope.statistics.interactions import analyze_interactions
    cfg, paths = fixture(tmp_path / "input")
    analyze_interactions(cfg, tmp_path, paths)
    rows = pd.read_parquet(
        tmp_path / "project_summary/statistics/interactions/replicate_interaction_metrics.parquet"
    )
    assert set(rows.metric) == {
        "occupancy", "event_rate_per_ns", "mean_lifetime_ps",
        "median_lifetime_ps", "total_lifetime_ps"
    }
    assert len(set(rows.interaction)) == 3
