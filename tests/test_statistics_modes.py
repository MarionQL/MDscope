"""Point-only modes, incomplete profiles, and non-destructive statistical reruns."""
from copy import deepcopy
import json
import hashlib
import numpy as np
import pandas as pd
import pytest
from mdscope.statistics import analyze_project
from mdscope.statistics.time_series import analyze_time_series
from mdscope.statistics.profiles import analyze_profiles, _prepare_profile, _common_profile_axis
from mdscope.statistics.interactions import analyze_interactions
from mdscope.reporting.figure_inventory import current_statistics
from test_interaction_statistics_reference import fixture


def _profiles(paths, incomplete=False):
    result = {}
    for system, reps in paths.items():
        for i, rep in enumerate(reps):
            times = np.arange(6+i, dtype=float)
            values = np.column_stack([times+i, times+i+10.]) + (5. if system == "b" else 0.)
            if incomplete and system == "b" and i == 1:
                values[2,0] = np.nan
            labels = np.array(["L", "R"])
            positions = np.array([0.,1.])
            if i == 1:
                values, labels, positions = values[:,::-1], labels[::-1], positions[::-1]
            result[rep] = {"profile": {"times": times, "values": values, "positions": positions, "labels": labels}}
    return result


@pytest.mark.parametrize("family", ["time_series", "profiles", "interactions"])
@pytest.mark.parametrize("mode", ["point_only", "no_hierarchy"])
def test_independent_bootstrap_hierarchy_and_test_controls(tmp_path, monkeypatch, family, mode):
    cfg, paths = fixture(tmp_path/"input")
    option = {"bootstrap": False, "hypothesis_tests": False} if mode == "point_only" else {"hierarchy": {"enabled": False}}
    cfg["statistics"]["families"] = {family: option}
    profiles = _profiles(paths)
    import mdscope.statistics.profiles as profile_module
    monkeypatch.setattr(profile_module, "collect_temporal_profiles", lambda rep, cfg: profiles[rep])
    for system, reps in paths.items():
        for i, rep in enumerate(reps):
            file = rep/"summary/timeseries.parquet"
            data = pd.read_parquet(file)
            data["value_A"] = np.arange(len(data)) + i + (5 if system == "b" else 0)
            data.to_parquet(file)
    original = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in (tmp_path/"input").rglob("*.parquet")}
    if mode == "point_only":
        def forbidden(*a, **k): raise AssertionError("disabled bootstrap kernel called")
        for module, name in (("time_series", "_block_bootstrap"), ("profiles", "_block_bootstrap"),
                             ("profiles", "_block_weights"), ("interactions", "_block_weights")):
            monkeypatch.setattr(f"mdscope.statistics.{module}.{name}", forbidden)
    {"time_series": analyze_time_series, "profiles": analyze_profiles, "interactions": analyze_interactions}[family](cfg, tmp_path, paths)
    folder = tmp_path/"project_summary/statistics"/family
    system_file, point, ci, compare = {
        "time_series": ("system_statistics.parquet", "mean", "mean_ci_low", "system_comparisons.parquet"),
        "profiles": ("system_profiles.parquet", "mean", "mean_ci_low", "system_comparisons.parquet"),
        "interactions": ("system_interaction_metrics.parquet", "mean", "mean_ci_low", "system_interaction_comparisons.parquet"),
    }[family]
    system = pd.read_parquet(folder/system_file)
    assert system[point].notna().any() and system[ci].isna().all()
    assert system.bootstrap_draws.eq(0).all()
    comp = pd.read_parquet(folder/compare)
    if mode == "point_only":
        assert comp[[c for c in comp if c.endswith("p_value")]].isna().all().all()
    else:
        assert comp[[c for c in comp if c.endswith("p_value")]].notna().any().any()
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == value for p,value in original.items())


@pytest.mark.parametrize("damage", ["gap", "missing_value", "infinite_value"])
def test_incomplete_profile_cannot_produce_temporal_confidence_intervals(damage):
    cfg = {"statistics": {"bootstrap_draws": 100}}
    times = np.arange(6, dtype=float)
    values = np.column_stack([times, times+10.])
    if damage == "gap": times[3:] += 2.
    else: values[2,0] = np.inf if damage == "infinite_value" else np.nan
    item = _prepare_profile({"times": times, "values": values, "positions": [0,1], "labels": ["L","R"]}, cfg, 100, np.random.default_rng(8))
    assert np.isnan(item["draws"][:,0]).all()
    if damage != "gap":
        assert np.isfinite(item["draws"][:,1]).all()
        assert item["point"][0] == pytest.approx(2.6)
    else: assert np.isnan(item["draws"]).all() and np.isnan(item["temporal"]["effective_n"])


def test_profile_alignment_and_incomplete_replica_are_not_silently_dropped(tmp_path, monkeypatch):
    cfg, paths = fixture(tmp_path/"input")
    profiles = _profiles(paths, incomplete=True)
    import mdscope.statistics.profiles as module
    monkeypatch.setattr(module, "collect_temporal_profiles", lambda rep,cfg: profiles[rep])
    analyze_profiles(cfg, tmp_path, paths)
    folder = tmp_path/"project_summary/statistics/profiles"
    system = pd.read_parquet(folder/"system_profiles.parquet")
    b = system[system.system.eq("b")].set_index("position_label")
    assert np.isnan(b.loc["L", "mean_ci_low"])
    assert np.isfinite(b.loc["R", "mean_ci_low"])
    assert b.loc["R", "mean"] == 19.  # replica means 17.5,19,20.5, equally weighted
    comp = pd.read_parquet(folder/"system_comparisons.parquet").set_index("position_label")
    assert np.isnan(comp.loc["L", "difference_ci_low"])
    assert np.isfinite(comp.loc["L", "difference_system2_minus_system1"])
    assert np.isfinite(comp.loc["R", "difference_ci_low"])
    first = {"positions": np.array([7.,9.]), "labels": np.array(["x","y"])}
    second = {"positions": np.array([9.,7.]), "labels": np.array(["y","x"])}
    pos, labels, indices = _common_profile_axis([first, second])
    np.testing.assert_array_equal(pos, [7.,9.]); np.testing.assert_array_equal(labels, ["x","y"])
    np.testing.assert_array_equal(indices[1], [1,0])


def test_paired_profiles_with_no_complete_pairs_are_unavailable_not_an_error(tmp_path, monkeypatch):
    cfg, paths = fixture(tmp_path/"input", design="paired")
    profiles = _profiles(paths)
    import mdscope.statistics.profiles as module
    def available(rep,cfg):
        return profiles[rep] if (rep.parent.name,rep.name) in {("a","rep1"),("b","rep2")} else {}
    monkeypatch.setattr(module,"collect_temporal_profiles",available)
    analyze_profiles(cfg,tmp_path,paths)
    folder = tmp_path/"project_summary/statistics/profiles"
    system = pd.read_parquet(folder/"system_profiles.parquet")
    assert system["mean"].notna().all() and system.mean_ci_low.isna().all()
    comp = pd.read_parquet(folder/"system_comparisons.parquet")
    assert comp.difference_ci_low.isna().all() and comp.p_value.isna().all()
    assert comp.difference_system2_minus_system1.isna().all()


def test_reruns_archive_empty_or_disabled_outputs_and_reject_partial_statistics(tmp_path, monkeypatch):
    cfg, paths = fixture(tmp_path/"input")
    cfg["statistics"]["families"] = {"profiles": False, "interactions": False, "solvation_rdfs": False}
    for reps in paths.values():
        for rep in reps:
            file = rep/"summary/timeseries.parquet"
            data = pd.read_parquet(file); data["value_A"] = np.arange(len(data)); data.to_parquet(file)
    analyze_project(cfg, tmp_path, paths)
    assert current_statistics(cfg,tmp_path,"time_series")
    cfg["statistics"]["families"]["time_series"] = False
    analyze_project(cfg,tmp_path,paths)
    assert not (tmp_path/"project_summary/statistics/time_series").exists()
    assert list((tmp_path/".archives/statistics/time_series").rglob("system_statistics.parquet"))
    cfg["statistics"]["families"]["time_series"] = True
    import mdscope.statistics as module
    def crash(*a, **k): raise RuntimeError("injected statistical failure")
    monkeypatch.setattr(module,"analyze_profiles",crash)
    with pytest.raises(RuntimeError,match="injected"): analyze_project(cfg,tmp_path,paths)
    assert not current_statistics(cfg,tmp_path,"time_series")
    assert json.loads((tmp_path/"project_summary/statistics/manifest.json").read_text())["status"] == "running"
