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


def test_interaction_selected_rng_invariant_to_other_requested_identities(tmp_path):
    """Selecting an unrelated pair cannot perturb an existing pair's CIs."""
    from test_interaction_statistics_reference import fixture
    from mdscope.statistics.interactions import analyze_interactions
    cfg, paths = fixture(tmp_path / "input")
    key = "ASP0 | LYS7 | interchain"
    cfg["statistics"]["metric_selection"] = {"interactions": [
        {"kind": "salt_bridges", "interaction": key, "metric": "occupancy"}
    ]}
    analyze_interactions(cfg, tmp_path / "one", paths)
    cfg["statistics"]["metric_selection"]["interactions"].insert(0, {
        "kind": "salt_bridges", "interaction": "ASP1 | LYS7 | interchain",
        "metric": "occupancy",
    })
    analyze_interactions(cfg, tmp_path / "two", paths)
    root = "project_summary/statistics/interactions"
    for filename in (
        "replicate_interaction_metrics.parquet",
        "system_interaction_metrics.parquet",
        "system_interaction_comparisons.parquet",
    ):
        original = pd.read_parquet(tmp_path / "one" / root / filename)
        expanded = pd.read_parquet(tmp_path / "two" / root / filename)
        expanded = expanded.loc[expanded.interaction.eq(key)].reset_index(drop=True)
        # BH q-values correctly depend on the selected multiple-testing
        # universe; point estimates, CIs and raw p-values must not.
        if "q_value" in original:
            original = original.drop(columns="q_value")
            expanded = expanded.drop(columns="q_value")
        pd.testing.assert_frame_equal(original, expanded, check_exact=True)


def test_merged_temporal_cache_and_selective_rng_can_execute_together():
    """Regression: merging #194 into selective inference must not lose selected_seed."""
    from mdscope.statistics.interactions import _interaction_bootstrap
    cache, zeros = {}, {}
    frames = np.arange(16, dtype=float) * 5
    cfg = deepcopy(DEFAULTS)
    rng = np.random.default_rng(34)
    one, meta_a = _interaction_bootstrap(
        frames, np.array([], float), pd.DataFrame(), cfg, 100,
        rng, cache, ("system", "rep1"), zeros, 2026,
    )
    two, meta_b = _interaction_bootstrap(
        frames, np.array([], float), pd.DataFrame(), cfg, 100,
        rng, cache, ("system", "rep1"), zeros, 2026,
    )
    np.testing.assert_equal(meta_a, meta_b)
    assert len(zeros) == 1
    for metric in one:
        np.testing.assert_array_equal(one[metric], two[metric])


def test_descriptive_only_time_series_retains_points_without_bootstrapping(tmp_path, monkeypatch):
    paths = _data(tmp_path / "inputs")
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["bootstrap_draws"] = 50
    cfg["statistics"]["metric_selection"] = {"time_series": [
        {"group": "PROA", "metric": "rmsd_A", "inference": False},
        {"metric": "membrane_thickness_A", "inference": True},
    ]}
    import mdscope.statistics.time_series as module
    original = module._block_bootstrap
    calls = []
    def observed(*args, **kwargs):
        calls.append(args[2])
        return original(*args, **kwargs)
    monkeypatch.setattr(module, "_block_bootstrap", observed)
    analyze_time_series(cfg, tmp_path, paths)
    folder = tmp_path / "project_summary/statistics/time_series"
    rows = pd.read_parquet(folder / "replicate_statistics.parquet")
    descriptive = rows.loc[rows.metric == "rmsd_A"]
    infer = rows.loc[rows.metric == "membrane_thickness_A"]
    assert descriptive.bootstrap_draws.eq(0).all()
    assert descriptive.mean_ci_low.isna().all()
    assert descriptive["mean"].notna().all()
    assert infer.bootstrap_draws.eq(50).all()
    assert infer.mean_ci_low.notna().all()
    assert len(calls) == len(infer) * 2
    systems = pd.read_parquet(folder / "system_statistics.parquet")
    assert systems.loc[systems.metric.eq("rmsd_A"), "mean_ci_low"].isna().all()
    contrasts = pd.read_parquet(folder / "system_comparisons.parquet")
    assert contrasts.loc[contrasts.metric.eq("rmsd_A"), "mean_p_value"].isna().all()


def test_descriptive_only_interaction_has_no_ci_or_p_value(tmp_path):
    from test_interaction_statistics_reference import fixture
    from mdscope.statistics.interactions import analyze_interactions
    cfg, reps = fixture(tmp_path / "inputs")
    cfg["statistics"]["metric_selection"] = {"interactions": [
        {"kind": "salt_bridges", "interaction": "ASP0 | LYS7 | interchain",
         "metric": "occupancy", "inference": False},
    ]}
    analyze_interactions(cfg, tmp_path, reps)
    folder = tmp_path / "project_summary/statistics/interactions"
    rows = pd.read_parquet(folder / "replicate_interaction_metrics.parquet")
    assert set(rows.metric) == {"occupancy"}
    assert rows.bootstrap_draws.eq(0).all()
    assert rows.ci_low.isna().all()
    assert rows.value.notna().any()
    comparisons = pd.read_parquet(folder / "system_interaction_comparisons.parquet")
    assert comparisons.p_value.isna().all()
    assert comparisons.difference_ci_low.isna().all()


@pytest.mark.parametrize("family,entry", [
    ("time_series", {"metric": "rmsd_A", "inference": "yes"}),
    ("interactions", {"kind": "salt_bridges", "inference": None}),
    ("profiles", {"profile": "rmsf", "inference": 1}),
    ("diffusion", {"metric": "DOPC", "inference": []}),
])
def test_per_metric_inference_flag_must_be_boolean(family, entry):
    from mdscope.statistics_selection import validate_metric_selection
    with pytest.raises(ValueError, match="metric_selection"):
        validate_metric_selection({"metric_selection": {family: [entry]}})


def test_selected_time_series_rng_is_stable_when_another_metric_is_added(tmp_path):
    """Same saved identity must keep its raw CIs and p-values across selections."""
    paths = _data(tmp_path / "inputs")
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["bootstrap_draws"] = 40
    cfg["statistics"]["hypothesis_tests"]["permutation_draws"] = 40
    primary = {"group": "PROA", "metric": "rmsd_A"}
    cfg["statistics"]["metric_selection"] = {"time_series": [primary]}
    analyze_time_series(cfg, tmp_path / "only", paths)
    cfg["statistics"]["metric_selection"]["time_series"].append(
        {"metric": "membrane_thickness_A"}
    )
    analyze_time_series(cfg, tmp_path / "both", paths)
    folder = "project_summary/statistics/time_series"
    for name in ("replicate_statistics.parquet", "system_statistics.parquet",
                 "system_comparisons.parquet"):
        one = pd.read_parquet(tmp_path / "only" / folder / name)
        both = pd.read_parquet(tmp_path / "both" / folder / name)
        a = one.loc[one.metric == "rmsd_A"].reset_index(drop=True)
        b = both.loc[both.metric == "rmsd_A"].reset_index(drop=True)
        q = [col for col in a if col.endswith("_q_value")]
        pd.testing.assert_frame_equal(a.drop(columns=q), b.drop(columns=q),
                                      check_exact=True)


@pytest.mark.parametrize("ci_enabled,tests_enabled", [(False, True), (True, False)])
def test_selected_time_series_estimator_independence(
        tmp_path, monkeypatch, ci_enabled, tests_enabled):
    paths = _data(tmp_path / "inputs")
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["bootstrap_draws"] = 50
    cfg["statistics"]["hypothesis_tests"]["permutation_draws"] = 50
    cfg["statistics"]["metric_selection"] = {"time_series": [
        {"group": "PROA", "metric": "rmsd_A",
         "ci": ci_enabled, "hypothesis_tests": tests_enabled},
    ]}
    if not ci_enabled:
        import mdscope.statistics.time_series as module
        monkeypatch.setattr(module, "_block_bootstrap",
                            lambda *args, **kwargs: (_ for _ in ()).throw(
                                AssertionError("CI-disabled metric must not bootstrap")))
    analyze_time_series(cfg, tmp_path, paths)
    folder = tmp_path / "project_summary/statistics/time_series"
    replicas = pd.read_parquet(folder / "replicate_statistics.parquet")
    comparisons = pd.read_parquet(folder / "system_comparisons.parquet")
    assert replicas.mean_ci_low.notna().any() == ci_enabled
    assert comparisons.mean_p_value.notna().any() == tests_enabled
    assert replicas.bootstrap_draws.eq(50 if ci_enabled else 0).all()


@pytest.mark.parametrize("family,entry", [
    ("solvation_rdfs", {"solvent": "water", "target_class": "protein",
                       "site_class": "acidic", "hypothesis_tests": True}),
    ("diffusion", {"metric": "DOPC", "hypothesis_tests": False}),
    ("coupling", {"mechanism": "membrane_coupling", "hypothesis_tests": True}),
    ("time_series", {"metric": "rmsd_A", "inference": False, "ci": True}),
])
def test_unsupported_or_contradictory_estimator_flags_rejected(family, entry):
    from mdscope.statistics_selection import validate_metric_selection
    with pytest.raises(ValueError, match="metric_selection"):
        validate_metric_selection({"metric_selection": {family: [entry]}})


def test_mixed_interaction_inference_draw_counts_do_not_share_incompatible_weights():
    from mdscope.statistics.interactions import _interaction_bootstrap
    times = np.arange(24, dtype=float) * 5
    cfg = deepcopy(DEFAULTS)
    cache, zero_cache = {}, {}
    rng = np.random.default_rng(2026)
    for draws in (0, 32, 0, 32):
        boot, point = _interaction_bootstrap(
            times, np.array([10., 25.]), pd.DataFrame(),
            cfg, draws, rng, cache, ("system", "rep1"), zero_cache, 2026,
        )
        assert len(boot["occupancy"]) == draws
        assert point["occupancy"] == pytest.approx(2 / len(times))
    assert sorted({key[-1] for key in cache}) == [0, 32]


def test_changing_metric_selection_rebuilds_statistics_without_touching_raw_data(tmp_path):
    """Re-finalization must not mistake old selected outputs for current ones."""
    import hashlib
    from mdscope.statistics import analyze_project
    paths = _data(tmp_path / "inputs")
    raw = {
        source: hashlib.sha256(source.read_bytes()).hexdigest()
        for source in (tmp_path / "inputs").rglob("*.parquet")
    }
    cfg = deepcopy(DEFAULTS)
    cfg["progress"] = {"enabled": False}
    cfg["statistics"]["bootstrap_draws"] = 20
    cfg["statistics"]["hypothesis_tests"]["enabled"] = False
    cfg["statistics"]["families"] = {
        name: False for name in
        ("profiles", "interactions", "solvation_rdfs", "diffusion", "coupling")
    }
    cfg["statistics"]["metric_selection"] = {"time_series": [
        {"group": "PROA", "metric": "rmsd_A"}
    ]}
    # The finalizer's root is the project containing these replicas.
    root = tmp_path / "inputs"
    analyze_project(cfg, root, paths)
    initial = pd.read_parquet(
        root / "project_summary/statistics/time_series/replicate_statistics.parquet"
    )
    assert set(initial.metric) == {"rmsd_A"}
    cfg["statistics"]["metric_selection"] = {"time_series": [
        {"metric": "membrane_thickness_A"}
    ]}
    analyze_project(cfg, root, paths)
    updated = pd.read_parquet(
        root / "project_summary/statistics/time_series/replicate_statistics.parquet"
    )
    assert set(updated.metric) == {"membrane_thickness_A"}
    assert all(hashlib.sha256(path.read_bytes()).hexdigest() == digest
               for path, digest in raw.items())
    import json
    manifest = json.loads(
        (root / "project_summary/statistics/manifest.json").read_text()
    )
    assert manifest["status"] == "complete"
    assert manifest["methods"]["metric_selection"]["time_series"][0]["metric"] == (
        "membrane_thickness_A"
    )


@pytest.mark.parametrize("raw_metric,derived_metric", [
    ("sdc4_local_thickness_A", "sdc4_local_thickness_A_available_fraction"),
    ("sdc4_motif_A_glycan_contacts", "sdc4_motif_A_hs_contacts_during_contact"),
    ("sdc4_motif_A_glycan_contacts", "sdc4_motif_A_hs_any_contact_state"),
])
def test_selected_derived_time_series_loads_only_source_columns(
        tmp_path, monkeypatch, raw_metric, derived_metric):
    from mdscope.statistics.time_series import collect_time_series
    path = tmp_path / "summary" / "timeseries.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame({
        "time_ps": [0., 1., 2.],
        raw_metric: [0., 1., 2.],
        "unrelated_expensive_observable": [9., 9., 9.],
    }).to_parquet(path)
    original = pd.read_parquet
    observed = []

    def tracked(source, *args, **kwargs):
        observed.append(kwargs.get("columns"))
        return original(source, *args, **kwargs)

    monkeypatch.setattr(pd, "read_parquet", tracked)
    result = collect_time_series(tmp_path, selected={(None, derived_metric)})
    assert (None, derived_metric) in result
    assert observed == [["time_ps", raw_metric]]
    assert "unrelated_expensive_observable" not in result


def test_selected_inference_respects_supplied_system_and_replica_population(tmp_path):
    paths = _data(tmp_path / "inputs")
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["bootstrap_draws"] = 30
    cfg["statistics"]["hypothesis_tests"]["enabled"] = False
    cfg["statistics"]["metric_selection"] = {
        "time_series": [{"group": "PROA", "metric": "rmsd_A"}]
    }
    analyze_time_series(cfg, tmp_path / "single", {"a": [paths["a"][0]]})
    folder = tmp_path / "single/project_summary/statistics/time_series"
    replicas = pd.read_parquet(folder / "replicate_statistics.parquet")
    assert len(replicas) == 1
    assert set(replicas.system) == {"a"}
    assert set(replicas.metric) == {"rmsd_A"}
    systems = pd.read_parquet(folder / "system_statistics.parquet")
    assert set(systems.system) == {"a"}
    assert not (folder / "system_comparisons.parquet").exists()
