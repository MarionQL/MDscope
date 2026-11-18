"""Issue 179: selected statistics families preserve canonical scientific outputs."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from mdscope.config import DEFAULTS, statistics_family_config
from mdscope.statistics_selection import validate_metric_selection
from mdscope.statistics.profiles import analyze_profiles
from mdscope.statistics.solvation import analyze_solvation_rdfs
from mdscope.statistics.diffusion import analyze_multi_origin_diffusion
from mdscope.statistics.coupling import analyze_coupling
from test_interaction_statistics_reference import fixture
from test_statistics_modes import _profiles


@pytest.mark.parametrize("family,entry", [
    ("profiles", {"profile": "profile"}),
    ("solvation_rdfs", {"solvent": "water", "target_class": "protein",
                       "site_class": "acidic"}),
    ("diffusion", {"metric": "DOPC"}),
    ("coupling", {"mechanism": "membrane_coupling"}),
])
def test_simple_selector_validation_and_explicit_empty(family, entry):
    validate_metric_selection({"metric_selection": {family: [entry]}})
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {family: []}
    assert statistics_family_config(cfg, family) is None
    validate_metric_selection({"metric_selection": {family: []}})


@pytest.mark.parametrize("family,entry", [
    ("profiles", {"profile": ""}),
    ("profiles", {"profile": "a", "extra": "b"}),
    ("diffusion", {"metric": 4}),
    ("coupling", {"mechanism": []}),
    ("solvation_rdfs", {"solvent": "water", "target_class": "protein"}),
])
def test_invalid_simple_selector_rejected(family, entry):
    with pytest.raises(ValueError, match="metric_selection"):
        validate_metric_selection({"metric_selection": {family: [entry]}})


def test_selected_profile_skips_other_temporal_kernels(tmp_path, monkeypatch):
    cfg, paths = fixture(tmp_path / "input")
    sources = _profiles(paths)
    for rep in sources:
        sources[rep]["other"] = sources[rep]["profile"]
    import mdscope.statistics.profiles as mod
    monkeypatch.setattr(mod, "collect_temporal_profiles",
                        lambda rep, config: sources[rep])
    cfg["statistics"]["metric_selection"] = {"profiles": [{"profile": "profile"}]}
    analyze_profiles(cfg, tmp_path, paths)
    dest = tmp_path / "project_summary/statistics/profiles"
    df = pd.read_parquet(dest / "replicate_profiles.parquet")
    assert set(df.profile) == {"profile"}
    assert len(df) > 0


def test_missing_profile_selector_fails_before_any_output(tmp_path, monkeypatch):
    cfg, paths = fixture(tmp_path / "input")
    sources = _profiles(paths)
    import mdscope.statistics.profiles as mod
    monkeypatch.setattr(mod, "collect_temporal_profiles",
                        lambda rep, config: sources[rep])
    cfg["statistics"]["metric_selection"] = {"profiles": [{"profile": "typo"}]}
    with pytest.raises(ValueError, match="Unknown saved profile"):
        analyze_profiles(cfg, tmp_path, paths)
    assert not (tmp_path / "project_summary/statistics/profiles").exists()


def test_selected_rdf_preserves_only_requested_chemical_identity(tmp_path):
    paths = {"s": []}
    for i in range(2):
        rep = tmp_path / "s" / f"rep{i}"
        source = rep / "extensions/solvation/rdf.parquet"
        source.parent.mkdir(parents=True)
        pd.DataFrame({
            "solvent": ["water"] * 3 + ["ions/K"] * 3,
            "target_class": ["protein"] * 6,
            "site_class": ["acidic"] * 6,
            "r_inner_A": [0., 1., 2.] * 2,
            "r_outer_A": [1., 2., 3.] * 2,
            "g_bulk": [1., 1.5, 2.] * 2,
        }).to_parquet(source)
        paths["s"].append(rep)
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {"solvation_rdfs": [
        {"solvent": "water", "target_class": "protein", "site_class": "acidic"}
    ]}
    analyze_solvation_rdfs(cfg, tmp_path, paths)
    table = pd.read_parquet(
        tmp_path / "project_summary/statistics/solvation/system_rdf.parquet"
    )
    assert set(table.solvent) == {"water"}
    assert table.complete_replica_coverage.all()


def test_unknown_rdf_identity_fails_without_products(tmp_path):
    rep = tmp_path / "s/rep"
    path = rep / "extensions/solvation/rdf.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame({
        "solvent": ["water"], "target_class": ["protein"],
        "site_class": ["acidic"], "r_inner_A": [0.],
        "r_outer_A": [1.], "g_bulk": [1.],
    }).to_parquet(path)
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {"solvation_rdfs": [
        {"solvent": "ions/K", "target_class": "protein", "site_class": "acidic"}
    ]}
    with pytest.raises(ValueError, match="Unknown saved solvation RDF"):
        analyze_solvation_rdfs(cfg, tmp_path, {"s": [rep]})
    assert not (tmp_path / "project_summary/statistics/solvation").exists()


def test_coupling_selection_filters_before_collecting_unrelated_columns(tmp_path, monkeypatch):
    cfg = deepcopy(DEFAULTS)
    cfg["analyses"]["coupling_screen"] = {"enabled": True, "block_frames": 50}
    cfg["statistics"]["metric_selection"] = {"coupling": [
        {"mechanism": "membrane_coupling"}]}
    rep = tmp_path / "s/rep1"
    path = rep / "summary/timeseries.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame({
        "time_ps": np.arange(500, dtype=float),
        "local_protein_lipid_contact_residues_count": np.arange(500, dtype=float),
        "membrane_thickness_A": 30 + np.arange(500, dtype=float) * .1,
    }).to_parquet(path)
    analyze_coupling(cfg, tmp_path, {"s": [rep]})
    result = pd.read_parquet(
        tmp_path / "project_summary/statistics/coupling/replicate_correlations.parquet"
    )
    assert set(result.mechanism) == {"membrane_coupling"}


def test_missing_coupling_is_rejected_before_read(tmp_path, monkeypatch):
    cfg = deepcopy(DEFAULTS)
    cfg["analyses"]["coupling_screen"] = {"enabled": True, "block_frames": 50}
    cfg["statistics"]["metric_selection"] = {"coupling": [
        {"mechanism": "nonexistent"}]}
    with pytest.raises(ValueError, match="Unknown configured coupling"):
        analyze_coupling(cfg, tmp_path, {"s": [tmp_path / "missing"]})


def test_diffusion_selection_filters_shared_lag_fits(tmp_path, monkeypatch):
    cfg = deepcopy(DEFAULTS)
    cfg["analyses"]["membrane"]["msd"].update(enabled=True, multi_origin=True)
    cfg["statistics"]["metric_selection"] = {"diffusion": [{"metric": "DOPC"}]}
    replicas = {"a": [], "b": []}
    for system in replicas:
        for i in range(2):
            rep = tmp_path / system / f"rep{i}"
            source = rep / "summary/multi_origin_msd.parquet"
            source.parent.mkdir(parents=True)
            lag = np.arange(30, dtype=float)
            pd.DataFrame({
                "metric": ["DOPC"] * 30 + ["POPC"] * 30,
                "lag_ps": list(lag) * 2,
                "msd_A2": list(lag * (2 + i)) + list(lag * 3),
                "n_time_origins": [100] * 60,
            }).to_parquet(source)
            replicas[system].append(rep)
    import mdscope.statistics.diffusion as mod
    original_fit = mod.fit_multi_origin_curve
    observed = []
    def counted(*args, **kwargs):
        observed.append(1)
        return original_fit(*args, **kwargs)
    monkeypatch.setattr(mod, "fit_multi_origin_curve", counted)
    analyze_multi_origin_diffusion(cfg, tmp_path, replicas)
    result = pd.read_parquet(
        tmp_path / "project_summary/statistics/diffusion/replicate_common_lag_slopes.parquet"
    )
    assert set(result.metric) == {"DOPC"}
    assert len(observed) == 4


def test_unknown_diffusion_metric_fails_before_products(tmp_path):
    cfg = deepcopy(DEFAULTS)
    cfg["analyses"]["membrane"]["msd"].update(enabled=True, multi_origin=True)
    cfg["statistics"]["metric_selection"] = {"diffusion": [{"metric": "NOT_SAVED"}]}
    rep = tmp_path / "a/rep1"
    source = rep / "summary/multi_origin_msd.parquet"
    source.parent.mkdir(parents=True)
    pd.DataFrame({
        "metric": ["DOPC"], "lag_ps": [0.], "msd_A2": [0.],
        "n_time_origins": [100],
    }).to_parquet(source)
    with pytest.raises(ValueError, match="Unknown saved diffusion metric"):
        analyze_multi_origin_diffusion(cfg, tmp_path, {"a": [rep]})
    assert not (tmp_path / "project_summary/statistics/diffusion").exists()


def test_studio_new_config_explicitly_disables_all_inference_families():
    from mdscope.statistics_selection import studio_metric_selection_defaults
    selection = studio_metric_selection_defaults()
    assert set(selection) == {
        "time_series", "profiles", "interactions",
        "solvation_rdfs", "diffusion", "coupling",
    }
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = selection
    validate_metric_selection(cfg["statistics"])
    assert all(statistics_family_config(cfg, family) is None for family in selection)
    assert statistics_family_config(deepcopy(DEFAULTS), "time_series") is not None


def test_studio_metric_selection_defaults_are_not_shared_mutable_lists():
    from mdscope.statistics_selection import studio_metric_selection_defaults
    selection = studio_metric_selection_defaults()
    selection["profiles"].append({"profile": "test"})
    assert studio_metric_selection_defaults()["profiles"] == []


def test_selected_diffusion_requires_multi_origin_enabled(tmp_path):
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {"diffusion": [{"metric": "DOPC"}]}
    with pytest.raises(ValueError, match="multi_origin"):
        analyze_multi_origin_diffusion(cfg, tmp_path, {"a": []})


def test_selected_coupling_requires_analysis_enabled(tmp_path):
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {"coupling": [
        {"mechanism": "membrane_coupling"}]}
    with pytest.raises(ValueError, match="coupling_screen.enabled"):
        analyze_coupling(cfg, tmp_path, {"a": []})


def test_descriptive_only_profile_retains_position_points_without_ci(tmp_path, monkeypatch):
    cfg, reps = fixture(tmp_path / "input")
    import mdscope.statistics.profiles as module
    sources = _profiles(reps)
    monkeypatch.setattr(module, "collect_temporal_profiles",
                        lambda rep, config: sources[rep])
    cfg["statistics"]["metric_selection"] = {
        "profiles": [{"profile": "profile", "inference": False}]
    }
    monkeypatch.setattr(module, "_block_bootstrap",
                        lambda *args, **kwargs: (_ for _ in ()).throw(
                            AssertionError("No bootstrap expected")))
    analyze_profiles(cfg, tmp_path, reps)
    folder = tmp_path / "project_summary/statistics/profiles"
    result = pd.read_parquet(folder / "system_profiles.parquet")
    assert result.bootstrap_draws.eq(0).all()
    assert result["mean"].notna().all()
    assert result.mean_ci_low.isna().all()
    if (folder / "system_comparisons.parquet").exists():
        comparisons = pd.read_parquet(folder / "system_comparisons.parquet")
        assert comparisons.p_value.isna().all()


def test_descriptive_only_rdf_skips_shell_comparison(tmp_path):
    rep = tmp_path / "s/rep1"
    path = rep / "extensions/solvation/rdf.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame({
        "solvent": ["water", "water"],
        "target_class": ["protein", "protein"],
        "site_class": ["acidic", "acidic"],
        "r_inner_A": [0., 1.], "r_outer_A": [1., 2.],
        "g_bulk": [1.2, 1.4],
    }).to_parquet(path)
    cfg = deepcopy(DEFAULTS)
    cfg["statistics"]["metric_selection"] = {"solvation_rdfs": [
        {"solvent": "water", "target_class": "protein",
         "site_class": "acidic", "inference": False}
    ]}
    analyze_solvation_rdfs(cfg, tmp_path, {"s": [rep]})
    directory = tmp_path / "project_summary/statistics/solvation"
    assert (directory / "system_rdf.parquet").is_file()
    assert not (directory / "rdf_shell_comparisons.parquet").exists()


def test_descriptive_only_coupling_skips_replica_contrasts(tmp_path):
    cfg = deepcopy(DEFAULTS)
    cfg["analyses"]["coupling_screen"] = {"enabled": True, "block_frames": 20}
    cfg["statistics"]["metric_selection"] = {"coupling": [
        {"mechanism": "membrane_coupling", "inference": False}
    ]}
    paths = {"s": []}
    for i in range(2):
        rep = tmp_path / "s" / f"rep{i}"
        path = rep / "summary/timeseries.parquet"
        path.parent.mkdir(parents=True)
        pd.DataFrame({
            "time_ps": np.arange(200, dtype=float),
            "local_protein_lipid_contact_residues_count": np.arange(200, dtype=float) + i,
            "membrane_thickness_A": np.linspace(30, 33, 200),
        }).to_parquet(path)
        paths["s"].append(rep)
    analyze_coupling(cfg, tmp_path, paths)
    folder = tmp_path / "project_summary/statistics/coupling"
    assert (folder / "replicate_correlations.parquet").exists()
    assert not (folder / "system_comparisons.parquet").exists()


def test_profile_selected_rng_stable_under_extra_profiles(tmp_path, monkeypatch):
    cfg, reps = fixture(tmp_path / "input")
    cfg["statistics"]["bootstrap_draws"] = 30
    import mdscope.statistics.profiles as module
    profiles = _profiles(reps)
    for rep in profiles:
        profiles[rep]["other"] = profiles[rep]["profile"]
    monkeypatch.setattr(module, "collect_temporal_profiles",
                        lambda rep, config: profiles[rep])
    cfg["statistics"]["metric_selection"] = {"profiles": [{"profile": "profile"}]}
    analyze_profiles(cfg, tmp_path / "first", reps)
    cfg["statistics"]["metric_selection"]["profiles"].append({"profile": "other"})
    analyze_profiles(cfg, tmp_path / "second", reps)
    folder = "project_summary/statistics/profiles"
    for file in ("replicate_profile_statistics.parquet",
                 "system_profiles.parquet", "system_comparisons.parquet"):
        path_a = tmp_path / "first" / folder / file
        if not path_a.is_file():
            continue
        a = pd.read_parquet(path_a)
        b = pd.read_parquet(tmp_path / "second" / folder / file)
        a = a.loc[a.profile == "profile"].reset_index(drop=True)
        b = b.loc[b.profile == "profile"].reset_index(drop=True)
        if "q_value" in a.columns:
            a = a.drop(columns=["q_value"])
            b = b.drop(columns=["q_value"])
        pd.testing.assert_frame_equal(a, b, check_exact=True)
