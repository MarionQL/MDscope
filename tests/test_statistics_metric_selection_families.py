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
    cfg["analyses"]["coupling_screen"]["enabled"] = True
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
    cfg["analyses"]["coupling_screen"]["enabled"] = True
    cfg["statistics"]["metric_selection"] = {"coupling": [
        {"mechanism": "nonexistent"}]}
    with pytest.raises(ValueError, match="Unknown configured coupling"):
        analyze_coupling(cfg, tmp_path, {"s": [tmp_path / "missing"]})
