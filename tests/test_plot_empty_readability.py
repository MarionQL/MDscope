"""Renderer regressions for unavailable data and label-heavy matrices.

No trajectory files or inferential estimators are used.
"""
from pathlib import Path

import numpy as np
import pandas as pd


def _cfg():
    return {"plots": {"enabled": True, "dpi": 72},
            "systems": {"control": {"label": "Control"}}}


def test_lipid_contact_bootstrap_skips_all_missing_matrix(tmp_path):
    from mdscope.reporting.lipid_contact_bootstrap import _draw

    group = pd.DataFrame([{"species_a": "POPC", "species_b": "POPC",
                           "mean": np.nan}])
    _draw(_cfg(), tmp_path, group, ordered=["POPC"], value="mean",
          supported=False, title="Missing", stem="lipid_contact_matrix_missing",
          metric="enrichment", limit=1., comparison=False)
    assert not list(tmp_path.glob("*.png"))


def test_lipid_contact_bootstrap_keeps_finite_zero_and_long_labels(tmp_path):
    from mdscope.reporting.lipid_contact_bootstrap import _draw

    group = pd.DataFrame([{
        "species_a": "phosphatidylethanolamine",
        "species_b": "phosphatidylethanolamine", "mean": 0.,
    }])
    _draw(_cfg(), tmp_path, group, ordered=["phosphatidylethanolamine"],
          value="mean", supported=False, title="Real measured zero",
          stem="lipid_contact_class_matrix_zero", metric="pair_count_undirected",
          limit=1., comparison=False)
    assert (tmp_path / "lipid_contact_class_matrix_zero.png").is_file()


def test_residue_response_skips_all_missing_effects(tmp_path):
    from mdscope.reporting.residue_response import plot_residue_response

    target = tmp_path / "project_summary/statistics/profiles"
    target.mkdir(parents=True)
    pd.DataFrame([
        dict(system_1="control", system_2="treated", chain="PROA",
             profile="water_coordination", resid=r,
             standardized_effect=np.nan, ci_excludes_zero=False, region_id=None)
        for r in (201, 202, 203)
    ]).to_parquet(target / "residue_response_matrix.parquet")
    plot_residue_response(_cfg(), tmp_path)
    assert not list((tmp_path / "project_summary/figures").rglob("*.png"))


def test_time_series_skips_all_nonfinite_metrics(tmp_path, monkeypatch):
    import mdscope.reporting.plotting as plotting

    monkeypatch.setattr(plotting, "collect_time_series", lambda *_args, **_kwargs: {
        (None, "absent_metric"): pd.DataFrame(
            {"time_ps": [0., 100.], "value": [np.nan, np.nan]})
    })
    rep = tmp_path / "control/rep1"
    plotting.plot_time_series(_cfg(), tmp_path, {"control": [rep]})
    assert not list((tmp_path / "project_summary/figures").rglob("*.png"))


def test_profiles_skip_all_missing_values(tmp_path):
    from mdscope.reporting.plotting import plot_profiles

    target = tmp_path / "project_summary/statistics/profiles"
    target.mkdir(parents=True)
    pd.DataFrame([dict(profile="water_coordination", system="control",
                       replicate="rep1", position=10, value=np.nan)]
                 ).to_parquet(target / "replicate_profiles.parquet")
    pd.DataFrame([dict(profile="water_coordination", system="control",
                       position=10, mean=np.nan)]
                 ).to_parquet(target / "system_profiles.parquet")
    plot_profiles(_cfg(), tmp_path, {"control": [tmp_path / "control/rep1"]})
    assert not list((tmp_path / "project_summary/figures").rglob("*.png"))


def test_diffusion_skips_all_missing_coefficients(tmp_path):
    from mdscope.reporting.membrane_figures import plot_project_diffusion

    target = tmp_path / "project_summary/diffusion"
    target.mkdir(parents=True)
    pd.DataFrame([dict(system="control", metric="POPC",
                       diffusion_nm2_ns=np.nan)]
                 ).to_parquet(target / "replicate_coefficients.parquet")
    plot_project_diffusion(_cfg(), tmp_path)
    assert not (target / "diffusion_coefficients.png").exists()
