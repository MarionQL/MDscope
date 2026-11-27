"""Renderer regressions for unavailable data and label-heavy matrices.

No trajectory files or inferential estimators are used.
"""

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
    _draw(_cfg(), tmp_path / "project_summary/figures/lipid_contacts", group, ordered=["phosphatidylethanolamine"],
          value="mean", supported=False, title="Real measured zero",
          stem="lipid_contact_class_matrix_zero_upper_pair_count_undirected", metric="pair_count_undirected",
          limit=1., comparison=False)
    assert (tmp_path / "project_summary/figures/lipid_contacts/secondary/lipid_contacts@lipid_contact_class_matrix_zero_upper_pair_count_undirected.png").is_file()


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


def test_solvation_rdf_skips_all_missing_curves(tmp_path):
    from mdscope.reporting.solvation_figures import solvation_rdf_figure

    path = tmp_path / "project_summary/statistics/solvation"
    path.mkdir(parents=True)
    pd.DataFrame([dict(system="control", solvent="water",
                       target_class="protein", site_class="acidic",
                       r_inner_A=2., r_outer_A=3.,
                       mean_g_bulk=np.nan, sd_g_bulk=np.nan,
                       n_replicates=3)]).to_parquet(path / "system_rdf.parquet")
    assert solvation_rdf_figure(tmp_path, {}) is None


def test_solvation_timeseries_skips_all_missing_counts(tmp_path):
    from mdscope.reporting.solvation_figures import solvation_timeseries_figure

    rep = tmp_path / "control/rep1"
    path = rep / "extensions/solvation"
    path.mkdir(parents=True)
    pd.DataFrame([dict(time_ps=100., solvent="water",
                       target_class="protein", site_class="acidic",
                       mean_n_solvent=np.nan)]).to_parquet(
                           path / "species_timeseries.parquet")
    assert solvation_timeseries_figure({"control": [rep]}) is None

def test_all_zero_lipid_footprint_is_suppressed_and_old_figure_archived(tmp_path, monkeypatch):
    from mdscope.reporting import plotting

    metric = "lipid_footprint_BGAL_lower_5_10A_fold_available_fraction"
    cfg = {
        "plots": {"enabled": True, "distributions": False,
                  "autocorrelation": False, "bootstrap_tables": False,
                  "time_series_rolling_window_ns": 100},
        "systems": {
            "synd4_glycolipids": {"replicates": [{"label": f"rep{i}"} for i in (1, 2, 3)]},
            "synd4_no_glycolipids": {"replicates": [{"label": f"rep{i}"} for i in (1, 2, 3)]},
        },
    }
    reps = {name: [tmp_path / name / f"rep{i}" for i in (1, 2, 3)]
            for name in cfg["systems"]}
    def collect(*args, **kwargs):
        return {(None, metric): pd.DataFrame({
            "time_ps": [0., 100_000., 200_000.],
            "value": [0., 0., 0.],
        })}
    monkeypatch.setattr(plotting, "collect_time_series", collect)
    dest = tmp_path / "project_summary/figures/time_series"
    dest.mkdir(parents=True)
    stale = dest / f"{metric}__timeseries.png"
    stale.write_bytes(b"older rendered figure")
    plotting.plot_time_series(cfg, tmp_path, reps)
    assert not stale.exists()
    archived = list((tmp_path / ".archives/zero_only_footprint_figures").glob(
        f"{metric}__timeseries*.png"))
    assert len(archived) == 1
    assert archived[0].read_bytes() == b"older rendered figure"


def test_lipid_footprint_one_nonzero_replica_keeps_figure(tmp_path, monkeypatch):
    from mdscope.reporting import plotting

    metric = "lipid_footprint_POPC_lower_5_10A_fold_available_fraction"
    cfg = {"plots": {"enabled": True, "distributions": False,
                     "autocorrelation": False, "bootstrap_tables": False},
           "systems": {"control": {"replicates": [{"label": "rep1"}, {"label": "rep2"}]}}}
    def collect(repout, **kwargs):
        values = [0., 0., 1.] if repout.name == "rep2" else [0., 0., 0.]
        return {(None, metric): pd.DataFrame({
            "time_ps": [0., 100., 200.], "value": values,
        })}
    monkeypatch.setattr(plotting, "collect_time_series", collect)
    reps = {"control": [tmp_path / "control/rep1", tmp_path / "control/rep2"]}
    plotting.plot_time_series(cfg, tmp_path, reps)
    from mdscope.registry import figure_descriptor
    assert (tmp_path / figure_descriptor(_cfg(), f"time_series:{metric}:timeseries")["destination"]).is_file()


def test_zero_only_non_footprint_metric_remains_eligible(tmp_path, monkeypatch):
    from mdscope.reporting import plotting

    monkeypatch.setattr(plotting, "collect_time_series", lambda *a, **k: {
        (None, "coordination_count"): pd.DataFrame(
            {"time_ps": [0., 100.], "value": [0., 0.]})
    })
    cfg = {"plots": {"enabled": True, "distributions": False,
                     "autocorrelation": False, "bootstrap_tables": False},
           "systems": {"control": {"replicates": [{"label": "rep1"}]}}}
    plotting.plot_time_series(cfg, tmp_path,
                              {"control": [tmp_path / "control/rep1"]})
    assert (tmp_path / "project_summary/figures/unmapped/secondary/time_series@coordination_count@timeseries.png").is_file()


def test_zero_footprint_filter_preserves_missing_vs_zero_semantics():
    from mdscope.reporting.plotting import _all_zero_lipid_footprint

    metric = "lipid_footprint_BGAL_lower_5_10A_fold_available_fraction"
    zero = pd.DataFrame({"value": [0., 0.]})
    missing = pd.DataFrame({"value": [np.nan]})
    positive = pd.DataFrame({"value": [0., 0.002]})
    assert _all_zero_lipid_footprint(metric, [zero, missing])
    assert not _all_zero_lipid_footprint(metric, [positive, zero])
    assert not _all_zero_lipid_footprint(metric, [missing])
    assert not _all_zero_lipid_footprint("rmsd_A", [zero])
