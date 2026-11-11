import numpy as np
import pandas as pd

from mdscope.statistics import analyze_project


def test_statistics_exports_intervals_and_explicit_coverage(tmp_path):
    root = tmp_path / "results"
    cfg = {
        "analyses": {"membrane": {"msd": {"enabled": True, "multi_origin": True}}},
        "statistics": {"bootstrap_draws": 64, "block": {"length": 2}},
    }
    repouts = {}
    for system, offset in (("control", 0), ("treatment", 1)):
        repouts[system] = []
        for index in range(2):
            rep = root / system / f"rep{index + 1}"
            repouts[system].append(rep)
            summary = rep / "summary"
            summary.mkdir(parents=True)
            time = np.arange(40, dtype=float)
            signal = np.full(len(time), 5.0 + offset)
            signal[::4] = np.nan
            frame = pd.DataFrame({
                "time_ps": time,
                "sdc4_local_upper_height_A": signal,
                "sdc4_local_thickness_A": np.where(time == 1, 3.0 + offset, np.nan),
            })
            if offset:
                frame["treatment_only"] = 1.0
            frame.to_parquet(summary / "timeseries.parquet", index=False)
            lag = np.arange(40 + index * 10, dtype=float)
            pd.DataFrame({
                "metric": "lateral", "lag_ps": lag,
                "msd_A2": (.4 + .4 * offset) * lag,
                "n_time_origins": len(lag) - lag.astype(int),
            }).to_parquet(summary / "multi_origin_msd.parquet", index=False)
            rdf = rep / "extensions" / "solvation"
            rdf.mkdir(parents=True)
            pd.DataFrame({
                "solvent": "water", "target_class": "HS", "site_class": "oxygen",
                "r_inner_A": [0.0, 2.0], "r_outer_A": [2.0, 4.0],
                "g_bulk": [1.0, 3.0 + offset],
                "g_aqueous_reference": [1.0, np.nan if index else 3.0 + offset],
            }).to_parquet(rdf / "rdf.parquet", index=False)

    analyze_project(cfg, root, repouts)
    stats = root / "project_summary" / "statistics"
    coverage_path = stats / "coverage" / "contrast_coverage.parquet"
    coverage = pd.read_parquet(coverage_path)
    geometry = coverage.loc[
        coverage.observable.str.contains("sdc4_local_upper_height_A", na=False)
        & ~coverage.observable.str.contains("_available_fraction", na=False)
    ]
    assert set(geometry.status) == {"above_zero"}
    np.testing.assert_allclose(geometry.ci_low, 1)
    assert set(coverage.loc[coverage.observable.str.contains("local_thickness_A\"", na=False), "status"]) == {
        "unavailable_temporal_support"
    }
    assert set(coverage.loc[coverage.observable.str.contains("treatment_only", na=False), "status"]) == {"one_system_only"}
    assert set(coverage.loc[coverage.observable.str.contains("available_fraction", na=False), "status"]) == {
        "includes_zero"
    }

    diffusion = pd.read_parquet(stats / "diffusion" / "common_lag_comparisons.parquet")
    np.testing.assert_allclose(diffusion.difference_ci_low, 1)
    np.testing.assert_allclose(diffusion.difference_ci_high, 1)
    fits = pd.read_parquet(stats / "diffusion" / "replicate_common_lag_slopes.parquet")
    assert fits.fit_start_ps.nunique() == fits.fit_end_ps.nunique() == 1
    assert fits.n_frames.nunique() == 2

    shell = pd.read_parquet(stats / "solvation" / "rdf_shell_comparisons.parquet").set_index("normalization")
    expected = (3.5 ** 3 - 2 ** 3) / 3.5 ** 3
    np.testing.assert_allclose(shell.loc["g_bulk", ["difference_ci_low", "difference_ci_high"]].astype(float), expected)
    assert shell.loc["g_aqueous_reference", "ci_status"] == "incomplete_replica_coverage"
    assert pd.isna(shell.loc["g_aqueous_reference", "difference_ci_low"])

    analyze_project(cfg, root, repouts)
    pd.testing.assert_frame_equal(coverage, pd.read_parquet(coverage_path))
