import numpy as np
import pandas as pd

from mdscope.statistics import analyze_project


def test_coupling_uses_shared_grid_without_compressing_gaps(tmp_path):
    cfg = {
        "analyses": {"coupling_screen": {
            "enabled": True, "include_defaults": False, "block_frames": 2,
            "pairs": [{"name": "contact_vs_shape", "x": "contact", "y": "shape"}],
        }},
        "statistics": {"bootstrap_draws": 32, "scalar_window_ps": [10, 49]},
    }
    repouts = {}
    for system, sign in (("control", -1), ("treatment", 1)):
        repouts[system] = []
        for i in range(2):
            rep = tmp_path / system / f"rep{i + 1}"
            summary = rep / "summary"
            summary.mkdir(parents=True)
            repouts[system].append(rep)
            time = np.arange(10 * i, 50 + 10 * i, dtype=float)
            x = time.copy()
            x[time == 17] = np.nan
            pd.DataFrame({
                "time_ps": time, "contact": x, "shape": sign * time,
                "window_check": np.where((time >= 10) & (time <= 49), 1.0, 100.0),
            }).to_parquet(summary / "timeseries.parquet", index=False)

    analyze_project(cfg, tmp_path, repouts)
    out = tmp_path / "project_summary" / "statistics"
    coeffs = pd.read_parquet(out / "coupling" / "replicate_correlations.parquet")
    assert coeffs.n_blocks.tolist() == [19, 19, 19, 19]
    assert coeffs.window_start_ps.tolist() == [10] * 4
    assert coeffs.window_end_ps.tolist() == [49] * 4
    contrasts = pd.read_parquet(out / "coupling" / "system_comparisons.parquet")
    np.testing.assert_allclose(contrasts[["difference_ci_low", "difference_ci_high"]], 2)
    coverage = pd.read_parquet(out / "coverage" / "contrast_coverage.parquet")
    assert coverage.loc[coverage.family.eq("coupling"), "status"].item() == "above_zero"
    scalar = pd.read_parquet(out / "time_series" / "replicate_statistics.parquet")
    np.testing.assert_allclose(scalar.loc[scalar.metric.eq("window_check"), "mean"], 1)
