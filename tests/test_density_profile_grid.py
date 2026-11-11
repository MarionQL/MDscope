import numpy as np
import pandas as pd

from mdscope.statistics.profiles import analyze_profiles


def test_density_profiles_use_shared_physical_bins(tmp_path):
    cfg = {
        "analyses": {"density_1d": {"enabled": True, "bin_width_A": 1.0}},
        "statistics": {
            "bootstrap_draws": 32,
            "block": {"length": 2},
        },
    }
    # One source grid straddles the shared bin edges. A 50/50 split of
    # neighboring slabs gives a different result than matching their labels.
    axes = (
        np.arange(-2.5, 3.0),
        np.arange(-2.0, 3.0),
    )
    reps = {}
    for system, shift in (("first", 0.0), ("second", 2.0)):
        reps[system] = []
        for index, centers in enumerate(axes):
            rep = tmp_path / system / f"rep{index + 1}"
            frames = rep / "profiles" / "density_1d" / "default" / "frames"
            frames.mkdir(parents=True)
            baseline = np.arange(len(centers), dtype=float) + 1
            if index == 1:
                baseline += 0.5
            np.savez(
                frames / "water__part-000000.npz",
                positions_A=centers,
                time_ps=np.arange(8, dtype=float),
                density=np.tile(baseline + shift, (8, 1)),
            )
            reps[system].append(rep)

    analyze_profiles(cfg, tmp_path, reps)
    out = tmp_path / "project_summary" / "statistics" / "profiles"
    comparison_path = out / "system_comparisons.parquet"
    assert comparison_path.is_file(), "Offset density grids lost the system comparison"
    comp = pd.read_parquet(comparison_path)
    np.testing.assert_allclose(comp.position, [-1.5, -0.5, 0.5, 1.5])
    np.testing.assert_allclose(comp.difference_system2_minus_system1, 2.0)
    np.testing.assert_allclose(comp.difference_ci_low, 2.0)
    np.testing.assert_allclose(comp.difference_ci_high, 2.0)
    profiles = pd.read_parquet(out / "system_profiles.parquet")
    reference = profiles.loc[profiles.system.eq("first"), "mean"]
    np.testing.assert_allclose(reference, [2.0, 3.0, 4.0, 5.0])
    # The independent overlap integral over [-2, 2] is 14 for both grids.
    per_replica = pd.read_parquet(out / "replicate_profiles.parquet")
    first = per_replica.loc[per_replica.system.eq("first")]
    np.testing.assert_allclose(first.groupby("replicate").value.sum(), [14.0, 14.0])
