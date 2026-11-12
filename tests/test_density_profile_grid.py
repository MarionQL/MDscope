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


def test_midplane_density_is_invariant_to_drift_and_periodic_wrap(tmp_path):
    import MDAnalysis as mda

    from mdscope.observables.density_1d import DensityProfileAccumulator, temporal_density_profiles

    u = mda.Universe.empty(7, n_residues=7, atom_resindex=np.arange(7), trajectory=True)
    u.add_TopologyAttr("names", ["P"] * 5 + ["O"] * 2)
    spec = {
        "enabled": True, "axis": "z", "bin_width_A": 1.0, "half_width_A": 20.0,
        "membrane_center_selection": "name P", "selections": {"probe": "name O"},
    }
    cfg = {"analyses": {"density_1d": {"enabled": True, "runs": {"membrane": spec}}}}
    # Unequal leaflet populations have the same midplane. Probe positions
    # land in the bins centered at -3.5 and +5.5 Å.
    relative = np.array([-10., -10., 10., 10., 10., -3.4, 5.4])
    u.dimensions = [20., 30., 100., 90., 90., 90.]
    observer = DensityProfileAccumulator(u, cfg, spec, "membrane")
    for i, origin in enumerate((45., 96.)):
        u.atoms.positions = np.column_stack([np.zeros(7), np.zeros(7), (origin + relative) % 100.])
        observer.update(u, float(i))
    observer.flush_frames(tmp_path, 0)
    observer.finalize(tmp_path)
    profile = temporal_density_profiles(tmp_path, cfg)["density_1d__membrane__probe"]
    np.testing.assert_array_equal(profile["values"][0], profile["values"][1])
    observed = profile["values"][0]
    np.testing.assert_array_equal(profile["positions"][observed > 0], [-3.5, 5.5])
    np.testing.assert_allclose(observed.sum() * 20 * 30, 2.)
