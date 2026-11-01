"""Independent tiny-system references for solvation counting, PBC and kinetics."""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.observables.solvation import (SolvationObserver, contact_episodes,
                                           rdf_from_counts)


def _universe():
    u = mda.Universe.empty(5, n_residues=4, atom_resindex=[0, 0, 1, 2, 3],
                           trajectory=True)
    u.add_TopologyAttr("names", ["OE1", "OE2", "OW", "OW", "SOD"])
    u.add_TopologyAttr("resnames", ["GLU", "SOL", "SOL", "SOD"])
    u.add_TopologyAttr("resids", [7, 1, 2, 3])
    xyz = np.array([[1, 1, 1], [1.5, 1, 1], [3, 1, 1],
                    [19.5, 1, 1], [3, 1, 2]], np.float32)
    shifted = xyz.copy()
    shifted[2] = [12, 1, 1]
    u.load_new(np.array([xyz, shifted]), format=MemoryReader,
               dimensions=np.tile([20, 20, 20, 90, 90, 90], (2, 1)))
    return u


def test_bulk_rdf_formula():
    edges = np.array([0., 1., 2.])
    expected = np.array([3., 21.])
    np.testing.assert_allclose(rdf_from_counts(expected, edges, 3 / (4*np.pi/3)),
                               [1, 1])


def test_contact_episodes_censor_and_gaps():
    events = [(0, "water", "glu", 4), (1, "water", "glu", 4),
              (3, "water", "glu", 4)]
    out = contact_episodes(events, [(i, i*5.) for i in range(4)])
    assert len(out) == 2
    assert out.duration_ps.tolist() == [5., 0.]
    assert out.left_censored.tolist() == [True, False]
    assert out.right_censored.tolist() == [False, True]


def test_periodic_unique_water_count_and_restart_histograms(tmp_path):
    cfg = {"analyses": {"solvation": {
        "enabled": True, "preset": "kinetics",
        "targets": {"protein": {"carboxylate": "resname GLU"}},
        "water_selection": "resname SOL and name OW",
        "ion_selection": "resname SOD",
        "search_radius_A": 5., "rdf_bin_width_A": 0.5,
        "shell_cutoffs_A": {"water": 3.5, "ions": 4.}
    }}}
    u = _universe()
    obs = SolvationObserver(u, cfg, tmp_path)
    obs.update(u.trajectory[0], 0, 0.)
    obs.flush(0)
    obs.update(u.trajectory[1], 1, 5.)
    obs.flush(1)
    obs.finalize()
    root = tmp_path / "extensions" / "solvation"
    counts = pd.read_parquet(root / "counts.parquet")
    water = counts[counts.solvent == "water"]
    assert water.n_solvent.tolist() == [2, 1]
    ions = counts[counts.solvent == "ions"]
    assert ions.n_solvent.tolist() == [1, 1]
    assert len(pd.read_parquet(root / "rdf.parquet")) == 20
    events = pd.read_parquet(root / "events.parquet")
    assert len(events[events.solvent == "water"]) == 3
    assert (root / "episodes.parquet").exists()
    assert (root / "metadata.json").exists()


def test_local_grid_uses_same_water_candidates(tmp_path):
    u = _universe()
    cfg = {"analyses": {"solvation": {
        "enabled": True, "preset": "standard",
        "targets": {"protein": {"oxygen": "resname GLU and name OE1 OE2"}},
        "water_selection": "resname SOL and name OW",
        "grid_frames": {"protein/oxygen/0": "index 0 1 4"},
        "grid_width_A": 10., "grid_spacing_A": 1.,
        "search_radius_A": 8.,
    }}}
    observer = SolvationObserver(u, cfg, tmp_path)
    observer.update(u.trajectory[0], 0, 0.)
    observer.flush(0)
    observer.finalize()
    density = pd.read_parquet(tmp_path/"extensions/solvation/water_density_3d.parquet")
    assert density["count"].sum() == 2
    assert np.all(np.isfinite(density["g_bulk"]))


def test_water_hbond_kernel_checks_angle_and_periodic_geometry():
    from mdscope.observables.solvation import water_donor_hbonds
    box = np.array([30., 30., 30., 90., 90., 90.])
    sites = np.array([[-2.7, 0., 0.], [2.7, 0., 0.],
                      [0., 2.7, 0.]])
    waters = np.array([[0., 0., 0.]])
    hydrogens = np.array([[[-1., 0., 0.], [1., 0., 0.]]])
    candidate = np.array([[0,0],[1,0],[2,0]])
    found = water_donor_hbonds(sites, waters, hydrogens, candidate, box)
    assert set(map(tuple, found)) == {(0,0,0), (0,1,1)}
    shifted = sites.copy()
    shifted[0, 0] += 30.
    wrapped = water_donor_hbonds(shifted, waters, hydrogens, candidate, box)
    assert set(map(tuple, wrapped)) == set(map(tuple, found))
    assert not len(water_donor_hbonds(sites, waters, hydrogens,
                                     candidate, box, min_angle_deg=181.))


def test_water_donor_bridge_persisted_with_topology_bonds(tmp_path):
    u = mda.Universe.empty(5, n_residues=3, atom_resindex=[0, 1, 2, 2, 2],
                           trajectory=True)
    u.add_TopologyAttr("names", ["OD1", "OE1", "OW", "H1", "H2"])
    u.add_TopologyAttr("resnames", ["ASP", "GLU", "SOL"])
    u.add_TopologyAttr("resids", [2, 5, 9])
    u.add_TopologyAttr("bonds", [(2, 3), (2, 4)])
    xyz = np.array([[-2.7, 0, 0], [2.7, 0, 0], [0, 0, 0],
                    [-1, 0, 0], [1, 0, 0]], dtype=np.float32) + [10, 10, 10]
    u.load_new(xyz[None, :, :], format=MemoryReader,
               dimensions=np.array([[30, 30, 30, 90, 90, 90]]))
    cfg = {"analyses": {"solvation": {
        "enabled": True, "preset": "detailed",
        "targets": {"protein": {"acceptor": "resname ASP GLU and name OD1 OE1"}},
        "water_selection": "resname SOL and name OW",
        "search_radius_A": 5.,
        "options": {"water_hbond_bridges": True},
    }}}
    observer = SolvationObserver(u, cfg, tmp_path)
    observer.update(u.trajectory[0], 0, 0.)
    observer.flush(0)
    observer.finalize()
    bridge = pd.read_parquet(
        tmp_path/"extensions/solvation/water_hbond_bridges.parquet")
    assert len(bridge) == 1
    assert {bridge.target_a.iloc[0], bridge.target_b.iloc[0]} == {
        "protein/acceptor/0", "protein/acceptor/1"}


def test_explicit_ion_species_labels_and_nonoverlap(tmp_path):
    u = _universe()
    cfg = {"analyses": {"solvation": {
        "enabled": True, "preset": "standard",
        "targets": {"protein": {"acidic": "resname GLU and name OE1 OE2"}},
        "ion_selections": {"Na": "resname SOD"},
        "search_radius_A": 5., "shell_cutoffs_A": {"ions": 4.},
    }}}
    obs = SolvationObserver(u, cfg, tmp_path)
    obs.update(u.trajectory[0], 0, 0.)
    obs.flush(0)
    obs.finalize()
    counts = pd.read_parquet(tmp_path/"extensions/solvation/counts.parquet")
    assert counts.solvent.tolist() == ["ions/Na"]
    assert counts.n_solvent.tolist() == [1]
    rdf = pd.read_parquet(tmp_path/"extensions/solvation/rdf.parquet")
    assert set(rdf.solvent) == {"ions/Na"}
    cfg["analyses"]["solvation"]["ion_selections"] = {
        "Na": "resname SOD", "Duplicate": "resname SOD"}
    import pytest
    with pytest.raises(ValueError, match="overlap"):
        SolvationObserver(u, cfg, tmp_path/"bad")


def test_mixed_ion_selection_configuration_rejected(tmp_path):
    import pytest
    u = _universe()
    cfg = {"analyses": {"solvation": {
        "enabled": True,
        "targets": {"protein": {"acidic": "resname GLU"}},
        "ion_selection": "resname SOD", "ion_selections": {"Na":"resname SOD"}
    }}}
    with pytest.raises(ValueError, match="OR"):
        SolvationObserver(u, cfg, tmp_path)


def test_hs_chain_grouping_reuses_covalent_topology_not_chain_ids():
    from mdscope.observables.solvation import _targets
    u = mda.Universe.empty(
        8, n_residues=6,
        atom_resindex=[0, 1, 2, 2, 3, 4, 4, 5], trajectory=True)
    u.add_TopologyAttr("names", ["OG", "OG", "O1", "C1", "O1", "O1", "C1", "O1"])
    u.add_TopologyAttr("resnames",
        ["SER", "SER", "AGLCNA", "BGLCA", "AGLCNA", "BGLCA"])
    u.add_TopologyAttr("resids", [39, 61, 1, 2, 3, 4])
    u.add_TopologyAttr("chainIDs", ["PROA"] * 8)
    u.add_TopologyAttr("bonds", [(0,2), (2,3), (3,4), (1,5), (5,6), (6,7)])
    target_spec = {
        "targets": {"hs": {"oxygen": {
            "selection": "resname AGLCNA BGLCA and name O1",
            "group_by": "hs_chain"}}},
        "hs_auto": {
            "glycan_selection": "resname AGLCNA BGLCA",
            "attachment_selection": "resname SER and name OG",
            "min_glycan_residues": 2, "expected_chains": 2}}
    grouped = _targets(u, target_spec)
    assert len(grouped) == 2
    assert {t[0] for t in grouped} == {
        "hs/oxygen/PROA_S39", "hs/oxygen/PROA_S61"}
    assert sorted(len(t[3]) for t in grouped) == [2,2]


def test_tracking_class_filter_preserves_counts_but_limits_events(tmp_path):
    u = _universe()
    config = {"analyses": {"solvation": {
        "enabled": True, "preset": "kinetics",
        "targets": {
            "protein": {"acceptor": "resname GLU and name OE1"},
            "lipid": {"other": "resname GLU and name OE2"}},
        "water_selection": "resname SOL and name OW",
        "search_radius_A": 5.,
        "options": {"track_target_classes": ["protein"]},
    }}}
    observer = SolvationObserver(u, config, tmp_path)
    observer.update(u.trajectory[0], 0, 0.)
    observer.flush(0)
    observer.finalize()
    result = tmp_path/"extensions/solvation"
    counts = pd.read_parquet(result/"counts.parquet")
    events = pd.read_parquet(result/"events.parquet")
    assert set(counts.target_class) == {"protein", "lipid"}
    assert set(events.target_id) == {"protein/acceptor/0"}
    assert pd.read_parquet(result/"bridges.parquet").empty if (
        result/"bridges.parquet").exists() else True


def test_fused_solvation_restart_rebuilds_same_outputs(tmp_path):
    from mdscope.pipeline.fused_extensions import FusedExtensions
    cfg = {"analyses": {"solvation": {
        "enabled": True, "preset": "kinetics",
        "targets": {"protein": {"acidic": "resname GLU and name OE1 OE2"}},
        "water_selection": "resname SOL and name OW",
        "ion_selections": {"Na": "resname SOD"},
        "search_radius_A": 5., "rdf_bin_width_A": .5,
        "shell_cutoffs_A": {"water": 3.5, "ions": 4.},
    }}, "time": {"dt_ps": 5., "origin_ps": 0.}}
    paths = []
    for restart in (False, True):
        root = tmp_path/("resumed" if restart else "clean")
        u = _universe()
        obs = FusedExtensions(u, cfg, root)
        obs.prune_from(0)
        obs.update(u.trajectory[0], 0)
        obs.flush(0)
        if restart:
            # Reconstruct observer at chunk boundary. Solvation histograms,
            # events and frame mapping must be rebuilt from saved chunks.
            u = _universe()
            obs = FusedExtensions(u, cfg, root)
            obs.prune_from(1)
        obs.update(u.trajectory[1], 1)
        obs.flush(1)
        obs.finalize()
        paths.append(root/"extensions/solvation")
    for name in ("counts.parquet", "rdf.parquet", "events.parquet",
                 "episodes.parquet", "survival.parquet",
                 "species_timeseries.parquet", "target_summary.parquet"):
        pd.testing.assert_frame_equal(
            pd.read_parquet(paths[0]/name),
            pd.read_parquet(paths[1]/name), check_dtype=False)


def test_target_water_overlap_is_rejected(tmp_path):
    import pytest
    u = _universe()
    cfg = {"analyses": {"solvation": {
        "enabled": True,
        "targets": {"water": {"self": "resname SOL and name OW"}},
        "water_selection": "resname SOL and name OW"}}}
    with pytest.raises(ValueError, match="overlap"):
        SolvationObserver(u, cfg, tmp_path)


def test_unique_coordination_matches_python_reference_with_duplicates():
    from collections import defaultdict
    from mdscope.observables.solvation import unique_target_solvent_pairs
    targets = np.array([2, 0, 2, 2, 1, 0, 2, 1, 0], dtype=np.int32)
    water = np.array([7, 3, 7, 9, 1, 3, 8, 1, 4], dtype=np.int32)
    reference = defaultdict(set)
    for i, j in zip(targets, water):
        reference[int(i)].add(int(j))
    counts, pairs = unique_target_solvent_pairs(targets, water, 5)
    assert counts.tolist() == [2, 1, 3, 0, 0]
    assert set(map(tuple, pairs)) == {
        (i, j) for i, values in reference.items() for j in values}
    empty_counts, empty_pairs = unique_target_solvent_pairs(
        np.empty(0, int), np.empty(0, int), 5)
    assert empty_counts.tolist() == [0]*5
    assert empty_pairs.shape == (0, 2)


def test_neighbor_methods_produce_identical_coordination(tmp_path):
    config = {"analyses": {"solvation": {
        "enabled": True, "preset": "standard",
        "targets": {"protein": {"acid": "resname GLU and name OE1 OE2"}},
        "water_selection": "resname SOL and name OW",
        "search_radius_A": 5.,
    }}}
    results = []
    for backend in ("auto", "nsgrid", "pkdtree", "bruteforce"):
        config["analyses"]["solvation"]["neighbor_method"] = backend
        u = _universe()
        observer = SolvationObserver(u, config, tmp_path/backend)
        ts = u.trajectory[0]
        # The shared fixture contains an EXACT RDF bin-edge separation
        # (2.0 Å). Distinct native PBC backends may round that distance
        # to opposite sides of the edge without changing coordination.
        # Move only this test's solvent slightly off that discontinuity.
        u.atoms[3].position = [19.43, 1., 1.]
        observer.update(ts, 0, 0.)
        observer.flush(0)
        observer.finalize()
        root = tmp_path/backend/"extensions/solvation"
        results.append((pd.read_parquet(root/"counts.parquet"),
                        pd.read_parquet(root/"rdf.parquet")))
    for counts, rdf in results[1:]:
        pd.testing.assert_frame_equal(counts, results[0][0])
        pd.testing.assert_frame_equal(rdf, results[0][1])
