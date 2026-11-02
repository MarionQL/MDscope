"""Synthetic PBC, output, restart and shared-state regression tests for SDC4."""
from pathlib import Path
from types import SimpleNamespace
import json

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.mechanisms.sdc4_mechanisms import SDC4Mechanisms, _contacts, _angle, _backbone_axis, _local_surface, _normal_order_setup, _order_by_residue
from mdscope.mechanisms.shared_states import fit_shared_states


def _synthetic():
    names, resnames, resids, chainids, resindices, coords = [], [], [], [], [], []
    def residue(name, resid, chain, atoms):
        ix = len(resnames)
        resnames.append(name)
        resids.append(resid)
        for atom, xyz in atoms:
            names.append(atom)
            chainids.append(chain)
            resindices.append(ix)
            coords.append(xyz)
    for chain, x in (("PROA", 0.), ("PROB", 8.)):
        for i, aa in enumerate(("ASN", "ALA", "ILE", "PRO", "ALA",
                                "GLY", "GLY", "ILE", "VAL", "GLY", "ALA", "ALA"), 1):
            z = 39. - 4*i if i <= 5 else 12. - 4*(i-6)
            residue(aa, i, chain, [("CA", [x, 0., z])])
    # 20 single-residue lipid anchors, balanced between upper/lower.
    for i in range(20):
        lipid = ("POPC", "GM1a", "BGAL", "CHL1")[i % 4]
        atom = "P" if lipid == "POPC" else "O3" if lipid == "CHL1" else "NF"
        x = ((i % 5)-2)*9.
        y = ((i // 5) % 2)*12. - 6.
        z = 17. if i < 10 else -17.
        residue(lipid, i+1, lipid, [(atom, [x, y, z])])
    for chain, x in (("HS_A", 0.), ("HS_B", 8.)):
        residue("AGLCNA", 1, chain, [("C1", [x, 1., 26.])])
        residue("AIDOA", 2, chain, [("C1", [x, 1., 31.])])
    for i in range(6):
        residue("TIP3", i+1, "SOL", [("OH2", [i*3., 1., 30.])])
    residue("SOD", 1, "ION", [("SOD", [1., 1., 29.])])
    u = mda.Universe.empty(
        len(names), n_residues=len(resnames),
        atom_resindex=resindices, trajectory=True,
    )
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("resids", resids)
    u.add_TopologyAttr("chainIDs", chainids)
    initial = np.asarray(coords, np.float32)
    frames = np.stack([initial.copy() for _ in range(4)])
    for i in range(4):
        frames[i, :24, 1] += .1 * i
        frames[i, -7:, 1] += .05 * i
    u.load_new(frames, format=MemoryReader,
               dimensions=np.tile([100., 100., 100., 90., 90., 90.], (4, 1)))
    return u


def _cfg():
    return {
        "selections": {"protein": "protein"},
        "trajectory": {"start": 0, "stop": 4, "stride": 1, "chunk_frames": 2},
        "time": {"dt_ps": 10., "origin_ps": 0.},
        "analyses": {"sdc4_mechanisms": {
            "enabled": True, "axis": "z",
            "tm_chains": {
                "A": "protein and chainID PROA and resid 6:12 and name CA",
                "B": "protein and chainID PROB and resid 6:12 and name CA",
            },
            "ectodomains": {
                "A": "protein and chainID PROA and resid 1:5 and name CA",
                "B": "protein and chainID PROB and resid 1:5 and name CA",
            },
            "hs_chains": {
                "A1": "chainID HS_A", "B1": "chainID HS_B",
            },
            "hs_attachment_sites": {
                "A1": "chainID PROA and resid 4 and name CA",
                "B1": "chainID PROB and resid 4 and name CA",
            },
            "lipid_anchors": {
                "POPC": "resname POPC and name P",
                "GM1a": "resname GM1a and name NF",
                "BGAL": "resname BGAL and name NF",
                "CHL1": "resname CHL1 and name O3",
            },
            "lipid_contact_selection": "resname POPC GM1a BGAL CHL1",
            "surface_anchor_species": ["POPC", "GM1a", "BGAL", "CHL1"],
            "water_oxygen_selection": "resname TIP3 and name OH2",
            "ion_selections": {"Na": "resname SOD"},
            "motif_selection": "protein and name CA",
            "motifs": {"NxIP": "N.IP"},
            "require_motifs": True, "require_gxxxg": True,
            "local_order": {"enabled": False},
            "local_radius_A": 30., "local_min_lipids": 8,
            "map_extent_A": 35., "map_bin_A": 5.,
            "contact_cutoff_A": 5., "hydration_cutoff_A": 3.5,
            "ion_cutoff_A": 4.,
        }},
    }


def _run(out, resume=False):
    u, cfg = _synthetic(), _cfg()
    obs = SDC4Mechanisms(u, cfg, out)
    if resume:
        for i in range(3):
            obs.update(u.trajectory[i], i)
            if i == 1:
                obs.flush(0)
        obs.flush(1)
        u, cfg = _synthetic(), _cfg()
        obs = SDC4Mechanisms(u, cfg, out)
        obs.prune_from(1)
        start, part = 2, 1
    else:
        obs.prune_from(0)
        start, part = 0, 0
    for i in range(start, 4):
        obs.update(u.trajectory[i], i)
        if (i+1) % 2 == 0:
            obs.flush(part)
            part += 1
    obs.finalize()
    return cfg


def test_sdc4_observer_and_checkpoint_recovery(tmp_path):
    _run(tmp_path / "clean")
    _run(tmp_path / "resume", resume=True)
    for folder in ("clean", "resume"):
        out = tmp_path / folder / "extensions/sdc4_mechanisms"
        data = pd.read_parquet(out / "timeseries.parquet")
        assert len(data) == 4
        assert np.all(np.isfinite(data["sdc4_tm_crossing_angle_deg"]))
        assert np.all(np.isfinite(data["sdc4_ectodomain_A_height_A"]))
        assert "sdc4_motif_NxIP_PROA_1_steric_occlusion_fraction" in data
        assert "sdc4_hs_A1_hydration_water_count" in data
        assert np.all(np.isfinite(data["sdc4_local_thickness_A"]))
        assert "sdc4_ectodomain_A_membrane_contact_state" in data
        with np.load(out / "lipid_density_maps.npz", allow_pickle=False) as maps:
            assert int(maps["nframes"]) == 4
            assert int(maps["counts"].sum()) <= 4*20
            assert maps["density_per_nm2"].shape[0] == 8
            assert np.all(np.isfinite(maps["bulk_density_per_nm2"]))
        audit = json.loads((out / "selection_audit.json").read_text())
        assert len(audit["motif_resindices"]) == 2
        assert all(len(audit["gxxxg_resindices"][chain]) == 5 for chain in ("A", "B"))
        assert audit["gxxxg_tm_sequences"] == {"A": "GGIVGAA", "B": "GGIVGAA"}
    a = tmp_path / "clean/extensions/sdc4_mechanisms"
    b = tmp_path / "resume/extensions/sdc4_mechanisms"
    pd.testing.assert_frame_equal(pd.read_parquet(a / "timeseries.parquet"),
                                  pd.read_parquet(b / "timeseries.parquet"))
    with np.load(a / "lipid_density_maps.npz") as left, np.load(b / "lipid_density_maps.npz") as right:
        np.testing.assert_array_equal(left["counts"], right["counts"])


def test_missing_ggivg_fails_with_tm_sequence(tmp_path):
    u, cfg = _synthetic(), _cfg()
    # Changing one residue in B must fail, even if a different GxxxG
    # pattern is present elsewhere in the selected TM segment.
    b = u.select_atoms("chainID PROB and resid 7 and name CA")
    assert len(b) == 1
    u.residues[b.resindices[0]].resname = "ALA"
    with pytest.raises(ValueError, match="No GGIVG found in TM chain B.*sequence="):
        SDC4Mechanisms(u, cfg, tmp_path)


def test_local_surface_recovers_quadratic_curvature():
    xy = np.array([[x, y] for x in range(-3, 4) for y in range(-3, 4)], float)
    z = .01*xy[:, 0]**2 + .02*xy[:, 1]**2 + .3*xy[:, 0]
    assert _local_surface(xy, z) == pytest.approx(.03, abs=1e-9)


def test_shared_state_model_replicate_balance_and_sampling(tmp_path):
    cfg = {"analyses": {
        "sdc4_mechanisms": {"enabled": True},
        "shared_states": {"enabled": True,
                          "features": ["sdc4_ectodomain_A_height_A",
                                       "sdc4_tm_crossing_angle_deg"],
                          "n_states": 2, "seed": 7,
                          "bootstrap_draws": 30,
                          "training_frames_per_replicate": 30,
                          "min_complete_fraction": 1.0},
    }}
    paths = {}
    for system, shift in (("carb", 0.), ("nocarb", 3.)):
        paths[system] = []
        for rep in ("rep1", "rep2"):
            folder = tmp_path / system / rep
            out = folder / "extensions/sdc4_mechanisms"
            out.mkdir(parents=True)
            x = np.r_[np.linspace(1, 2, 20), np.linspace(8, 9, 20)]
            df = pd.DataFrame({
                "time_ps": np.arange(40)*100.,
                "sdc4_ectodomain_A_height_A": x + shift,
                "sdc4_tm_crossing_angle_deg": x * 2 + shift,
                "sdc4_hs_A1_membrane_contact_state": np.r_[[0]*20, [1]*20],
                "sdc4_hs_A1_hydration_water_per_heavy_atom": x,
                "sdc4_hs_A1_ion_Na_count": x*0 + 2,
            })
            df.to_parquet(out / "timeseries.parquet", index=False)
            paths[system].append(folder)
    folder = fit_shared_states(cfg, tmp_path, paths)
    pop = pd.read_parquet(folder / "replicate_populations.parquet")
    assert len(pop) == 8
    assert np.allclose(pop.groupby(["system", "replicate"]).fraction.sum(), 1)
    assignments = pd.read_parquet(folder / "frame_assignments.parquet")
    assert len(assignments) == 160
    conditioned = pd.read_parquet(folder / "state_conditioned_hs_solvation.parquet")
    assert "hydration_water_per_heavy_atom_mean" in conditioned
    assert len(pd.read_parquet(folder / "sampling_diagnostics.parquet")) == 4
    model = json.loads((folder / "model.json").read_text())
    assert model["n_replicates"] == 4


def test_local_lipid_order_matches_minus_sch_convention():
    u = mda.Universe.empty(
        8, n_residues=2,
        atom_resindex=[0, 0, 0, 0, 1, 1, 1, 1],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["C21", "H21", "C22", "H22"] * 2)
    u.add_TopologyAttr("resnames", ["POPC", "POPC"])
    u.add_TopologyAttr("resids", [1, 2])
    u.add_TopologyAttr("bonds", [(0, 1), (2, 3), (4, 5), (6, 7)])
    xyz = np.asarray([
        [1, 1, 1], [1, 1, 2], [3, 1, 1], [3, 1, 2],
        [8, 1, 1], [9, 1, 1], [8, 3, 1], [9, 3, 1],
    ], np.float32)
    u.load_new(xyz[None], format=MemoryReader,
               dimensions=np.array([[100, 100, 100, 90, 90, 90]]))
    setup = _normal_order_setup(u, {
        "enabled": True, "lipid_selection": "resname POPC",
        "min_carbons_per_tail": 2,
        "tails": {"sn1": {"carbon_pattern": r"^C2(\d+)$"}},
    })
    ids, values = _order_by_residue(
        setup, np.array([100., 100., 100., 90., 90., 90.]), 2
    )
    assert len(ids) == 2
    np.testing.assert_allclose(values, [-1.0, 0.5], atol=1e-8)


def test_tm_contacts_independent_distance_oracle_and_periodic_boundary():
    """Heavy-atom contact pair indices must match brute-force minimum image."""
    u = mda.Universe.empty(6, n_residues=4,
                           atom_resindex=[0, 0, 1, 2, 2, 3],
                           trajectory=True)
    u.add_TopologyAttr("names", ["CA", "CB", "CA", "CA", "CB", "CA"])
    u.add_TopologyAttr("resnames", ["GLY", "ILE", "GLY", "VAL"])
    u.add_TopologyAttr("resids", [1, 2, 10, 11])
    u.atoms.positions = np.asarray([
        [1.0, 1, 1], [1.2, 1, 1], [6, 1, 1],
        [99, 1, 1], [98.5, 1, 1], [50, 1, 1],
    ], dtype=np.float32)
    box = np.array([100., 100., 100., 90., 90., 90.])
    left, right = u.atoms[:3], u.atoms[3:]
    observed = set(map(tuple, _contacts(left, right, 3., box)))
    expected = set()
    for i, a in enumerate(left.positions):
        for j, b in enumerate(right.positions):
            delta = a - b
            delta -= np.rint(delta / box[:3]) * box[:3]
            if np.linalg.norm(delta) <= 3.:
                expected.add((i, j))
    assert observed == expected == {(0, 0), (0, 1), (1, 0), (1, 1)}
    residue_pairs = {(int(left.resindices[i]), int(right.resindices[j]))
                     for i, j in observed}
    assert len(residue_pairs) == 1
    assert len(_contacts(left, right, 1.0, box)) == 0
    right.positions += [0, 20, 0]
    assert len(_contacts(left, right, 3., box)) == 0


def test_tm_axis_orientation_and_crossing_angle_oracle():
    u = mda.Universe.empty(6, n_residues=6, trajectory=True)
    u.add_TopologyAttr("names", ["CA"] * 6)
    u.add_TopologyAttr("resnames", ["GLY"] * 6)
    u.add_TopologyAttr("resids", list(range(1, 7)))
    u.atoms.positions = np.asarray([
        [99, 0, 0], [0, 0, 0], [1, 0, 0],
        [0, 0, 0], [0, 1, 0], [0, 2, 0],
    ], dtype=np.float32)
    box = np.array([100., 100., 100., 90., 90., 90.])
    a = _backbone_axis(u.atoms[:3], box)
    b = _backbone_axis(u.atoms[3:], box)
    np.testing.assert_allclose(a, [1, 0, 0], atol=1e-6)
    np.testing.assert_allclose(b, [0, 1, 0], atol=1e-6)
    assert _angle(a, b) == pytest.approx(90., abs=1e-6)


def test_lipid_residence_tracks_molecules_and_censors():
    from mdscope.mechanisms.motif_residence import summarize_lipid_residence
    times = np.array([0., 10., 20., 30., 40., 50.])
    events = pd.DataFrame([
        {"time_ps": t, "motif": "proline_hydrophobic_A_1",
         "molecule_id": 13, "species": "GM1a"}
        for t in [0., 10., 30., 40., 50.]
    ])
    episodes = summarize_lipid_residence(events, times)
    assert len(episodes) == 2
    assert episodes["duration_ps"].tolist() == [20., 30.]
    assert episodes["n_frames"].tolist() == [2, 3]
    assert episodes["left_censored"].tolist() == [True, False]
    assert episodes["right_censored"].tolist() == [False, True]
    assert episodes["occupancy_fraction"].iloc[0] == pytest.approx(5 / 6)
    with pytest.raises(ValueError, match="Duplicate"):
        summarize_lipid_residence(pd.concat([events, events.iloc[:1]]), times)
    with pytest.raises(ValueError, match="align"):
        summarize_lipid_residence(events.assign(time_ps=events.time_ps + 1), times)
    assert summarize_lipid_residence(events.iloc[:0], times).empty


def test_shared_states_exports_replica_conditioned_nxip(tmp_path):
    from mdscope.mechanisms.shared_states import fit_shared_states
    cfg = {"analyses": {"sdc4_mechanisms": {"enabled": True},
                         "shared_states": {"enabled": True, "features": [
                             "sdc4_ectodomain_A_height_A",
                             "sdc4_ectodomain_B_height_A"],
                             "n_states": 2, "seed": 5,
                             "bootstrap_draws": 30,
                             "training_frames_per_replicate": 20,
                             "min_complete_fraction": 1.0}}}
    paths = {}
    for sys, shift in (("with", 0.0), ("without", 1.0)):
        paths[sys] = []
        for r in (1, 2):
            folder = tmp_path / sys / f"rep{r}"
            dest = folder / "extensions/sdc4_mechanisms"
            dest.mkdir(parents=True)
            height = np.r_[np.linspace(1, 2, 10),
                           np.linspace(7, 8, 10)] + shift
            pd.DataFrame({
                "time_ps": np.arange(20) * 10.,
                "sdc4_ectodomain_A_height_A": height,
                "sdc4_ectodomain_B_height_A": height * 1.5,
                "sdc4_motif_NxIP_PROB_1_steric_occlusion_fraction":
                    np.r_[np.full(10, .8), np.full(10, .2)],
                "sdc4_motif_NxIP_PROB_1_hs_PROB_S309_contact_state": np.r_[np.ones(10), np.zeros(10)],
                "sdc4_ectodomain_B_membrane_contact_state": np.r_[np.ones(10), np.zeros(10)],
                "sdc4_hs_A1_membrane_contact_state": np.zeros(20, dtype=int),
                "sdc4_hs_A1_hydration_water_per_heavy_atom": np.ones(20),
            }).to_parquet(dest / "timeseries.parquet", index=False)
            pd.DataFrame({
                "time_ps": np.arange(20) * 10.,
                "region": ["TM_A"] * 20,
                "species": ["GM1a" if sys == "with" else "CER160"] * 20,
                "leaflet": ["upper"] * 20,
                "moiety": ["any"] * 20,
                "contact_fraction": np.full(20, 0.25 * r),
                "composition_enrichment": np.full(20, 1.5 * r),
                "local_composition_enrichment": np.full(20, float(r)),
                "available_molecules": np.full(20, 4),
            }).to_parquet(dest / "lipid_region_specificity.parquet", index=False)
            paths[sys].append(folder)
    output = fit_shared_states(cfg, tmp_path, paths)
    cond = pd.read_parquet(output / "state_conditioned_nxip.parquet")
    assert set(cond["replicate"]) == {"rep1", "rep2"}
    assert len(cond) == 8
    assert (cond["frames_in_state"] > 0).all()
    comp = pd.read_parquet(output / "state_conditioned_nxip_comparisons.parquet")
    assert len(comp) == 2
    assert (comp["n_replicates_0"] == 2).all()
    assert set(comp["inference"]) == {"exploratory; no multiplicity-adjusted test"}
    hs = pd.read_parquet(output / "nxip_hs_contact_conditioning.parquet")
    assert len(hs) == 4
    assert set(hs["hs_chain"]) == {"PROB_S309"}
    membrane = pd.read_parquet(
        output / "nxip_membrane_contact_conditioning.parquet"
    )
    assert len(membrane) == 8
    subset = membrane[membrane.chain.eq("B")]
    assert (subset.frames == 10).all()
    expected = {0: .2, 1: .8}
    for row in subset.itertuples():
        assert row.mean_nxip_occlusion == pytest.approx(
            expected[row.membrane_contact_state]
        )

    assert (hs["frames_contact"] == 10).all()
    assert (hs["frames_no_contact"] == 10).all()
    np.testing.assert_allclose(hs["occlusion_contact_mean"], .8)
    np.testing.assert_allclose(hs["occlusion_no_contact_mean"], .2)

    exposure = pd.read_parquet(
        output / "nxip_state_exposure_threshold_sensitivity.parquet"
    )
    assert set(exposure.threshold) == {.25, .5, .75}
    null = pd.read_parquet(output / "nxip_state_circular_shift_null.parquet")
    assert len(null) == 8
    assert (null.n_shifts == 19).all()
    assert null.observed_difference.notna().all()
    assert (null.shift_null_ci_low <= null.shift_null_ci_high).all()

    assert exposure["exposed_fraction_given_state"].between(0, 1).all()
    mid = exposure[exposure.threshold.eq(.5)]
    assert set(mid["exposed_fraction_given_state"]) == {0., 1.}

    lipid_rep = pd.read_parquet(output / "lipid_specificity_replicate_statistics.parquet")
    assert len(lipid_rep) == 4
    assert set(lipid_rep.species) == {"GM1a", "CER160"}
    sys_stats = pd.read_parquet(output / "lipid_specificity_system_statistics.parquet")
    row = sys_stats.loc[
        sys_stats.system.eq("with")
        & sys_stats.species.eq("GM1a")
        & sys_stats.measure.eq("mean_contact_fraction")
    ].iloc[0]
    assert row.n_replicates == 2
    assert row.replicate_mean == pytest.approx(.375)
    assert row.replicate_bootstrap_ci_low >= .25
    assert row.replicate_bootstrap_ci_high <= .5



def test_intermittent_lipid_residence_gap_and_censoring():
    from mdscope.mechanisms.motif_residence import (
        summarize_lipid_residence, summarize_intermittent_residence
    )
    times = np.array([0., 10., 20., 30., 40., 50.])
    events = pd.DataFrame({
        "time_ps": [0., 10., 30., 40., 50.],
        "motif": ["region_TM_A"] * 5,
        "molecule_id": [42] * 5,
        "species": ["GM1a"] * 5,
    })
    continuous = summarize_lipid_residence(events, times)
    intermittent = summarize_intermittent_residence(events, times)
    assert continuous.duration_ps.tolist() == [20., 30.]
    assert intermittent.duration_ps.tolist() == [60.]
    assert intermittent.n_frames.tolist() == [5]
    assert intermittent.left_censored.tolist() == [True]
    assert intermittent.right_censored.tolist() == [True]
    assert intermittent.occupancy_fraction.iloc[0] == pytest.approx(5 / 6)
    assert len(summarize_intermittent_residence(events, times, gap_frames=0)) == 2
    with pytest.raises(ValueError, match="gap_frames"):
        summarize_intermittent_residence(events, times, gap_frames=-1)


def test_intermittent_residence_no_assumed_uniform_sampling():
    from mdscope.mechanisms.motif_residence import summarize_intermittent_residence
    times = np.array([0., 4., 10., 25., 30.])
    events = pd.DataFrame({
        "time_ps": [0., 10., 25.],
        "motif": ["region_TM_B"] * 3,
        "molecule_id": [17] * 3,
        "species": ["BGAL"] * 3,
    })
    result = summarize_intermittent_residence(events, times)
    assert result.n_frames.tolist() == [3]
    assert result.duration_ps.tolist() == [30.]
    assert result.right_censored.tolist() == [False]


def test_lipid_composition_enrichment_independent_toy_denominators():
    from mdscope.mechanisms.sdc4_mechanisms import _composition_enrichment
    # 1 GM1 out of 2 contacts, against 1 GM1 in 10 available lipids:
    # 0.5 / 0.1 == 5.0 enrichment. Other lipid: 0.5 / 0.9.
    table = pd.DataFrame({
        "time_ps": [0., 0., 0., 10., 10.],
        "region": ["TM_A"] * 5,
        "leaflet": ["upper"] * 5,
        "species": ["GM1a", "POPC", "GM1a", "GM1a", "POPC"],
        "moiety": ["any", "any", "glycan", "any", "any"],
        "contacted_molecules": [1, 1, 1, 0, 0],
        "available_molecules": [1, 9, 1, 0, 10],
        "locally_available_molecules": [1, 3, 1, 0, 4],
    })
    global_result = _composition_enrichment(table, "available_molecules")
    local_result = _composition_enrichment(table, "locally_available_molecules")
    assert global_result[0] == pytest.approx(5.0)
    assert global_result[1] == pytest.approx(5 / 9)
    assert np.isnan(global_result[2])
    assert local_result[0] == pytest.approx(2.)
    assert local_result[1] == pytest.approx(2 / 3)
    assert np.isnan(global_result[3]) and np.isnan(global_result[4])
    assert np.isnan(local_result[3]) and np.isnan(local_result[4])
