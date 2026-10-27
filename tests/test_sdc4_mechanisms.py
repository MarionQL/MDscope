"""Synthetic PBC, output, restart and shared-state regression tests for SDC4."""
from pathlib import Path
from types import SimpleNamespace
import json

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.analysis.sdc4_mechanisms import SDC4Mechanisms, _local_surface
from mdscope.analysis.shared_states import fit_shared_states


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
                                "GLY", "ALA", "ALA", "ALA", "GLY", "ALA", "ALA"), 1):
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
        assert all(audit["gxxxg_resindices"][chain] for chain in ("A", "B"))
    a = tmp_path / "clean/extensions/sdc4_mechanisms"
    b = tmp_path / "resume/extensions/sdc4_mechanisms"
    pd.testing.assert_frame_equal(pd.read_parquet(a / "timeseries.parquet"),
                                  pd.read_parquet(b / "timeseries.parquet"))
    with np.load(a / "lipid_density_maps.npz") as left, np.load(b / "lipid_density_maps.npz") as right:
        np.testing.assert_array_equal(left["counts"], right["counts"])


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
