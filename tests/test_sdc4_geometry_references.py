"""Analytic SDC4 geometry, mapping failure, and persisted-output references.

Coordinates and expected values are specified in docs/developer/sdc4-fixtures.md.
No MDscope helper constructs the expected physical values.
"""
import json

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from mdscope.mechanisms.sdc4_mechanisms import (
    SDC4Mechanisms, _backbone_axis, _local_surface,
)


def _fixture(separation=3., transform="identity"):
    names, residues, chains, xyz = [], [], [], []
    # Deliberately put B before A and use unrelated residue numbers.
    for chain, x in (("B", separation), ("A", 0.)):
        for aa, z in zip(
            ["ASN", "ALA", "ILE", "PRO", "ALA", "GLY", "GLY", "ILE", "VAL", "GLY", "ALA", "ALA"],
            [36., 32., 28., 24., 20., 12., 8., 4., 0., -4., -8., -12.],
        ):
            names.append("CA"); residues.append(aa); chains.append(chain); xyz.append([x, 0., z])
    for z in (16., -16.):
        for x in (-10., 0., 10.):
            for y in (-10., 0., 10.):
                names.append("P"); residues.append("POPC"); chains.append("L"); xyz.append([x, y, z])
    for z in (26., 30.):
        names.append("C1"); residues.append("AGLCNA"); chains.append("HS"); xyz.append([0., 0., z])
    names.append("OH2"); residues.append("TIP3"); chains.append("W"); xyz.append([0., 0., 29.])
    n = len(names)
    u = mda.Universe.empty(n, n_residues=n, atom_resindex=np.arange(n), trajectory=True)
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("resnames", residues)
    u.add_TopologyAttr("resids", list(range(301, 313)) + list(range(71, 83)) + [1]*(n-24))
    u.add_TopologyAttr("chainIDs", chains)
    coords = np.asarray(xyz, float)
    if transform == "rotate_z":
        coords = coords[:, [1, 0, 2]] * [-1., 1., 1.]
    elif transform == "wrap_translate":
        # Bilayer crosses z=0; atoms independently shifted by lattice vectors.
        coords += [47., -23., 47.]
        coords += ((np.arange(n*3).reshape(n, 3) % 3)-1) * 100.
    u.atoms.positions = coords
    u.dimensions = [100., 100., 100., 90., 90., 90.]
    spec = {
        "enabled": True,
        "tm_chains": {c: f"chainID {c} and resindex {lo}:{lo+6}" for c, lo in (("A", 17), ("B", 5))},
        "ectodomains": {"A": "resindex 12:16", "B": "resindex 0:4"},
        "hs_chains": {"A1": "chainID HS"},
        "hs_attachment_sites": {"A1": "resindex 15"},
        "lipid_anchors": {"POPC": "resname POPC and name P"},
        "surface_anchor_species": ["POPC"], "lipid_contact_selection": "resname POPC",
        "water_oxygen_selection": "resname TIP3 and name OH2",
        "motif_selection": "protein and name CA", "motifs": {"NxIP": "N.IP"},
        "require_motifs": True, "require_gxxxg": True,
        "local_radius_A": 30., "local_min_lipids": 8, "contact_cutoff_A": 5.,
    }
    return u, {"selections": {"protein": "protein"}, "analyses": {"sdc4_mechanisms": spec}}


def _observe(u, cfg, out):
    obs = SDC4Mechanisms(u, cfg, out)
    obs.update(u.trajectory.ts, 0)
    direct = pd.DataFrame(obs.rows["timeseries"])
    obs.flush(0); obs.finalize()
    persisted = pd.read_parquet(out/"extensions/sdc4_mechanisms/timeseries.parquet")
    pd.testing.assert_frame_equal(persisted, direct)
    return persisted.iloc[0], obs


@pytest.mark.parametrize("separation,tm_pairs,motif_pairs,other_pairs", [(3., 19, 13, 14), (8., 0, 0, 0)])
@pytest.mark.parametrize("transform", ["identity", "rotate_z", "wrap_translate"])
def test_analytic_dimer_and_membrane_geometry(tmp_path, separation, tm_pairs, motif_pairs, other_pairs, transform):
    u, cfg = _fixture(separation, transform)
    row, obs = _observe(u, cfg, tmp_path)
    assert row.sdc4_tm_crossing_angle_deg == pytest.approx(0., abs=1e-6)
    assert row.sdc4_tm_interface_residue_pairs == tm_pairs
    assert row.sdc4_gxxxg_interface_residue_pairs == motif_pairs
    assert row.sdc4_gxxxg_interface_atom_contacts == motif_pairs
    for chain in ("A", "B"):
        assert row[f"sdc4_gxxxg_{chain}_to_other_tm_residue_pairs"] == other_pairs
        assert row[f"sdc4_ectodomain_{chain}_height_A"] == pytest.approx(12., abs=1e-5)
        assert row[f"sdc4_ectodomain_{chain}_tilt_deg"] == pytest.approx(0., abs=1e-6)
    assert row.sdc4_local_thickness_A == pytest.approx(32., abs=1e-5)
    for leaf in ("upper", "lower"):
        assert row[f"sdc4_local_{leaf}_curvature_Ainv"] == pytest.approx(0., abs=1e-8)
    assert row.sdc4_hs_A1_extension_A == pytest.approx(6., abs=1e-5)
    assert row.sdc4_hs_A1_orientation_cos == pytest.approx(1.)
    assert row.sdc4_hs_A1_radius_of_gyration_A == pytest.approx(2.)
    assert row.sdc4_hs_A1_attachment_lipid_contact_pairs == 0
    assert row.sdc4_hs_A1_membrane_contact_state == 0
    assert obs.audit["gxxxg_resindices"] == {"A": list(range(17, 22)), "B": list(range(5, 10))}


def test_optional_missing_motif_is_unavailable_not_physical_zero(tmp_path):
    u, cfg = _fixture(separation=8.)
    u.residues[6].resname = "ALA"
    cfg["analyses"]["sdc4_mechanisms"]["require_gxxxg"] = False
    row, _ = _observe(u, cfg, tmp_path)
    assert row.sdc4_tm_interface_residue_pairs == 0  # valid mapped chains
    assert row.sdc4_gxxxg_A_to_other_tm_residue_pairs == 0
    assert np.isnan(row.sdc4_gxxxg_B_to_other_tm_residue_pairs)
    assert np.isnan(row.sdc4_gxxxg_interface_atom_contacts)
    assert np.isnan(row.sdc4_gxxxg_interface_residue_pairs)
    audit = json.loads((tmp_path/"extensions/sdc4_mechanisms/selection_audit.json").read_text())
    assert audit["gxxxg_status"] == {"A": "valid", "B": "unavailable_missing_motif"}


def test_ambiguous_tm_motif_selection_is_rejected(tmp_path):
    u, cfg = _fixture()
    cfg["analyses"]["sdc4_mechanisms"]["tm_chains"]["A"] = "protein and resindex 5:11 17:23"
    with pytest.raises(ValueError, match="ambiguous.*GGIVG"):
        SDC4Mechanisms(u, cfg, tmp_path)


@pytest.mark.parametrize("coords", [np.zeros((4, 3)), [[1,0,0], [-1,0,0], [0,1,0], [0,-1,0]]])
def test_degenerate_tm_axis_is_unavailable(coords):
    u = mda.Universe.empty(4, trajectory=True)
    u.atoms.positions = coords
    assert np.isnan(_backbone_axis(u.atoms, [100.,100.,100.,90.,90.,90.])).all()


def test_sparse_or_rank_deficient_local_surface_is_unavailable(tmp_path):
    u, cfg = _fixture()
    cfg["analyses"]["sdc4_mechanisms"]["local_radius_A"] = 1.
    row, _ = _observe(u, cfg, tmp_path)
    assert np.isnan(row.sdc4_local_thickness_A)
    assert np.isnan(row.sdc4_local_upper_curvature_Ainv)
    assert np.isnan(_local_surface(np.zeros((9, 2)), np.zeros(9)))


def test_motif_sasa_requires_every_mapped_residue_at_each_time(tmp_path):
    u, cfg = _fixture()
    obs = SDC4Mechanisms(u, cfg, tmp_path)
    for i in range(3):
        u.trajectory.ts.time = float(i)
        obs.update(u.trajectory.ts, i)
    obs.flush(0)
    folder = tmp_path/"profiles/sasa"
    folder.mkdir(parents=True)
    # NxIP_A comprises resindices 12:15. A complete motif totals 10 nm²;
    # a missing or nonfinite residue is unavailable, never a smaller SASA.
    pd.DataFrame({"time_ps": [0.]*4 + [1.]*3 + [2.]*4,
                  "topology_residue_index": [12,13,14,15,12,13,14,12,13,14,15],
                  "sasa_nm2": [1.,2.,3.,4.,1.,2.,3.,1.,2.,3.,np.nan]}).to_parquet(folder/"residue_sasa.parquet")
    obs.finalize()
    data = pd.read_parquet(tmp_path/"extensions/sdc4_mechanisms/timeseries.parquet")
    values = data.sdc4_motif_NxIP_A_1_sasa_nm2.to_numpy()
    assert values[0] == 10.
    assert np.isnan(values[1:]).all()
