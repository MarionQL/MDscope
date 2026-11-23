"""Complete interprotomer contact and local-normal geometry references."""
import numpy as np
import pandas as pd
import pytest

from mdscope.core.lipid_geometry import local_leaflet_normal
from mdscope.mechanisms.protomer_contacts import (
    InterprotomerContacts, summarize_interprotomer_contacts,
)
from mdscope.mechanisms.ceramide_scaffold import pool_ceramide_scaffolds

import MDAnalysis as mda


def _synthetic():
    universe = mda.Universe.empty(4, n_residues=4,
                                  atom_resindex=[0, 1, 2, 3], trajectory=True)
    universe.add_TopologyAttr("names", ["CA"] * 4)
    universe.add_TopologyAttr("resnames", ["ALA", "GLY", "SER", "LEU"])
    universe.add_TopologyAttr("resids", [1, 2, 1, 2])
    universe.add_TopologyAttr("chainIDs", ["PROA", "PROA", "PROB", "PROB"])
    universe.atoms.positions = np.array([[0., 0., 0.], [0., 0., 10.],
                                         [8., 0., 0.], [8., 0., 10.]])
    universe.dimensions = np.array([100., 100., 100., 90., 90., 90.])
    return universe


def test_complete_interprotomer_interface_has_zero_residues_and_pair_occupancy(tmp_path):
    u = _synthetic()
    u.trajectory[0]
    a = u.select_atoms("protein and chainID PROA")
    b = u.select_atoms("protein and chainID PROB")
    tm = {"A": a.select_atoms("name CA"), "B": b.select_atoms("name CA")}
    contacts = InterprotomerContacts(u.select_atoms("protein"), tm, cutoff=5.)
    # Same-z positions are 4 A apart. Every protein residue is included.
    b.positions = b.positions - np.array([4., 0., 0.])
    near_pairs, near_residues = contacts.frame(u.dimensions, 0.)
    assert len(near_pairs) > 0
    assert len(near_residues) == len(a.residues) + len(b.residues)
    assert set(near_residues.chain) == {"A", "B"}
    assert near_pairs.residue_A.str.startswith("PROA:").all()
    assert near_pairs.residue_B.str.startswith("PROB:").all()
    # No contacts in a deliberately separated second frame.
    b.positions = b.positions + np.array([30., 0., 0.])
    far_pairs, far_residues = contacts.frame(u.dimensions, 10.)
    assert far_pairs.empty
    assert (far_residues.any_contact == 0).all()
    pair_parts, residue_parts = [], []
    for i, (pairs, residues) in enumerate(
        [(near_pairs, near_residues), (far_pairs, far_residues)]
    ):
        pair_path = tmp_path / f"pairs-{i}.parquet"
        residue_path = tmp_path / f"residues-{i}.parquet"
        pairs.to_parquet(pair_path, index=False)
        residues.to_parquet(residue_path, index=False)
        pair_parts.append(pair_path)
        residue_parts.append(residue_path)
    summarize_interprotomer_contacts(tmp_path / "results", pair_parts, residue_parts,
                                     [0., 10.])
    pair_summary = pd.read_parquet(tmp_path / "results/interprotomer_pair_summary.parquet")
    residue_summary = pd.read_parquet(tmp_path / "results/interprotomer_residue_summary.parquet")
    assert (pair_summary.contact_occupancy == .5).all()
    assert len(residue_summary) == len(a.residues) + len(b.residues)
    assert (residue_summary.frames == 2).all()
    assert (residue_summary.contact_occupancy < 1).all()


def test_local_membrane_normal_tilt_and_periodic_translation():
    box = np.array([100., 100., 100., 90., 90., 90.])
    xy = np.array([(x, y) for x in (-10., 0., 10.) for y in (-10., 0., 10.)])
    anchors = np.column_stack((xy, 20. + .1 * xy[:, 0]))
    fitted = local_leaflet_normal(np.zeros(3), anchors, box, 20.)
    target = np.array([-.1, 0., 1.])
    target /= np.linalg.norm(target)
    np.testing.assert_allclose(fitted, target, atol=1e-6)
    translated = anchors + [100., 0., 0.]
    np.testing.assert_allclose(
        local_leaflet_normal(np.zeros(3), translated, box, 20.),
        fitted, atol=1e-6,
    )
    assert np.isnan(local_leaflet_normal(
        np.zeros(3), anchors[:2], box, 20.,
    )).all()


def test_ceramide_scaffold_pools_counts_before_ratios_and_leaves_absent_unavailable():
    raw = pd.DataFrame([
        dict(time_ps=0., region="TM_A", leaflet="upper", moiety="ceramide_head",
             species=sp, available_molecules=total,
             locally_available_molecules=local, contacted_molecules=hit,
             locally_contacted_molecules=l_hit)
        for sp, total, local, hit, l_hit in [
            ("GM1a", 2, 2, 1, 1), ("BGAL", 8, 4, 2, 1),
        ]
    ])
    pooled = pool_ceramide_scaffolds(raw)
    assert len(pooled) == 1
    assert pooled.cohort.iloc[0] == "GM1a+BGAL"
    assert pooled.contact_fraction.iloc[0] == pytest.approx(3./10.)
    assert pooled.matched_local_contact_fraction.iloc[0] == pytest.approx(2./6.)
    assert not (pooled.cohort == "CER160").any()
    raw.loc[2] = [0., "TM_A", "upper", "ceramide_head", "CER160", 5, 5, 0, 0]
    with_cer = pool_ceramide_scaffolds(raw)
    assert with_cer.loc[with_cer.cohort == "CER160", "contact_fraction"].iloc[0] == 0.
