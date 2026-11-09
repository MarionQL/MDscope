"""Hand-built topology, periodic steric reference and RDF interpretation fixtures."""
import json

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from mdscope.config import DEFAULTS, Replicate, merge
from mdscope.mechanisms.potassium_bridge_audit import peak_rows, topology_sites
from mdscope.pipeline.runner import RunState, run_replicate


def test_audit_requires_bonded_charged_sites_and_tracks_periodic_steric_reference(tmp_path, monkeypatch):
    # Two ASP/Ser protomers, one attached carboxylated and sulfated sugar on
    # each, one K and one water. PRO* HS names sort before acid* patch names.
    names, masses, charges, atom_resindices, bonds, chains = [], [], [], [], [], []
    resnames = []

    def add(chain, resname, atoms, edges):
        start = len(names)
        for name, mass, charge in atoms:
            names.append(name)
            masses.append(mass)
            charges.append(charge)
            atom_resindices.append(len(resnames))
            chains.append(chain)
        resnames.append(resname)
        bonds.extend((start+a, start+b) for a, b in edges)
        return start

    for chain in ("PROA", "PROB"):
        add(chain, "ASP", [("CA",12,0),("CB",12,0),("CG",12,.7),
                           ("OD1",16,-.85),("OD2",16,-.85)], [(0,1),(1,2),(2,3),(2,4)])
        ser = add(chain, "SER", [("OG",16,-.2)], [])
        sugar = add(chain, "BGLCA", [("C5",12,0),("C6",12,.7),
                                    ("O61",16,-.85),("O62",16,-.85)], [(0,1),(1,2),(1,3)])
        sulfate = add(chain, "AGLCNA", [("C1",12,0),("OS",16,-.3),("S",32,1.4),
                                        ("O1",16,-.7),("O2",16,-.7),("O3",16,-.7)],
                      [(0,1),(1,2),(2,3),(2,4),(2,5)])
        bonds.extend([(ser,sugar),(sugar,sulfate)])
    potassium = add("K", "POT", [("POT",39.1,1.)], [])
    oxygen = add("W", "TIP3", [("OH2",16,-.8),("H1",1,.4),("H2",1,.4)], [(0,1),(0,2)])
    u = mda.Universe.empty(len(names), n_residues=len(resnames), atom_resindex=atom_resindices,
                           trajectory=True)
    for key, value in dict(names=names, masses=masses, charges=charges, chainIDs=chains,
                           resnames=resnames, resids=np.arange(len(resnames)), bonds=bonds).items():
        u.add_TopologyAttr(key, value)
    xyz = np.full((len(names), 3), 14., dtype=np.float32)
    xyz[3:5] = [1, 3, 3]
    xyz[8:10] = [3, 1, 3]
    xyz[potassium] = [3, 3, 3]
    # A periodic shift by 1/3 box reaches a distinct K location with a clash.
    xyz[12:16] = [13, 3, 3]
    xyz[4] = [12, 3, 3]
    xyz[9] = [13, 2, 3]
    xyz[oxygen:oxygen+3] = [[4, 3, 3], [3.1, 3, 3], [4, 4, 3]]
    u.load_new(xyz[None, ...], order="fac", dimensions=[30, 30, 30, 90, 90, 90])
    config = {"analyses": {
        "solvation": {"water_selection": "resname TIP3 and name OH2", "targets": {
            "motif": {"acidic_patch": {"selection": "protein and name CA", "sequence": "D"}}}},
        "sdc4_mechanisms": {"hs_auto": {"glycan_selection": "resname BGLCA AGLCNA",
            "attachment_selection": "resname SER and name OG", "min_glycan_residues": 2,
            "expected_chains": 2}},
    }}
    cfg = merge(DEFAULTS, config)
    for spec in cfg["analyses"].values():
        spec["enabled"] = False
    cfg["analyses"]["potassium_bridge_audit"]["enabled"] = True
    cfg["cache"]["aligned_coordinates"] = False
    cfg["time"]["dt_ps"] = 100.
    cfg["time"]["origin_ps"] = 200.
    rep = Replicate("test", 1, tmp_path / "fixture.tpr", tmp_path / "fixture.xtc", "rep1")
    monkeypatch.setattr(mda, "Universe", lambda *args: u)
    state = RunState(tmp_path)
    run_replicate(rep, cfg, tmp_path, state)
    output = tmp_path / "test/rep1/extensions/potassium_bridge_audit"
    topology = json.loads((output / "topology_audit.json").read_text())
    assert topology["site_counts"] == {"protein_carboxylate": 4, "hs_carboxylate": 4,
                                       "hs_sulfate": 6}
    assert topology["potassium_ions"] == 1
    frames = pd.read_parquet(output / "bridge_frames.parquet")
    assert frames.time_ps.unique().tolist() == [200.]
    first = frames.query("acid_patch == 'acid_PROA_1' and hs_chain == 'PROA_S1' "
                         "and hs_kind == 'hs_carboxylate' and contact_cutoff_A == 3")
    observed = first.query("shift_x == 0 and shift_y == 0")
    assert observed.bridging_ions.tolist() == [1, 1, 0, 0]
    shifted = first[np.isclose(first.shift_x, 1/3) & first.shift_y.eq(0)]
    assert shifted.eligible_ions.tolist() == [1, 0, 0, 0]
    assert shifted.bridging_ions.tolist() == [1, 0, 0, 0]
    rdf = pd.read_parquet(output / "rdf.parquet")
    assert rdf.query("kind == 'water_oxygen'").pair_count.sum() == 1
    u.atoms[potassium].charge = 0.
    with pytest.raises(ValueError, match="unbonded"):
        topology_sites(u, config)
    # A complete directory remains reusable, even with invalid source data.
    run_replicate(rep, cfg, tmp_path, state)
    (output / "rdf.parquet").unlink()
    with pytest.raises(RuntimeError, match="invalid required potassium_bridge_audit artifact"):
        run_replicate(rep, cfg, tmp_path, state)


def test_peak_candidates_have_intervening_minimum_and_allow_missing_second_peak():
    r = np.arange(.05, 8., .1)
    first = np.exp(-((r-2.8)/.3)**2)
    second = .5*np.exp(-((r-5.)/.4)**2)
    rdf = pd.DataFrame({"kind": "test", "g_box": first + second})
    for _, _, _, low, minimum, high in peak_rows(rdf, 100):
        assert low == pytest.approx(2.8, abs=.051)
        assert high == pytest.approx(5., abs=.051)
        assert low < minimum < high
    rdf["g_box"] = first
    assert all(np.isnan(row[-1]) and np.isnan(row[-2])
               for row in peak_rows(rdf, 100))
