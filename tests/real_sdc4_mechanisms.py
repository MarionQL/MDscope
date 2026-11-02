"""Four-frame SDC4 geometry/lipid-map smoke on the published real TPR/XTCs.

HS attachment mapping is intentionally *not* guessed. The synthetic suite
covers HS hydration/contact semantics; real HS analysis requires verified
per-chain attachment-site selections from the user's actual topology.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd

from mdscope.mechanisms.sdc4_mechanisms import SDC4Mechanisms

DATA = Path(os.environ.get("MDSCOPE_REAL_DATA", "real_test_data")).resolve()
ROOT = Path(os.environ.get("MDSCOPE_REAL_OUTPUT", "real_test_results")).resolve()

SYSTEMS = {
    "carb_sdc4": ("carb-sdc4-step7_production.tpr", "carb-sdc4-rep1_test.xtc"),
    "nocarb_sdc4": ("nocarb-sdc4-step7_production.tpr", "nocarb-sdc4-rep1_test.xtc"),
}


def run():
    for system, (top, traj) in SYSTEMS.items():
        u = mda.Universe(str(DATA / top), str(DATA / traj))
        u.trajectory[0]
        species = {}
        for name, selection in (
            ("POPC", "resname POPC and name P"),
            ("POPE", "resname POPE and name P"),
            ("POPI", "resname POPI and name P"),
            ("DOPC", "resname DOPC and name P"),
            ("DPPC", "resname DPPC and name P"),
            ("PSM", "resname PSM and name P"),
            ("CHL1", "resname CHL1 and name O3"),
            ("GM1a", "chainID GM1a and resname CER160 and name NF"),
            ("BGAL", "chainID BGAL and resname CER160 and name NF"),
            ("CER160", "chainID CER160 and resname CER160 and name NF"),
        ):
            if len(u.select_atoms(selection)):
                species[name] = selection
        assert len(species) >= 3, species
        tm = {
            "A": "protein and chainID PROA and resid 150:170 and name CA",
            "B": "protein and chainID PROB and resid 419:439 and name CA",
        }
        ecto = {
            "A": "protein and chainID PROA and resid 1:149 and name CA",
            "B": "protein and chainID PROB and resid 270:418 and name CA",
        }
        for selection in list(tm.values()) + list(ecto.values()):
            assert len(u.select_atoms(selection)) >= 4, (
                f"{system}: missing real SDC4 selection {selection}"
            )
        for chain in ("A", "B"):
            ecto_atoms = u.select_atoms(ecto[chain])
            tm_atoms = u.select_atoms(tm[chain])
            assert len(tm_atoms) == 21, (system, chain, len(tm_atoms))
            assert len(ecto_atoms) >= 140, (system, chain, len(ecto_atoms))
            assert not np.intersect1d(ecto_atoms.resindices, tm_atoms.resindices).size, (
                system, chain, "overlapping ectodomain / TM selections"
            )
        # Check the independent residue-membrane distance selection uses
        # precisely the configured ectodomain residues and no TM residues.
        rmd = u.select_atoms(
            "protein and name CA and ((chainID PROA and resid 1:149) "
            "or (chainID PROB and resid 270:418))"
        )
        ecto_ids = np.union1d(
            u.select_atoms(ecto["A"]).resindices,
            u.select_atoms(ecto["B"]).resindices,
        )
        np.testing.assert_array_equal(
            np.unique(rmd.resindices), ecto_ids,
            err_msg=f"{system}: residue_membrane_distance target drift",
        )
        # Release-data diagnostic: report the complete SDC4 chain sequence and
        # selected TM segment independently, so a missing GGIVG cannot be
        # confused with zero physical interhelical contacts.
        for chain, sel in tm.items():
            ca_tm = u.select_atoms(sel)
            ca_full = u.select_atoms(
                f"protein and chainID PRO{chain} and name CA"
            )
            from mdscope.mechanisms.mechanistic_common import AA3
            tm_seq = "".join(AA3.get(a.resname.upper(), "X") for a in ca_tm)
            full_seq = "".join(AA3.get(a.resname.upper(), "X") for a in ca_full)
            print(
                f"REAL GGIVG {system}/{chain}: TM={tm_seq!r} "
                f"TM_resids={ca_tm.resids.tolist()} "
                f"full_GGIVG_offsets={[i for i in range(len(full_seq)-4) if full_seq[i:i+5] == 'GGIVG']}",
                flush=True,
            )
            ca_ecto = u.select_atoms(ecto[chain])
            ecto_seq = "".join(AA3.get(a.resname.upper(), "X") for a in ca_ecto)
            if not ecto_seq.endswith("TEVLAA"):
                raise AssertionError(
                    f"{system}/{chain}: ectodomain does not end in TEVLAA; "
                    f"observed ending={ecto_seq[-15:]!r}"
                )
            if tm_seq != "LIVGGIVGILFAVFLILLLMY":
                raise AssertionError(
                    f"{system}/{chain}: TM sequence {tm_seq!r} does not match "
                    "the user-verified domain boundary LIVGGIVGILFAVFLILLLMY"
                )
        spec = {
            "enabled": True, "axis": "z",
            "tm_chains": tm, "ectodomains": ecto,
            "hs_chains": {}, "hs_attachment_sites": {},
            "hs_auto": {
                "enabled": True,
                "glycan_selection": "resname AGLCNA AIDOA BGLCA BGAL BXYL and chainID PROA PROB",
                "attachment_selection": "protein and resname SER and name OG",
                "expected_chains": 6,
            },
            "water_oxygen_selection": "resname TIP3 and name OH2",
            "ion_selections": {"K": "resname POT", "Cl": "resname CLA"},
            "lipid_anchors": species,
            # Use only headgroup anchors for this small smoke, not a full
            # membrane-contact definition.
            "lipid_contact_selection": "(resname POPC POPE POPI POPS DOPC DPPC PSM CHL1) or (chainID GM1a BGAL CER160)",
            "motif_selection": "protein and name CA",
            "motifs": {"NxIP": "N.IP"},
            "require_motifs": False, "require_gxxxg": True,
            "map_extent_A": 40., "map_bin_A": 4.,
            "local_radius_A": 25., "local_min_lipids": 8,
            "local_order": {"enabled": False},
        }
        cfg = {
            "analyses": {"sdc4_mechanisms": spec},
            "selections": {"protein": "protein"},
            "trajectory": {"start": 0, "stop": 4, "stride": 1, "chunk_frames": 2},
            "time": {"dt_ps": None, "origin_ps": 0.},
        }
        out = ROOT / "sdc4_mechanisms_smoke" / system
        observer = SDC4Mechanisms(u, cfg, out)
        observer.prune_from(0)
        for i in range(4):
            observer.update(u.trajectory[i], i)
            if i % 2 == 1:
                observer.flush(i // 2)
        result = observer.finalize()
        df = pd.read_parquet(out / "extensions/sdc4_mechanisms/timeseries.parquet")
        assert len(df) == 4
        assert np.all(np.isfinite(df["sdc4_tm_crossing_angle_deg"]))
        assert np.all(np.isfinite(df["sdc4_ectodomain_A_height_A"]))
        assert result["n_hs_chains"] == 6
        assert "sdc4_hs_PROA_S39_hydration_water_count" in df
        assert "sdc4_hs_PROB_S333_lipid_contact_lipids" in df
        with np.load(out / "extensions/sdc4_mechanisms/lipid_density_maps.npz") as maps:
            assert int(maps["nframes"]) == 4
            assert maps["counts"].sum() > 0
        print(f"REAL SDC4 SMOKE {system}: PASS {result}", flush=True)
        print(f"REAL SDC4 AUDIT {system}: " +
              (out / "extensions/sdc4_mechanisms/selection_audit.json").read_text(),
              flush=True)


if __name__ == "__main__":
    run()
