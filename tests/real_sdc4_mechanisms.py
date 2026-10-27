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

from mdscope.analysis.sdc4_mechanisms import SDC4Mechanisms

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
            "A": "protein and chainID PROA and resid 146:170 and name CA",
            "B": "protein and chainID PROB and resid 344:368 and name CA",
        }
        ecto = {
            "A": "protein and chainID PROA and resid 1:145 and name CA",
            "B": "protein and chainID PROB and resid 270:343 and name CA",
        }
        for selection in list(tm.values()) + list(ecto.values()):
            assert len(u.select_atoms(selection)) >= 4, (
                f"{system}: missing real SDC4 selection {selection}"
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
            "ion_selections": {"Na": "resname SOD"},
            "lipid_anchors": species,
            # Use only headgroup anchors for this small smoke, not a full
            # membrane-contact definition.
            "lipid_contact_selection": "(resname POPC POPE POPI POPS DOPC DPPC PSM CHL1) or (chainID GM1a BGAL CER160)",
            "motif_selection": "protein and name CA",
            "motifs": {"NxIP": "N.IP"},
            "require_motifs": False, "require_gxxxg": False,
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
