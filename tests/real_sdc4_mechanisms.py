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
import yaml

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
            "B": "protein and chainID PROB and resid 420:440 and name CA",
        }
        ecto = {
            "A": "protein and chainID PROA and resid 1:149 and name CA",
            "B": "protein and chainID PROB and resid 270:419 and name CA",
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
            "or (chainID PROB and resid 270:419))"
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
        from mdscope.mechanisms.mechanistic_common import sequence_matches
        for chain in ("A", "B"):
            ca = u.select_atoms(f"protein and chainID PRO{chain} and name CA")
            pro_hits = sequence_matches(ca, {"proline_hydrophobic": "VVHPLVPL"})
            nxip_hits = sequence_matches(ca, {"NxIP": "N.IP"})
            assert len(pro_hits) == 1, (
                system, chain, "VVHPLVPL mapping", list(pro_hits)
            )
            assert len(nxip_hits) >= 1, (
                system, chain, "NxIP mapping", list(nxip_hits)
            )
            print(
                f"REAL EXTRACELLULAR MOTIFS {system}/{chain}: "
                f"VVHPLVPL={list(pro_hits.values())[0].tolist()} "
                f"NxIP={[v.tolist() for v in nxip_hits.values()]}",
                flush=True,
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
            "motifs": {"NxIP": "N.IP", "proline_hydrophobic": "VVHPLVPL"},
            "require_motifs": False, "require_gxxxg": True,
            "motif_lipid_residence": True,
            "map_extent_A": 40., "map_bin_A": 4.,
            "local_radius_A": 25., "local_min_lipids": 8,
            "local_order": {"enabled": False},
        }
        production = yaml.safe_load(
            Path("synd4_full_mechanisms.yaml").read_text()
        )["analyses"]["sdc4_mechanisms"]
        spec["lipid_moiety_selections"] = production["lipid_moiety_selections"]
        chemical = {
            label: u.select_atoms(selection)
            for label, selection in spec["lipid_moiety_selections"].items()
        }
        assert all(len(group) for group in chemical.values()), (
            system, {k: len(v) for k, v in chemical.items()}
        )
        for left, left_group in chemical.items():
            for right, right_group in chemical.items():
                if left < right:
                    assert not np.intersect1d(
                        left_group.indices, right_group.indices
                    ).size, (system, left, right)
        assert set(chemical["ceramide_tail"].resnames) == {"CER160"}
        assert set(chemical["ceramide_head"].resnames) == {"CER160"}
        assert not (set(chemical["glycan"].resnames) & {"CER160"})
        cfg = {
            "analyses": {"sdc4_mechanisms": spec},
            "selections": {"protein": "protein"},
            "trajectory": {"start": 0, "stop": 4, "stride": 1, "chunk_frames": 2},
            "time": {"dt_ps": None, "origin_ps": 0.},
        }
        out = ROOT / "sdc4_mechanisms_smoke" / system
        # Inspect force-field-specific CER160 chemistry before choosing
        # tail versus polar-head selections; avoid name-prefix guesses.
        ceramide = u.select_atoms("resname CER160")
        if len(ceramide):
            first_ceramide = ceramide.residues[0]
            print(
                f"REAL CER160 ATOM NAMES {system}: "
                f"{first_ceramide.atoms.names.tolist()}",
                flush=True,
            )
        observer = SDC4Mechanisms(u, cfg, out)
        observer.prune_from(0)
        def independent_residue_pairs(left, right, cutoff=5.0):
            # Independent brute-force orthorhombic minimum-image oracle for
            # the real released TPR; not the production neighbor search.
            pairs = set()
            xyz_b = right.positions.astype(float)
            box = u.dimensions[:3].astype(float)
            for atom, resid in zip(left.positions, left.resindices):
                delta = xyz_b - atom
                delta -= np.rint(delta / box) * box
                neighbors = np.flatnonzero(np.einsum("ij,ij->i", delta, delta) <= cutoff**2)
                pairs.update((int(resid), int(right.resindices[k])) for k in neighbors)
            return len(pairs)

        reference = []
        for i in range(4):
            u.trajectory[i]
            reference.append((
                independent_residue_pairs(observer.tm_heavy["A"], observer.tm_heavy["B"]),
                independent_residue_pairs(observer.gxxxg["A"], observer.gxxxg["B"]),
            ))
            observer.update(u.trajectory[i], i)
            if i % 2 == 1:
                observer.flush(i // 2)
        result = observer.finalize()
        df = pd.read_parquet(out / "extensions/sdc4_mechanisms/timeseries.parquet")
        assert len(df) == 4
        assert np.all(np.isfinite(df["sdc4_tm_crossing_angle_deg"]))
        assert np.all(np.isfinite(df["sdc4_ectodomain_A_height_A"]))
        assert result["n_hs_chains"] == 6
        from mdscope.mechanisms.motif_residence import COLUMNS
        events_path = out / "extensions/sdc4_mechanisms/motif_lipid_contact_events.parquet"
        residence_path = out / "extensions/sdc4_mechanisms/motif_lipid_residence.parquet"
        assert events_path.exists() and residence_path.exists()
        intermittent_path = (
            out / "extensions/sdc4_mechanisms/motif_lipid_residence_intermittent.parquet"
        )
        assert intermittent_path.exists()
        intermittent = pd.read_parquet(intermittent_path)
        assert set(COLUMNS).issubset(intermittent.columns)
        assert (intermittent["n_frames"] > 0).all()

        assert set(COLUMNS).issubset(pd.read_parquet(residence_path).columns)
        residence_audit = json.loads(
            (out / "extensions/sdc4_mechanisms/selection_audit.json").read_text()
        )
        assert residence_audit["motif_lipid_residence_unmapped_fragments"] >= 0
        molecule_counts = residence_audit["motif_lipid_residence_molecules_by_species"]
        if system == "carb_sdc4":
            assert molecule_counts.get("GM1a", 0) > 0
            assert molecule_counts.get("BGAL", 0) > 0
        else:
            assert molecule_counts.get("CER160", 0) > 0
        contact_events = pd.read_parquet(events_path)
        assert set(contact_events["species"]).issubset(molecule_counts)


        specificity = pd.read_parquet(
            out / "extensions/sdc4_mechanisms/lipid_region_specificity.parquet"
        )
        assert not specificity.empty
        assert {"GM1a", "BGAL"}.issubset(
            set(specificity.species)
        ) if system == "carb_sdc4" else "CER160" in set(specificity.species)
        assert (specificity.contacted_molecules <=
                specificity.available_molecules).all()
        available = specificity.available_molecules > 0
        assert specificity.loc[~available, "contact_fraction"].isna().all()
        assert specificity.loc[available, "contact_fraction"].between(0, 1).all()
        assert len(specificity.time_ps.unique()) == 4
        assert {"locally_available_molecules", "min_contact_distance_A",
                "local_composition_enrichment"}.issubset(specificity.columns)
        assert (specificity["locally_available_molecules"] >= 0).all()
        assert (specificity["min_contact_distance_A"].dropna() <= 5.0 + 1e-5).all()
        zero_local = specificity.loc[
            specificity["moiety"].eq("any")
            & specificity["locally_available_molecules"].eq(0)
        ]
        assert zero_local["local_composition_enrichment"].isna().all()
        assert "composition_enrichment" in specificity
        absent = specificity.available_molecules.eq(0)
        assert specificity.loc[absent, "composition_enrichment"].isna().all()
        assert (specificity["moiety"] == "any").any()
        assert set(specificity.loc[specificity.species.eq("GM1a"), "moiety"]) == (
            {"any", "glycan", "ceramide"} if system == "carb_sdc4" else set()
        )

        assert any(col.endswith("_P4_phi_deg") for col in df.columns)
        assert any(col.endswith("_P7_psi_deg") for col in df.columns)

        actual_tm = df["sdc4_tm_interface_residue_pairs"].to_numpy()
        actual_motif = df["sdc4_gxxxg_interface_residue_pairs"].to_numpy()
        np.testing.assert_array_equal(actual_tm, np.asarray(reference)[:, 0])
        np.testing.assert_array_equal(actual_motif, np.asarray(reference)[:, 1])
        print(f"REAL TM GEOMETRY {system}: TM pairs={actual_tm.tolist()} "
              f"GGIVG pairs={actual_motif.tolist()}", flush=True)
        assert "sdc4_motif_proline_hydrophobic_PROA_1_steric_occlusion_fraction" in df
        assert "sdc4_motif_proline_hydrophobic_PROB_1_steric_occlusion_fraction" in df
        # Issue #101: chain-resolved NxIP and ectodomain signals must be
        # present in the same frame-indexed table for conditional analyses.
        assert "sdc4_ectodomain_B_height_A" in df
        assert "sdc4_motif_NxIP_PROB_1_steric_occlusion_fraction" in df
        assert df["time_ps"].is_unique
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
