"""Protein–lipid/HS local environment estimator shared by execution paths.

All measurements use source coordinates, Angstrom cutoffs and topology residue
indices. Scientific geometry is defined here; fused and standalone runners
only own trajectory traversal and persistence.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import distance_array

from ..biophysical_common import heavy, box_of, residue_pairs, sequence_motifs
from ..config import atomic_parquet, ensure_dir, time_ps


def prepare_local_environment(universe, cfg):
    """Resolve topology selections once, before frame iteration."""
    settings = cfg["analyses"]["local_environment"]
    protein = heavy(universe.select_atoms(
        settings.get("protein_selection") or cfg["selections"]["protein"]
    ))
    lipids = heavy(universe.select_atoms(settings.get(
        "lipid_selection",
        "resname POPC POPE POPI POPS DOPC DPPC PSM CHL1 CER160 CER1 BGAL GM1a",
    )))
    if not len(protein) or not len(lipids):
        raise ValueError("No protein/lipid atoms matched local_environment selections")
    glycans = heavy(universe.select_atoms(settings.get(
        "glycan_selection", "resname AGLCNA AIDOA BGLCA BGAL BXYL"
    )))
    ions = heavy(universe.select_atoms(settings.get(
        "ion_selection", "resname SOD NA CLA CL POT K CAL CA"
    )))
    hits = sequence_motifs(protein.residues, settings.get("motifs", {
        "nxip": "N.IP", "gxxxg": "G...G", "hs_candidate": "SG",
    }))
    groups = {}
    for hit in hits:
        key = "motif_" + hit["motif"] + "_" + hit["chain"] + "_" + str(hit["start_resid"])
        groups[key] = heavy(universe.atoms[
            np.isin(universe.atoms.resindices, hit["resindices"])
        ])
    grouped = {
        str(name): lipids[lipids.resnames == name]
        for name in np.unique(lipids.resnames)
    }
    metadata = [
        {"motif_id": key, "sequence": "",
         "resids": ",".join(map(str, group.residues.resids))}
        for key, group in groups.items()
    ]
    return {
        "prot": protein, "lip": lipids, "glycan": glycans, "ions": ions,
        "groups": groups, "grouped": grouped, "meta": metadata,
        "cutoff": float(settings.get("contact_cutoff_A", 5.)),
    }


def local_environment_frame(state, universe, time):
    """Calculate one frame without writing outputs or modifying the Universe.

    Returns (scalar row, residue-contact rows); absent species retain the
    established missing-column semantics, not manufactured zeros.
    """
    protein, lipids, cutoff = state["prot"], state["lip"], state["cutoff"]
    box = box_of(universe)
    row = {"time_ps": time}
    contacts = residue_pairs(protein, lipids, cutoff, box)
    touched = np.unique(contacts[:, 0]) if len(contacts) else np.empty(0, int)
    row["local_protein_lipid_contact_residues_count"] = len(touched)
    per_residue = []
    if len(contacts):
        residues, counts = np.unique(contacts[:, 0], return_counts=True)
        per_residue = [
            (time, int(resid), int(count))
            for resid, count in zip(residues, counts)
        ]
    for species, atoms in state["grouped"].items():
        pairs = residue_pairs(protein, atoms, cutoff, box)
        row["local_" + species + "_contact_residues_count"] = (
            len(np.unique(pairs[:, 0])) if len(pairs) else 0
        )
    for name, group in state["groups"].items():
        if not len(group):
            continue
        distances = distance_array(group.positions, lipids.positions, box=box)
        row[name + "_min_lipid_distance_A"] = float(distances.min())
        row[name + "_lipid_contact_fraction"] = float(
            np.mean(distances.min(axis=1) <= cutoff)
        )
        others = protein[~np.isin(protein.indices, group.indices)]
        if len(others):
            internal = distance_array(group.positions, others.positions, box=box)
            # Geometric occlusion, not solvent-accessible surface burial.
            row[name + "_protein_occlusion_fraction"] = float(
                np.mean(np.any(internal < cutoff, axis=1))
            )
    if len(state["glycan"]):
        row["local_hs_lipid_contact_pairs_count"] = len(
            residue_pairs(state["glycan"], lipids, cutoff, box)
        )
        if len(state["ions"]):
            row["local_hs_ion_contacts_count"] = len(
                residue_pairs(state["glycan"], state["ions"], 3.5, box)
            )
    return row, per_residue


def run_local_environment(rep, cfg, repout):
    """Standalone trajectory reader; reuses the same per-frame scientific kernel."""
    import MDAnalysis as mda

    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    state = prepare_local_environment(universe, cfg)
    rows, per_residue = [], []
    tr = cfg["trajectory"]
    for frame_index, ts in enumerate(universe.trajectory[
        slice(tr.get("start"), tr.get("stop"), int(tr.get("stride", 1)))
    ]):
        row, contacts = local_environment_frame(
            state, universe, time_ps(ts, frame_index, cfg)
        )
        rows.append(row)
        per_residue.extend(contacts)
    out = ensure_dir(Path(repout) / "extensions" / "local_environment")
    frame = pd.DataFrame(rows)
    atomic_parquet(out / "timeseries.parquet", frame)
    if per_residue:
        atomic_parquet(out / "residue_contacts.parquet", pd.DataFrame(
            per_residue, columns=["time_ps", "resindex", "n_lipid_residues"],
        ))
    if state["meta"]:
        atomic_parquet(out / "motifs.parquet", pd.DataFrame(state["meta"]))
    return {"frames": len(rows), "metrics": len(frame.columns) - 1}
