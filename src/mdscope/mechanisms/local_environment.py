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

from ..biophysical_common import box_of, heavy, residue_pairs
from ..config import atomic_parquet, ensure_dir, time_ps
from ..core.motifs import resolve_motifs
from ..core.sequence import sequence_chain_label

DEFAULT_MOTIFS = {"nxip": "N.IP", "gxxxg": "G...G", "hs_candidate": "SG"}


def motif_groups(universe, settings, shared, protein_selection):
    """Heavy-atom groups for local-environment motifs, via the shared matcher.

    A regex string matches within ``protein_selection`` (non-overlapping, zero
    matches allowed); ``{motif: name}`` uses that ``selections.motifs``
    definition with its own scope and policies. Keys stay
    ``motif_<label>_<segid>_<first resid>`` (the saved column identities; the
    former matcher read a residue-level chainID that does not exist and so
    always fell back to segid), ordered by chain, motif, position.
    """
    configured = settings.get("motifs", DEFAULT_MOTIFS)
    definitions = {
        label: shared[value["motif"]] if isinstance(value, dict)
        else {"pattern": value, "selection": protein_selection,
              "overlap": False, "missing": "allow"}
        for label, value in configured.items()
    }
    if not definitions:
        return {}
    first = {}
    for residue in universe.select_atoms(protein_selection).residues:
        first.setdefault((sequence_chain_label(residue.atoms[0]), residue.segindex), residue.ix)
    rank = {label: i for i, label in enumerate(definitions)}
    hits = sorted(
        ((label, identity, hit) for label, found in resolve_motifs(universe, definitions).items()
         for identity, hit in found.items()),
        key=lambda item: (first.get((item[2]["chain"], item[2]["segindex"]), item[2]["resindices"][0]),
                          rank[item[0]], item[2]["resindices"][0]),
    )
    groups = {}
    for label, _, hit in hits:
        key = f"motif_{label}_{hit['segid'].strip() or 'unknown'}_{hit['resids'][0]}"
        if key in groups:
            raise ValueError(f"local_environment motif hits share the identity {key}; "
                             "give each chain a distinct segid")
        groups[key] = heavy(universe.residues[hit["resindices"]].atoms)
    return groups


def prepare_local_environment(universe, cfg):
    """Resolve topology selections once, before frame iteration."""
    settings = cfg.get("analyses", {}).get("local_environment", {})
    protein_selection = settings.get("protein_selection") or cfg["selections"]["protein"]
    protein = heavy(universe.select_atoms(protein_selection))
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
    groups = motif_groups(universe, settings, cfg.get("selections", {}).get("motifs", {}),
                          protein_selection)
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
