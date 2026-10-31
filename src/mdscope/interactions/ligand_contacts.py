"""Ligand-protein pairwise contacts and closest-distance detector."""
from __future__ import annotations
import numpy as np
from MDAnalysis.lib.distances import capped_distance,distance_array
from .interaction_tables import atom_label,residue_label
def initialize_ligand_contacts(universe, protein, ligand, cfg: dict):
    spec = cfg["analyses"]["ligand_contacts"]
    if not spec.get("enabled", False) or ligand is None or not len(ligand):
        return None
    return {
        "protein": protein,
        "ligand": ligand,
        "cutoff_A": float(spec.get("cutoff_A", 4.0)),
        "save_atom_events": bool(spec.get("save_atom_events", True)),
        "ligand_labels": np.asarray(
            [atom_label(atom) for atom in ligand],
            dtype=object,
        ),
        "protein_atom_labels": np.asarray(
            [atom_label(atom) for atom in protein],
            dtype=object,
        ),
        "protein_residue_labels": np.asarray(
            [residue_label(atom) for atom in protein],
            dtype=object,
        ),
    }


def _minimum_distance(left, right, box=None, block_atoms: int = 4096) -> float:
    """Return the exact minimum distance without a full giant matrix."""
    if not len(left) or not len(right):
        return np.nan
    minimum = np.inf
    block_atoms = max(1, int(block_atoms))
    for start in range(0, len(right), block_atoms):
        stop = min(start + block_atoms, len(right))
        distances = distance_array(
            left.positions,
            right.positions[start:stop],
            box=box,
        )
        minimum = min(minimum, float(np.min(distances)))
    return minimum


def ligand_contact_events(state, time_ps: float, box=None):
    if not state:
        return np.nan, 0, []
    ligand, protein = state["ligand"], state["protein"]
    pairs, distances = capped_distance(
        ligand.positions,
        protein.positions,
        max_cutoff=state["cutoff_A"],
        box=box,
        return_distances=True,
    )
    if not len(pairs):
        return _minimum_distance(ligand, protein, box), 0, []

    li = pairs[:, 0].astype(int)
    pi = pairs[:, 1].astype(int)
    protein_residues = state["protein_residue_labels"]
    rows = []
    if state.get("save_atom_events", True):
        rows = list(
            zip(
                np.full(len(pairs), time_ps),
                state["ligand_labels"][li],
                state["protein_atom_labels"][pi],
                protein_residues[pi],
                distances.astype(float),
            )
        )
    return (
        float(np.min(distances)),
        int(len(np.unique(protein_residues[pi]))),
        rows,
    )
