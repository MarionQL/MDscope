"""Acidic/basic protein residue salt-bridge detector."""
from __future__ import annotations
import numpy as np
from MDAnalysis.lib.distances import capped_distance
from .interaction_tables import residue_label
# Interaction detectors
NEGATIVE_ATOMS = "(resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2)"
POSITIVE_ATOMS = "(resname LYS and name NZ) or (resname ARG and name NH1 NH2)"


def initialize_salt_bridges(universe, cfg: dict):
    spec = cfg["analyses"]["salt_bridges"]
    if not spec.get("enabled", False):
        return None
    return {
        "negative": universe.select_atoms(f"protein and ({NEGATIVE_ATOMS})"),
        "positive": universe.select_atoms(f"protein and ({POSITIVE_ATOMS})"),
        "cutoff_A": float(spec.get("cutoff_A", 4.0)),
    }


def salt_bridge_events(state, time_ps: float, box=None) -> list[tuple]:
    # capped_distance returns only pairs inside the cutoff, avoiding
    # construction of the full acidicxbasic distance matrix.
    if not state or not len(state["negative"]) or not len(state["positive"]):
        return []
    neg, pos = state["negative"], state["positive"]
    pairs, distances = capped_distance(
        neg.positions,
        pos.positions,
        max_cutoff=state["cutoff_A"],
        box=box,
        return_distances=True
    )
    best = {}
    for (i, j), distance in zip(pairs, distances):
        pair = tuple(sorted((residue_label(neg[int(i)]), residue_label(pos[int(j)]))))
        best[pair] = min(float(distance), best.get(pair, np.inf))
    return [
        (
            time_ps, r1, r2, r1.split(":", 1)[0], r2.split(":", 1)[0],
            "intrachain" if r1.split(":", 1)[0] == r2.split(":", 1)[0] else "interchain",
            (
                f"intra_{r1.split(':', 1)[0]}"
                if r1.split(":", 1)[0] == r2.split(":", 1)[0]
                else "interchain"
            ),
            distance,
        )
        for (r1, r2), distance in best.items()
    ]
