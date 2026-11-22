"""Shared, PBC-aware primitives for the optional molecular-mechanism analyses.

Coordinates are never aligned: membrane-relative measurements must use the
physical trajectory. All distances are Angstrom; time is ps.
"""
from __future__ import annotations

import numpy as np
from MDAnalysis.lib.distances import capped_distance, distance_array

from ..core.atom_selection import select_heavy_atoms
from ..core.geometry import minimum_image_displacements


def selected(universe, selection, *, required=False):
    if not selection:
        if required:
            raise ValueError("Required MDAnalysis selection is missing")
        return universe.atoms[[]]
    atoms = universe.select_atoms(str(selection))
    if required and not len(atoms):
        raise ValueError(f"Selection matched no atoms: {selection}")
    return atoms


def heavy(atoms):
    return select_heavy_atoms(atoms, policy="mechanistic")


def nearest(a, b, box):
    if not len(a) or not len(b):
        return np.nan
    # Batching avoids large all-pairs allocations for long glycan chains.
    return float(min(np.min(distance_array(a[i:i+128], b, box=box))
                     for i in range(0, len(a), 128)))


def contact_count(a, b, box, cutoff=5.0):
    if not len(a) or not len(b):
        return 0
    pairs = capped_distance(a, b, max_cutoff=cutoff, box=box, return_distances=False)
    return int(len(pairs))


def planar_distances(centers, anchor, box, normal=2):
    offsets = minimum_image_displacements(np.asarray(centers) - np.asarray(anchor), box)
    offsets[:, normal] = 0.0
    return np.linalg.norm(offsets, axis=1)


def unit(vec):
    arr = np.asarray(vec, dtype=float)
    norm = float(np.linalg.norm(arr))
    return arr / norm if norm > 1e-8 else np.full(3, np.nan)


def orientation_vector(atoms):
    # Robust enough for short helices when atom ordering follows chain topology.
    pos = atoms.positions.astype(float)
    if len(pos) < 3:
        return np.full(3, np.nan)
    _, _, v = np.linalg.svd(pos - pos.mean(axis=0), full_matrices=False)
    axis = v[0]
    if np.dot(axis, pos[-1] - pos[0]) < 0:
        axis = -axis
    return unit(axis)


def describe_timeseries(rows, group_rows):
    """Return wide tables, preserving one row per time/group for shared MDscope stats."""
    import pandas as pd
    whole = pd.DataFrame(rows)
    grouped = pd.DataFrame(group_rows)
    if not grouped.empty:
        grouped = grouped.groupby(["group", "time_ps"], sort=False, as_index=False).first()
    return whole, grouped
