"""Shared, PBC-aware primitives for the optional molecular-mechanism analyses.

Coordinates are never aligned: membrane-relative measurements must use the
physical trajectory. All distances are Angstrom; time is ps.
"""
from __future__ import annotations

import re
import numpy as np
from MDAnalysis.exceptions import NoDataError
from ..core.atom_selection import select_heavy_atoms
from MDAnalysis.lib.distances import capped_distance, distance_array
from ..core.geometry import minimum_image_displacements

AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))
AA3.update({name: "H" for name in ("HSD", "HSE", "HSP", "HID", "HIE", "HIP")})


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


def sequence_matches(ca_atoms, patterns, *, allow_overlaps=True):
    """Match protein motifs within contiguous topology-ordered chain fragments.

    Returns motif IDs mapped to ordered topology residue indices (resindices).
    Labels retain the legacy motif_chain_ordinal form when chain names are
    unambiguous. A CA atom must represent each residue exactly once.

    Missing/modified amino acids are represented as X, not silently treated
    as their closest standard residue. Hits containing X are rejected.
    """
    groups = {}
    for atom in ca_atoms:
        if str(getattr(atom, "name", "CA")) != "CA":
            raise ValueError("Motif matching requires CA-only protein atoms")
        try:
            chain = str(getattr(atom, "chainID", "") or "").strip()
        except NoDataError:
            chain = ""
        try:
            segid = str(getattr(atom, "segid", "") or "").strip()
        except NoDataError:
            segid = ""
        if not chain:
            chain = segid
        if not chain:
            raise ValueError("Motif matching requires chainID or segid")
        group_id = (int(getattr(atom, "segindex", 0)), chain)
        groups.setdefault(group_id, []).append(atom)

    # Two distinct segments may carry the same chain label; avoid collisions
    # without changing existing IDs for normal single-segment SDC4 chains.
    counts = {}
    for _, chain in groups:
        counts[chain] = counts.get(chain, 0) + 1
    matches = {}
    for (segindex, chain), atoms in groups.items():
        atoms.sort(key=lambda atom: int(atom.resindex))
        indices = [int(atom.resindex) for atom in atoms]
        if len(indices) != len(set(indices)):
            raise ValueError(f"Duplicate CA residue in motif selection: {chain}")
        label = chain if counts[chain] == 1 else f"{chain}_seg{segindex}"
        ordinal = {motif: 0 for motif in patterns}
        # A selection gap is a fragment boundary, never a motif bridge.
        fragments = []
        for atom in atoms:
            if not fragments or int(atom.resindex) != int(fragments[-1][-1].resindex) + 1:
                fragments.append([])
            fragments[-1].append(atom)
        for fragment in fragments:
            seq = "".join(AA3.get(str(a.resname).upper(), "X") for a in fragment)
            for motif, pattern in patterns.items():
                try:
                    compiled = re.compile(str(pattern))
                except re.error as exc:
                    raise ValueError(f"Invalid motif regex {motif}: {pattern}") from exc
                # Lookaheads provide overlapping hits; a normal finditer
                # provides nonoverlapping hits. Both must consume residues.
                expression = re.compile(f"(?=({pattern}))") if allow_overlaps else compiled
                for hit in expression.finditer(seq):
                    start = hit.start()
                    end = start + len(hit.group(1)) if allow_overlaps else hit.end()
                    if end <= start or "X" in seq[start:end]:
                        continue
                    ordinal[motif] += 1
                    key = f"{motif}_{label}_{ordinal[motif]}"
                    matches[key] = np.asarray(
                        [a.resindex for a in fragment[start:end]], dtype=int
                    )
    return matches


def describe_timeseries(rows, group_rows):
    """Return wide tables, preserving one row per time/group for shared MDscope stats."""
    import pandas as pd
    whole = pd.DataFrame(rows)
    grouped = pd.DataFrame(group_rows)
    if not grouped.empty:
        grouped = grouped.groupby(["group", "time_ps"], sort=False, as_index=False).first()
    return whole, grouped
