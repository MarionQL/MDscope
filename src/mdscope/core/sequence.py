"""Residue identity: one-letter codes, chain labels and sequence-motif matching.

Topology-only helpers shared by core selection logic, observables and mechanisms.
"""
from __future__ import annotations

import re

import numpy as np
from MDAnalysis.core.groups import Residue
from MDAnalysis.exceptions import NoDataError

AA3 = dict(zip(
    ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"],
    "ARNDCQEGHILKMFPSTWYV",
))
AA3.update({name: "H" for name in ("HSD", "HSE", "HSP", "HID", "HIE", "HIP")})


def chain_label(atom) -> str:
    """Atom chain identity, with the historical segment/unknown fallback."""
    return str(getattr(atom, "chainID", "") or getattr(atom, "segid", "") or "?")


def residue_label(entity) -> str:
    """Use the same atom-derived identity for atoms and topology residues."""
    atom = entity.atoms[0] if isinstance(entity, Residue) else entity
    return f"{chain_label(atom)}:{entity.resname}{entity.resid}"


def sequence_chain_label(atom):
    """Use a nonblank chainID, otherwise a nonblank segment label."""
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
        raise ValueError("Motif and glycan chain labels require chainID or segid")
    return chain


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
        chain = sequence_chain_label(atom)
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
                # Match the original expression at each start: wrapping a
                # capture would break backreferences and inline flags.
                hits = (compiled.match(seq, start) for start in range(len(seq))) if allow_overlaps else compiled.finditer(seq)
                for hit in hits:
                    if hit is None:
                        continue
                    start, end = hit.span()
                    if end <= start:
                        raise ValueError(f"Motif pattern {motif} must consume residues")
                    if "X" in seq[start:end]:
                        continue
                    ordinal[motif] += 1
                    key = f"{motif}_{label}_{ordinal[motif]}"
                    matches[key] = np.asarray(
                        [a.resindex for a in fragment[start:end]], dtype=int
                    )
    return matches
