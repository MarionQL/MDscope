"""Shared geometry and output contracts for the opt-in extended membrane pass.

All positions are in Å, all times in ps.  Every per-frame scalar is written
with the same wide time-series contract used by analysis.statistics.
"""
from __future__ import annotations
from dataclasses import dataclass
import re
import numpy as np
from MDAnalysis.lib.distances import capped_distance
from ..core.geometry import minimum_image_displacements

def valid_box(u):
    box = np.asarray(u.dimensions, dtype=float)
    if box.shape != (6,) or np.any(box[:3] <= 0) or not np.all(np.isfinite(box)):
        raise ValueError("Extended analyses require valid periodic box dimensions")
    return box

def heavy(group):
    return group[~np.asarray([str(x).upper().startswith("H") for x in group.names])]

def pairs(left, right, cutoff, box, exclude_self=False):
    """Sparse atom contacts without constructing a dense atom×atom matrix."""
    if len(left) == 0 or len(right) == 0:
        return np.empty((0, 2), dtype=np.intp)
    ij = capped_distance(left.positions.astype(np.float32),
                         right.positions.astype(np.float32),
                         max_cutoff=float(cutoff), box=box,
                         return_distances=False)
    if exclude_self and len(ij):
        ij = ij[left.indices[ij[:, 0]] != right.indices[ij[:, 1]]]
    return ij

def centers(atoms):
    """One equal-weight anchor center per residue, in stable resindex order.

    Anchor selections should not span periodic boundaries *within* a residue.
    Prefer a single lipid headgroup anchor when possible.
    """
    if len(atoms) == 0:
        return np.empty((0,), int), np.empty((0, 3), float)
    indices, inv = np.unique(atoms.resindices, return_inverse=True)
    n = len(indices)
    out = np.column_stack([np.bincount(inv, weights=atoms.positions[:, j],
                              minlength=n) for j in range(3)])
    out /= np.bincount(inv, minlength=n)[:, None]
    return indices, out

def minimum_image(delta, box):
    return minimum_image_displacements(delta, box, dtype=np.float32)

def lateral_distance(xyz, reference, box):
    delta = minimum_image(xyz - reference, box)
    return np.linalg.norm(delta[:, :2], axis=1)

def cosine_angle(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b)/n,-1,1)))) if n > 1e-8 else np.nan

def axis_from_ca(atoms):
    """Oriented principal axis; sign follows C-alpha sequence order."""
    ca = atoms.select_atoms("name CA")
    if len(ca) < 4:
        return None
    pos = ca.positions.astype(float)
    x = pos-pos.mean(axis=0)
    _, _, vh = np.linalg.svd(x, full_matrices=False)
    axis = vh[0]
    return axis*np.sign(np.dot(axis, pos[-1]-pos[0]) or 1)

def sequence_hits(protein, expression):
    """Find motif residue objects per chain, handling noncontiguous global resid.

    regex syntax supports G...G.  Returned labels reference actual topology
    resids, not canonical UniProt numbering.
    """
    chains = {}
    for res in protein.residues:
        ag = res.atoms
        chain = str(ag.chainIDs[0]) if hasattr(ag, "chainIDs") else str(ag.segids[0])
        chains.setdefault(chain, []).append(res)
    from MDAnalysis.lib.util import convert_aa_code
    hits = []
    for chain, residues in chains.items():
        letters=[]
        for res in residues:
            try: letters.append(convert_aa_code(str(res.resname).upper()))
            except Exception: letters.append("X")
        sequence="".join(letters)
        for match in re.finditer("(?="+expression+")", sequence):
            start=match.start()
            length=len(re.match(expression, sequence[start:]).group())
            subset=residues[start:start+length]
            hits.append((chain, int(subset[0].resid), int(subset[-1].resid),
                         subset, sequence[start:start+length]))
    return hits

def safe_id(s):
    return re.sub(r"[^0-9A-Za-z_]+", "_", str(s)).strip("_")

@dataclass
class LipidAnchors:
    names: np.ndarray
    resindices: np.ndarray
    upper: np.ndarray
    coords: np.ndarray

def load_lipid_anchors(universe, lipid_specs):
    names, indices, coords = [], [], []
    for species, selector in lipid_specs.items():
        atoms=universe.select_atoms(str(selector))
        if not len(atoms):
            continue
        ix, xy=centers(atoms)
        names.extend([str(species)]*len(ix))
        indices.extend(ix.tolist())
        coords.extend(xy.tolist())
    if not names:
        raise ValueError("No lipid anchors matched. Check analyses.extended.lipids")
    if len(set(indices)) != len(indices):
        raise ValueError("Lipid anchor selections overlap on a residue")
    xyz=np.asarray(coords,float)
    middle=np.median(xyz[:,2])
    upper=xyz[:,2] >= middle
    if not np.any(upper) or not np.any(~upper):
        raise ValueError("Could not assign both leaflets from lipid anchors")
    return LipidAnchors(np.asarray(names),np.asarray(indices),upper,xyz)

def update_lipid_anchors(u, lipid_specs, known):
    coords=np.empty_like(known.coords)
    mapping={int(i):n for n,i in enumerate(known.resindices)}
    for _,selector in lipid_specs.items():
        ix,xyz=centers(u.select_atoms(str(selector)))
        for j,k in enumerate(ix):
            coords[mapping[int(k)]]=xyz[j]
    return coords
