"""PBC-safe membrane anchor and protein-center primitives.

Species selections must pick one anchor *residue* per physical lipid. For
CHARMM-GUI glycosphingolipids use a ceramide anchor (e.g. chainID GM1a and
resname CER160 and name NF), not carbohydrate residues.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from MDAnalysis.lib.distances import minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors
from .geometry import minimum_image_displacements


@dataclass
class AnchorTracker:
    """Topology-fixed atom-to-lipid mapping; no per-frame AtomGroup intersections."""
    atoms: object
    owner: np.ndarray
    counts: np.ndarray

    def __len__(self):
        return len(self.counts)

    def positions(self):
        xyz = self.atoms.positions.astype(float, copy=False)
        # Reconstruct selected atoms relative to an anchor atom of each lipid.
        # This avoids an incorrect Cartesian midpoint when a lipid crosses PBC.
        first = np.full(len(self.counts), len(xyz), dtype=int)
        np.minimum.at(first, self.owner, np.arange(len(xyz)))
        if np.any(first == len(xyz)):
            raise ValueError("Lipid anchor has no selected atoms")
        reference = xyz[first]
        box = self.atoms.universe.dimensions
        if box is None or not np.all(np.asarray(box[:3]) > 0):
            raise ValueError("Lipid anchors require valid periodic box dimensions")
        unwrapped = reference[self.owner] + minimize_vectors(
            xyz - reference[self.owner], box
        )
        return np.column_stack([
            np.bincount(self.owner, weights=unwrapped[:, axis], minlength=len(self.counts))
            / self.counts for axis in range(3)
        ])


def classify_leaflets(z, length):
    """Circular two-cluster partition; bilayers may cross the periodic z edge."""
    z = np.asarray(z, float)
    if len(z) < 4 or length <= 0:
        raise ValueError("Need >=4 lipid anchors and positive box length")
    angle = (2 * np.pi * z / length) % (2 * np.pi)
    unit = np.stack([np.cos(angle), np.sin(angle)], axis=1)
    seed0 = unit[0]
    seed1 = unit[np.argmin(unit @ seed0)]
    centers = np.stack([seed0, seed1])
    assignment = np.zeros(len(z), dtype=int)
    for _ in range(30):
        assignment = np.argmax(unit @ centers.T, axis=1)
        if np.any(np.bincount(assignment, minlength=2) < 2):
            raise ValueError("Cannot resolve two populated bilayer leaflets")
        fresh = np.stack([unit[assignment == i].mean(axis=0) for i in range(2)])
        norms = np.linalg.norm(fresh, axis=1)
        if np.any(norms < 1e-8):
            raise ValueError("Ambiguous leaflet circular centers")
        fresh /= norms[:, None]
        if np.allclose(fresh, centers, rtol=0, atol=1e-7):
            break
        centers = fresh
    # Orientation is conventional; for temporal continuity align subsequent
    # classifications to the initial leaflet labels by maximum overlap.
    angle_centers = np.arctan2(centers[:, 1], centers[:, 0])
    delta = (angle_centers[1] - angle_centers[0] + np.pi) % (2 * np.pi) - np.pi
    upper_label = 1 if delta > 0 else 0
    return assignment == upper_label


def aligned_leaflets(z, length, reference):
    """Independent current-frame classification aligned to reference lipid IDs.

    Returns current labels and the number of lipids whose leaflet changed.
    Alignment prevents an arbitrary upper/lower label swap when the bilayer
    drifts across a periodic boundary.
    """
    reference = np.asarray(reference, bool)
    current = classify_leaflets(z, length)
    if current.shape != reference.shape:
        raise ValueError("Leaflet reference size differs from current anchors")
    if np.count_nonzero(current != reference) > np.count_nonzero(~current != reference):
        current = ~current
    return current, int(np.count_nonzero(current != reference))


def periodic_center(positions, dimensions):
    """Circular center of a wrapped atom group in a periodic triclinic box.

    Cartesian positions are converted to fractional lattice coordinates,
    averaged on the unit circle along each lattice direction, then mapped
    back. The result is invariant to independent whole-box translations of
    individual atoms. Requires a compact group (not spread uniformly around
    a periodic axis); for diffuse groups use a local reference selection.
    """
    xyz = np.asarray(positions, float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) == 0:
        raise ValueError("positions must be a nonempty (N,3) array")
    box = np.asarray(dimensions, float)
    if box.shape != (6,) or np.any(~np.isfinite(box)) or np.any(box[:3] <= 0):
        raise ValueError("A valid periodic box with six dimensions is required")
    lattice = np.asarray(triclinic_vectors(box), float)
    frac = xyz @ np.linalg.inv(lattice)
    angle = 2 * np.pi * frac
    sine = np.sin(angle).mean(axis=0)
    cosine = np.cos(angle).mean(axis=0)
    strength = np.hypot(sine, cosine)
    if np.any(strength < 1e-5):
        raise ValueError("Protein atoms span an ambiguous periodic center; narrow protein_selection")
    mean_frac = (np.arctan2(sine, cosine) / (2 * np.pi)) % 1.0
    return mean_frac @ lattice


def anchor_groups(u, spec):
    """Vectorized per-residue lipid anchors and initial-frame leaflet labels."""
    selections = spec.get("species_selections") or {}
    if not selections:
        raise ValueError("Define species_selections with one anchor per lipid")
    names, resids, indices, owners = [], [], [], []
    used_residues = set()
    for name, selection in selections.items():
        atoms = u.select_atoms(str(selection))
        if not len(atoms):
            if spec.get("allow_absent_species", False):
                continue
            raise ValueError(f"Lipid species {name!r} selection matched zero atoms: {selection}")
        unique_resids, inverse = np.unique(atoms.resindices, return_inverse=True)
        overlap = used_residues.intersection(unique_resids.tolist())
        if overlap:
            raise ValueError(f"species_selections overlap on {len(overlap)} lipid residues")
        used_residues.update(unique_resids.tolist())
        offset = len(names)
        names.extend([str(name)] * len(unique_resids))
        resids.extend(unique_resids.tolist())
        indices.append(atoms.indices)
        owners.append(inverse + offset)
    if not indices:
        raise ValueError("No lipid anchors matched; check species selections")
    owner = np.concatenate(owners).astype(np.int32, copy=False)
    tracker = AnchorTracker(
        u.atoms[np.concatenate(indices)],
        owner,
        np.bincount(owner, minlength=len(names)),
    )
    coords = tracker.positions()
    if not np.all(np.isfinite(coords)):
        raise ValueError("Non-finite anchor coordinates")
    normal = {"x": 0, "y": 1, "z": 2}[str(spec.get("axis", "z")).lower()]
    box = np.asarray(u.dimensions, float)
    if box.shape != (6,) or np.any(box[:3] <= 0):
        raise ValueError("Lipid anchors require a periodic simulation box")
    leaf = classify_leaflets(coords[:, normal], float(box[normal]))
    return np.asarray(names, str), np.asarray(resids, int), tracker, leaf


def centers(groups):
    return groups.positions() if isinstance(groups, AnchorTracker) else np.stack(
        [g.center_of_geometry() for g in groups]
    )


def planar_offsets(points, reference, box, normal=2):
    delta = minimum_image_displacements(
        np.asarray(points, float) - np.asarray(reference, float), box,
    )
    return np.delete(delta, normal, axis=-1)


def lateral_distances(points, reference, box, normal=2):
    return np.linalg.norm(planar_offsets(points, reference, box, normal), axis=-1)
