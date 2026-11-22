"""Shared rigid-body frame alignment and selection validation (not an observable)."""
from __future__ import annotations

import numpy as np


def _valid_box(universe):
    dims = np.asarray(
        universe.dimensions,
        float
    ) if universe.dimensions is not None else np.array([])
    return universe.dimensions if len(dims) >= 3 and np.all(dims[:3] > 0) else None


def _fit_transform(
    mobile_fit: np.ndarray,
    reference_fit: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the unweighted rigid-body transform mobile -> reference."""
    mobile_fit = np.asarray(mobile_fit, dtype=np.float64)
    reference_fit = np.asarray(reference_fit, dtype=np.float64)
    if mobile_fit.shape != reference_fit.shape:
        raise ValueError(
            "Alignment coordinate arrays must have identical shapes: "
            f"mobile={mobile_fit.shape}, reference={reference_fit.shape}"
        )
    if mobile_fit.ndim != 2 or mobile_fit.shape[1] != 3:
        raise ValueError(
            f"Alignment coordinates must have shape (N, 3), got {mobile_fit.shape}"
        )
    if len(mobile_fit) < 3:
        raise ValueError("Alignment requires at least three atoms")
    mobile_center = mobile_fit.mean(axis=0)
    reference_center = reference_fit.mean(axis=0)
    mobile_centered = mobile_fit - mobile_center
    reference_centered = reference_fit - reference_center
    covariance = mobile_centered.T @ reference_centered
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0.0:
        u[:, -1] *= -1.0
        rotation = u @ vt
    return mobile_center, rotation, reference_center


def _apply_transform(
    coordinates: np.ndarray,
    transform: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    """Apply a previously calculated rigid-body transform."""
    mobile_center, rotation, reference_center = transform
    coordinates = np.asarray(coordinates, dtype=np.float64)
    return (coordinates - mobile_center) @ rotation + reference_center


def selection_details(atoms):
    """Topology-derived preview with a bounded atom-identity working set.

    Hash bytes match json.dumps(full_identity_list, ensure_ascii=False) so
    saved preflight identity digests retain their existing meaning.
    """
    import hashlib
    import json

    digest = hashlib.sha256()
    digest.update(b"[")
    residue_labels = {}
    segments = {}
    for index, atom in enumerate(atoms):
        segid = str(atom.segid)
        resid = int(atom.resid)
        resname = str(atom.resname)
        name = str(atom.name)
        if index:
            digest.update(b", ")
        digest.update(json.dumps(
            (segid, resid, resname, name), ensure_ascii=False
        ).encode("utf-8"))
        residue_labels.setdefault(f"{segid}:{resname}{resid}", None)
        segments.setdefault(segid, None)
    digest.update(b"]")
    return {
        "atoms": len(atoms), "residues": len(atoms.residues),
        "residue_labels": list(residue_labels),
        "segments": list(segments),
        "first_atoms": [{"index": int(atom.index), "segid": str(atom.segid),
                         "resid": int(atom.resid), "resname": str(atom.resname),
                         "name": str(atom.name)} for atom in atoms[:5]],
        "ordered_identity_sha256": digest.hexdigest(),
    }


def resolve_structural_atoms(universe, cfg, requirements):
    """Do not parse unused defaults; preflight and streaming share these roles."""
    selections = cfg["selections"]
    alignment = (cfg["analyses"]["alignment"].get("selection")
                 or selections.get("alignment") or "protein and name CA")
    expressions = {"protein": selections["protein"],
                   "features": selections.get("features") or alignment,
                   "alignment": alignment}
    return ({role: universe.select_atoms(expression) if requirements[role] else universe.atoms[[]]
             for role, expression in expressions.items()}, requirements)
