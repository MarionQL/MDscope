"""Shared rigid-body frame alignment and selection validation (not an observable)."""
from __future__ import annotations
import numpy as np
from .config import analysis_groups
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


def resolved_selection_counts(universe, cfg: dict) -> dict:
    selections = cfg["selections"]
    alignment_spec = cfg["analyses"].get("alignment", {})
    alignment_enabled = bool(alignment_spec.get("enabled", True))
    alignment_selection = (
        str(
            alignment_spec.get("selection") or selections.get("alignment") or "protein and name CA"
        )
    )
    protein = universe.select_atoms(selections["protein"])
    alignment = (
        universe.select_atoms(alignment_selection) if alignment_enabled else universe.atoms[[]]
    )
    features = universe.select_atoms(selections.get("features") or alignment_selection)
    out = {
        "protein": len(protein),
        "alignment": len(alignment),
        "features": len(features),
        "groups": {}
    }
    for group_id, group in analysis_groups(cfg).items():
        atoms = universe.select_atoms(group["selection"])
        measured = features.intersection(atoms)
        align_atoms = (
            universe.select_atoms(group["alignment_selection"]) if alignment_enabled
            and group.get("alignment_selection") else alignment
        )
        out["groups"][group_id] = {
            "group_atoms": len(atoms),
            "feature_atoms": len(measured),
            "alignment_atoms": len(align_atoms)
        }
    return out
