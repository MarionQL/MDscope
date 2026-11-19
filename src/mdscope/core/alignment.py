"""Shared rigid-body frame alignment and selection validation (not an observable)."""
from __future__ import annotations
import numpy as np
from ..config import analysis_groups
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
    """Topology-derived identity/preview; never infer counts from YAML text."""
    identities = [(str(atom.segid), int(atom.resid), str(atom.resname),
                   str(atom.name)) for atom in atoms]
    import hashlib
    import json
    return {
        "atoms": len(atoms), "residues": len(atoms.residues),
        "residue_labels": list(dict.fromkeys(
            f"{segid}:{resname}{resid}" for segid, resid, resname, _ in identities)),
        "segments": list(dict.fromkeys(item[0] for item in identities)),
        "first_atoms": [{"index": int(atom.index), "segid": str(atom.segid),
                         "resid": int(atom.resid), "resname": str(atom.resname),
                         "name": str(atom.name)} for atom in atoms[:5]],
        "ordered_identity_sha256": hashlib.sha256(
            json.dumps(identities, ensure_ascii=False).encode()).hexdigest(),
    }


def resolved_selection_counts(universe, cfg: dict) -> dict:
    alignment_spec = cfg["analyses"].get("alignment", {})
    alignment_enabled = bool(alignment_spec.get("enabled", True))
    atoms, requirements = resolve_structural_atoms(universe, cfg)
    protein, alignment, features = (atoms[role] for role in ("protein", "alignment", "features"))
    out = {
        "protein": len(protein) if requirements["protein"] else None,
        "alignment": len(alignment) if requirements["alignment"] else None,
        "features": len(features) if requirements["features"] else None,
        "required_by": requirements,
        "groups": {},
        "details": {role: selection_details(atoms[role]) if requirements[role]
                    else {"status": "not_required"} for role in atoms},
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
            "feature_status": "required" if requirements["features"] else "not_required",
            "alignment_atoms": len(align_atoms),
            "measurement": selection_details(atoms),
            "features": selection_details(measured),
            "alignment": selection_details(align_atoms),
        }
    return out


def structural_selection_requirements(cfg):
    """Global structural roles actually consumed by enabled production paths."""
    from ..config import analysis_runs
    analyses = cfg["analyses"]
    enabled = lambda name: bool(analyses.get(name, {}).get("enabled", False))
    features = [name for name in ("alignment", "rmsd", "rmsf") if enabled(name)]
    protein = [name for name in ("rmsd", "rmsf", "dssp", "salt_bridges", "ligand_contacts", "sdc4_mechanisms") if enabled(name)]
    if enabled("radius_of_gyration") and not analyses["radius_of_gyration"].get("selection"):
        protein.append("radius_of_gyration")
    if enabled("hydrogen_bonds") and analyses["hydrogen_bonds"].get("between") is None:
        protein.append("hydrogen_bonds")
    for name in ("local_environment", "hydration"):
        if enabled(name) and not analyses[name].get("protein_selection"):
            protein.append(name)
    if any(not (run.spec.get("target_selection") or run.spec.get("selection")) for run in analysis_runs(cfg, "sasa")):
        protein.append("sasa")
    if any(not run.spec.get("selection") for run in analysis_runs(cfg, "residue_membrane_distance")):
        protein.append("residue_membrane_distance")
    return {"protein": protein, "features": features,
            "alignment": ["alignment"] if enabled("alignment") else []}


def resolve_structural_atoms(universe, cfg):
    """Do not parse unused defaults; preflight and streaming share these roles."""
    requirements = structural_selection_requirements(cfg)
    selections = cfg["selections"]
    alignment = (cfg["analyses"]["alignment"].get("selection")
                 or selections.get("alignment") or "protein and name CA")
    expressions = {"protein": selections["protein"],
                   "features": selections.get("features") or alignment,
                   "alignment": alignment}
    return ({role: universe.select_atoms(expression) if requirements[role] else universe.atoms[[]]
             for role, expression in expressions.items()}, requirements)
