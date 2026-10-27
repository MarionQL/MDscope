"""Membrane APL, thickness and the original minimum-image cumulative MSD implementation.

One shared frame pass updates this state; no duplicate diffusion/MSD trajectory reader.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from MDAnalysis.lib.distances import minimize_vectors
)

AXIS_INDEX = {"x": 0, "y": 1, "z": 2}

# Membrane and density helpers live here because they are updated in the
# same frame loop as the other time series.
@dataclass
class MembraneState:
    enabled: bool = False
    axis: str = "z"
    ai: int = 2
    headgroups: object = None
    upper_resix: np.ndarray | None = None
    lower_resix: np.ndarray | None = None
    upper_atoms: object = None
    lower_atoms: object = None
    n_upper: int = 0
    n_lower: int = 0
    apl: bool = False
    thickness: bool = False
    msd: bool = False
    msd_include_3d: bool = False
    msd_residue_atoms: list | None = None
    msd_resindices: np.ndarray | None = None
    msd_upper_mask: np.ndarray | None = None
    msd_lower_mask: np.ndarray | None = None
    msd_last_positions: np.ndarray | None = None
    msd_displacements: np.ndarray | None = None


def _residue_centers(residue_atoms) -> np.ndarray:
    # Residue centers are calculated once per lipid/residue instead of
    # tracking every atom for membrane MSD.
    centers = []
    for atoms in residue_atoms:
        try:
            center = atoms.center_of_mass()
        except (AttributeError, ValueError, TypeError):
            center = atoms.center_of_geometry()
        if not np.all(np.isfinite(center)):
            center = atoms.center_of_geometry()
        centers.append(np.asarray(center, float))
    return np.asarray(centers, float)


def initialize_membrane(universe, cfg: dict) -> MembraneState:
    spec = cfg.get("analyses", {}).get("membrane", {})
    state = MembraneState(enabled=bool(spec.get("enabled", False)))
    if not state.enabled:
        return state
    state.axis = str(spec.get("axis", "z")).lower()
    state.ai = AXIS_INDEX[state.axis]
    state.apl = bool(spec.get("area_per_lipid", {}).get("enabled", True))
    state.thickness = bool(spec.get("thickness", {}).get("enabled", True))
    state.headgroups = universe.select_atoms(str(spec["headgroup_selection"]))
    if not len(state.headgroups):
        raise ValueError(
            "membrane.headgroup_selection matched zero atoms: "
            f"{spec['headgroup_selection']}"
        )
    centers = [(
        int(res.ix),
        float((state.headgroups & res.atoms).center_of_mass()[state.ai])
    ) for res in state.headgroups.residues]
    midpoint = float(np.median([value for _, value in centers]))
    state.upper_resix = np.asarray([ix for ix, value in centers if value >= midpoint], int)
    state.lower_resix = np.asarray([ix for ix, value in centers if value < midpoint], int)
    state.upper_atoms = state.headgroups[np.isin(state.headgroups.resindices, state.upper_resix)]
    state.lower_atoms = state.headgroups[np.isin(state.headgroups.resindices, state.lower_resix)]
    state.n_upper, state.n_lower = len(state.upper_resix), len(state.lower_resix)
    if not state.n_upper or not state.n_lower:
        raise ValueError("first-frame leaflet assignment produced an empty leaflet")

    msd = spec.get("msd", {}) or {}
    state.msd = bool(msd.get("enabled", False))
    state.msd_include_3d = bool(msd.get("include_3d", False))
    if state.msd:
        if universe.dimensions is None or not np.all(np.asarray(
            universe.dimensions[:3],
            float
        ) > 0):
            raise ValueError("membrane MSD requires valid periodic box dimensions")
        atoms = universe.select_atoms(str(msd.get("selection") or spec["headgroup_selection"]))
        if not len(atoms):
            raise ValueError("membrane.msd.selection matched zero atoms")
        state.msd_residue_atoms = [atoms & residue.atoms for residue in atoms.residues]
        state.msd_resindices = np.asarray([int(residue.ix) for residue in atoms.residues], int)
        state.msd_upper_mask = np.isin(state.msd_resindices, state.upper_resix)
        state.msd_lower_mask = np.isin(state.msd_resindices, state.lower_resix)
        if not np.all(state.msd_upper_mask | state.msd_lower_mask):
            raise ValueError(
                "membrane.msd.selection includes residues outside the leaflet classification"
            )
        state.msd_last_positions = _residue_centers(state.msd_residue_atoms)
        state.msd_displacements = np.zeros_like(state.msd_last_positions)
    return state


def membrane_metrics(universe, state: MembraneState) -> dict[str, float]:
    if not state.enabled:
        return {}
    dims = np.asarray(universe.dimensions[:3], float)
    vectors = np.asarray(getattr(universe.trajectory.ts, "triclinic_dimensions", None))
    other = [i for i in range(3) if i != state.ai]
    area = (
        float(
            np.linalg.norm(np.cross(vectors[other[0]], vectors[other[1]]))
        ) if vectors.shape == (3, 3) else float(dims[other[0]] * dims[other[1]])
    )
    out = {}
    if state.apl:
        out["area_per_lipid_upper_A2"] = area / state.n_upper
        out["area_per_lipid_lower_A2"] = area / state.n_lower
        out["area_per_lipid_A2"] = (
            0.5 * (out["area_per_lipid_upper_A2"] + out["area_per_lipid_lower_A2"])
        )
    if state.thickness:
        delta = abs(
            float(state.upper_atoms.center_of_mass()[state.ai])
            - float(state.lower_atoms.center_of_mass()[state.ai])
        )
        length = float(dims[state.ai])
        if length > 0:
            delta %= length
            delta = min(delta, length - delta)
        out["membrane_thickness_A"] = delta
    if state.msd:
        current = _residue_centers(state.msd_residue_atoms)
        # Accumulate minimum-image step vectors so lipids crossing
        # periodic boundaries do not create artificial MSD jumps.
        state.msd_displacements += minimize_vectors(
            current - state.msd_last_positions,
            universe.dimensions
        )
        state.msd_last_positions = current
        lateral = [i for i in range(3) if i != state.ai]
        lateral_sq = np.sum(state.msd_displacements[:, lateral] ** 2, axis=1)
        out["membrane_msd_lateral_A2"] = float(np.mean(lateral_sq))
        out["membrane_msd_lateral_upper_A2"] = float(np.mean(lateral_sq[state.msd_upper_mask]))
        out["membrane_msd_lateral_lower_A2"] = float(np.mean(lateral_sq[state.msd_lower_mask]))
        if state.msd_include_3d:
            total_sq = np.sum(state.msd_displacements ** 2, axis=1)
            out["membrane_msd_3d_A2"] = float(np.mean(total_sq))
            out["membrane_msd_3d_upper_A2"] = float(np.mean(total_sq[state.msd_upper_mask]))
            out["membrane_msd_3d_lower_A2"] = float(np.mean(total_sq[state.msd_lower_mask]))
    return out
