"""Membrane APL, thickness and the original minimum-image cumulative MSD implementation.

One shared frame pass updates this state; no duplicate diffusion/MSD trajectory reader.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .observables.msd import initialize_msd,msd_metrics
from .core.leaflet_partition import median_leaflet_mask
AXIS_INDEX = {'x':0,'y':1,'z':2}


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
    msd_remove_drift: bool = False
    msd_atoms: object = None
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
    upper_mask = median_leaflet_mask([value for _, value in centers])
    state.upper_resix = np.asarray([ix for (ix, _), upper in zip(centers, upper_mask) if upper], int)
    state.lower_resix = np.asarray([ix for (ix, _), upper in zip(centers, upper_mask) if not upper], int)
    state.upper_atoms = state.headgroups[np.isin(state.headgroups.resindices, state.upper_resix)]
    state.lower_atoms = state.headgroups[np.isin(state.headgroups.resindices, state.lower_resix)]
    state.n_upper, state.n_lower = len(state.upper_resix), len(state.lower_resix)
    if not state.n_upper or not state.n_lower:
        raise ValueError("first-frame leaflet assignment produced an empty leaflet")

    initialize_msd(universe,state,spec)
    return state


def membrane_metrics(universe, state: MembraneState, msd_universe=None) -> dict[str, float]:
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
        out.update(msd_metrics(msd_universe if msd_universe is not None else universe, state))
    return out
