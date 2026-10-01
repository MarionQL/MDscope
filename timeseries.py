"""
Time-series analyses: shared trajectory pass, protein/membrane metrics,
and 1D density.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis import transformations
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import minimize_vectors

from .config import (
    DEFAULT_INSTANCE_ID,
    analysis_groups,
    analysis_runs,
    atomic_parquet,
    ensure_dir,
    group_label,
    profile_dir,
    time_ps,
)
from .interactions import (
    chain_label,
    initialize_ligand_contacts,
    initialize_salt_bridges,
    ligand_contact_events,
    residue_label,
    salt_bridge_events,
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
    )for res in state.headgroups.residues]
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
        delta = (
            abs(
                float(
                    state.upper_atoms.center_of_mass()[state.ai]
                ) - float(state.lower_atoms.center_of_mass()[state.ai])
            )
        )
        length = float(dims[state.ai])
        out["membrane_thickness_A"] = min(delta, length - delta) if length > 0 else delta
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


@dataclass
class DensityTarget:
    key: str
    label: str
    atoms: object
    mode: str
    counts: np.ndarray
    density_sum: np.ndarray
    frames: int = 0
    frame_density: list[np.ndarray] = field(default_factory=list)


class DensityProfileAccumulator:
    """
    Accumulate mean density while retaining chunked frame profiles for
    block bootstrap statistics.
    """

    def __init__(
        self,
        universe,
        cfg: dict,
        spec: dict | None = None,
        run_id: str = DEFAULT_INSTANCE_ID
    ):
        spec = spec if spec is not None else cfg.get("analyses", {}).get("density_1d", {})
        self.run_id = str(run_id)
        self.enabled = bool(spec.get("enabled", False))
        self.frame_times: list[float] = []
        if not self.enabled:
            return
        self.axis = str(spec.get("axis", "z")).lower()
        if self.axis not in AXIS_INDEX:
            raise ValueError("density_1d.axis must be x, y, or z")
        self.ai = AXIS_INDEX[self.axis]
        self.bin_width = float(spec.get("bin_width_A", 1.0))
        self.center_atoms = (
            universe.select_atoms(
                str(spec["center_selection"])
            ) if spec.get("center_selection") else None
        )
        if self.center_atoms is not None and not len(self.center_atoms):
            raise ValueError(
                "density_1d.center_selection matched zero atoms: "
                f"{spec['center_selection']}"
            )
        length = float(universe.dimensions[self.ai]) if universe.dimensions is not None else 0.0
        if length <= 0:
            raise ValueError("density_1d requires valid periodic box dimensions")
        nbins = max(1, int(np.ceil(length / self.bin_width)))
        self.edges = np.linspace(-length / 2.0, length / 2.0, nbins + 1)
        self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.targets: list[DensityTarget] = []
        for key, raw in (spec.get("selections", {}) or {}).items():
            target = {"selection": raw} if isinstance(raw, str) else raw
            selection = target.get("selection")
            if not selection:
                raise ValueError(f"density_1d.selections.{key} is missing selection")
            atoms = universe.select_atoms(str(selection))
            if not len(atoms):
                raise ValueError(f"density_1d selection {key!r} matched zero atoms: {selection}")
            mode = str(target.get("mode", "number")).lower()
            if mode not in {"number", "mass"}:
                raise ValueError("density_1d mode must be number or mass")
            self.targets.append(DensityTarget(
                str(key),
                str(target.get("label", key)),
                atoms,
                mode,
                np.zeros(nbins),
                np.zeros(nbins)
            ))
            if target.get("groups", False):
                for group_id, group in analysis_groups(cfg).items():
                    subset = atoms.intersection(universe.select_atoms(group["selection"]))
                    if len(subset):
                        self.targets.append(DensityTarget(
                            f"{key}__{group_id}",
                            f"{target.get('label', key)} * {group_label(cfg, group_id)}",
                            subset, mode, np.zeros(nbins), np.zeros(nbins),
                        ))

    def update(self, universe, tps: float | None = None):
        if not self.enabled:
            return
        length = float(universe.dimensions[self.ai])
        dims = np.asarray(universe.dimensions[:3], float)
        if length <= 0 or np.any(dims <= 0):
            return
        origin = (
            float(
                self.center_atoms.center_of_mass()[self.ai]
            ) if self.center_atoms is not None else length / 2.0
        )
        cross_area = float(getattr(universe.trajectory.ts, "volume", np.prod(dims))) / length
        volume = cross_area * np.diff(self.edges)
        if tps is not None:
            self.frame_times.append(float(tps))
        for target in self.targets:
            x = np.asarray(target.atoms.positions[:, self.ai], float) - origin
            x = (x + length / 2.0) % length - length / 2.0
            weights = None if target.mode == "number" else np.asarray(target.atoms.masses, float)
            counts, _ = np.histogram(x, bins=self.edges, weights=weights)
            density = counts / volume
            target.counts += counts
            target.density_sum += density
            target.frames += 1
            if tps is not None:
                # Retain frame profiles only until the next chunk flush;
                # float32 + compression keeps hierarchical profile
                # bootstrap practical.
                target.frame_density.append(np.asarray(density, np.float32))

    def flush_frames(self, repout: Path, part: int):
        if not self.enabled or not self.frame_times:
            return
        out = ensure_dir(profile_dir(repout, "density_1d", self.run_id) / "frames")
        times = np.asarray(self.frame_times, float)
        for target in self.targets:
            if target.frame_density:
                path = out / f"{target.key}__part-{part:06d}.npz"
                tmp = path.with_suffix(".tmp.npz")
                np.savez_compressed(
                    tmp,
                    time_ps=times,
                    density=np.stack(target.frame_density).astype(np.float32),
                    positions_A=self.centers.astype(np.float32)
                )
                tmp.replace(path)
                target.frame_density.clear()
        self.frame_times.clear()

    def clear_frame_parts(self, repout: Path):
        if not self.enabled:
            return
        directory = profile_dir(repout, "density_1d", self.run_id) / "frames"
        if directory.exists():
            for path in directory.glob("*.npz"):
                path.unlink(missing_ok=True)

    def checkpoint_payload(self) -> dict:
        if not self.enabled:
            return {}
        prefix = f"density_{self.run_id}_"
        out = {prefix + "edges": self.edges}
        for i, target in enumerate(self.targets):
            out[f"{prefix}{i}_counts"] = target.counts
            out[f"{prefix}{i}_sum"] = target.density_sum
            out[f"{prefix}{i}_frames"] = np.asarray(target.frames)
        return out

    def restore_checkpoint(self, checkpoint) -> None:
        if not self.enabled:
            return
        prefix = f"density_{self.run_id}_"
        if prefix + "edges" not in checkpoint or not np.allclose(
            np.asarray(checkpoint[prefix + "edges"], float),
            self.edges
        ):
            raise ValueError(
                f"density_1d[{self.run_id}] checkpoint bins do not match "
                "configuration"
            )
        for i, target in enumerate(self.targets):
            for suffix in ("counts", "sum", "frames"):
                if f"{prefix}{i}_{suffix}" not in checkpoint:
                    raise ValueError(
                        f"density_1d[{self.run_id}] checkpoint target set does "
                        "not match configuration"
                    )
            target.counts = np.asarray(checkpoint[f"{prefix}{i}_counts"], float)
            target.density_sum = np.asarray(checkpoint[f"{prefix}{i}_sum"], float)
            target.frames = int(checkpoint[f"{prefix}{i}_frames"])

    def finalize(self, repout: Path):
        if not self.enabled:
            return
        out = ensure_dir(profile_dir(repout, "density_1d", self.run_id))
        for target in self.targets:
            atomic_parquet(out / f"{target.key}.parquet", pd.DataFrame({
                f"{self.axis}_A": self.centers,
                "density": target.density_sum / max(1, target.frames),
                "analysis_instance": self.run_id,
                "selection_key": target.key,
                "selection_label": target.label,
                "mode": target.mode,
                "density_unit": "1/A^3" if target.mode == "number" else "amu/A^3",
                "n_frames": target.frames,
            }))


def temporal_density_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    """
    Load chunked per-frame density profiles for hierarchical block
    bootstrap statistics.
    """
    out = {}
    for run in analysis_runs(cfg, "density_1d"):
        directory = profile_dir(repout, "density_1d", run.run_id) / "frames"
        if not directory.exists():
            continue
        keys = sorted({path.name.split(
            "__part-",
            1
        )[0] for path in directory.glob("*__part-*.npz")})
        for key in keys:
            parts = [np.load(path) for path in sorted(directory.glob(f"{key}__part-*.npz"))]
            if not parts:
                continue
            times = np.concatenate([part["time_ps"] for part in parts]).astype(float)
            values = np.concatenate([part["density"] for part in parts]).astype(float)
            positions = np.asarray(parts[0]["positions_A"], float)
            out[f"density_1d__{run.run_id}__{key}"] = {
                "times": times,
                "values": values,
                "positions": positions,
                "labels": np.asarray([f"{x:.2f}" for x in positions], object),
                "reducer": "mean",
            }
            for part in parts:
                part.close()
    return out


# Shared trajectory pass and protein-derived time series.
class Welford:
    """
    Online coordinate variance used for RMSF without storing every frame
    in memory.
    """
    def __init__(self, shape):
        self.n = 0
        self.mean = np.zeros(shape, float)
        self.M2 = np.zeros(shape, float)

    def update(self, values):
        self.n += 1
        delta = values - self.mean
        self.mean += delta / self.n
        self.M2 += delta * (values - self.mean)

    def rmsf(self):
        if self.n < 2:
            return np.full(self.mean.shape[0], np.nan)
        return np.sqrt((self.M2 / (self.n - 1)).sum(axis=1))


def _valid_box(universe):
    dims = np.asarray(
        universe.dimensions,
        float
    ) if universe.dimensions is not None else np.array([])
    return universe.dimensions if len(dims) >= 3 and np.all(dims[:3] > 0) else None


def _install_alignment(mobile, reference, selection: str):
    moving = mobile.select_atoms(selection)
    fixed = reference.select_atoms(selection)
    if not len(moving) or len(moving) != len(fixed):
        raise ValueError(
            f"Invalid alignment selection {selection!r}: mobile={len(moving)}, "
            f"reference={len(fixed)} atoms"
        )
    mobile.trajectory.add_transformations(
        transformations.center_in_box(moving, center="mass"),
        transformations.fit_rot_trans(moving, fixed),
    )
    return moving


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
    out ={
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
        out["groups"][group_id] ={
            "group_atoms": len(atoms),
            "feature_atoms": len(measured),
            "alignment_atoms": len(align_atoms)
        }
    return out


def _checkpoint_payload(path, cfg_hash, processed, next_part, rmsf, groups, densities, membrane):
    payload = {
        "config_hash": np.asarray(cfg_hash), "processed": np.asarray(processed),
        "next_part": np.asarray(next_part),
        "w_n": np.asarray(rmsf.n), "w_mean": rmsf.mean, "w_M2": rmsf.M2,
    }
    for group_id, group in groups.items():
        safe = group_id.replace("/", "_")
        payload[f"g_{safe}_n"] = np.asarray(group["rmsf"].n)
        payload[f"g_{safe}_mean"] = group["rmsf"].mean
        payload[f"g_{safe}_M2"] = group["rmsf"].M2
    for density in densities:
        payload.update(density.checkpoint_payload())
    if membrane.enabled and membrane.msd:
        payload.update({
            "membrane_msd_last_positions": membrane.msd_last_positions,
            "membrane_msd_displacements": membrane.msd_displacements,
            "membrane_msd_resindices": membrane.msd_resindices,
        })
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **payload)
    tmp.replace(path)


def _restore_checkpoint(checkpoint, cfg_hash, rmsf, groups, densities, membrane):
    cp = np.load(checkpoint, allow_pickle=False)
    if str(cp["config_hash"].item()) != cfg_hash or cp["w_mean"].shape != rmsf.mean.shape:
        raise ValueError("checkpoint does not match current run")
    rmsf.n, rmsf.mean, rmsf.M2 = (
        int(cp["w_n"]), cp["w_mean"].astype(float), cp["w_M2"].astype(float)
    )
    for group_id, group in groups.items():
        safe = group_id.replace("/", "_")
        group["rmsf"].n = int(cp[f"g_{safe}_n"])
        group["rmsf"].mean = cp[f"g_{safe}_mean"].astype(float)
        group["rmsf"].M2 = cp[f"g_{safe}_M2"].astype(float)
    for density in densities:
        density.restore_checkpoint(cp)
    if membrane.enabled and membrane.msd:
        if not np.array_equal(
            np.asarray(cp["membrane_msd_resindices"], int),
            membrane.msd_resindices
        ):
            raise ValueError("membrane MSD checkpoint residue selection changed")
        membrane.msd_last_positions = np.asarray(cp["membrane_msd_last_positions"], float)
        membrane.msd_displacements = np.asarray(cp["membrane_msd_displacements"], float)
    return int(cp["processed"]), int(cp["next_part"])


def run_streaming_pass(rep, cfg: dict, outroot: Path, cfg_hash: str):
    # Protein, membrane, density, salt-bridge, and ligand-contact
    # calculations share one trajectory pass.
    # Reading a large XTC once is substantially cheaper than reopening
    # it independently for each simple metric.
    repout = ensure_dir(outroot / rep.system / rep.label)
    raw, cache = ensure_dir(repout / "raw"), ensure_dir(repout / "cache")
    checkpoint = cache / "stream_checkpoint.npz"
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    selections = cfg["selections"]
    protein = universe.select_atoms(selections["protein"])
    if not len(protein):
        raise ValueError(f"Protein selection matched zero atoms: {selections['protein']}")
    alignment_spec = cfg["analyses"].get("alignment", {})
    alignment_enabled = bool(alignment_spec.get("enabled", True))
    alignment_selection = (
        str(
            alignment_spec.get("selection") or selections.get("alignment") or "protein and name CA"
        )
    )
    features = universe.select_atoms(selections.get("features") or alignment_selection)
    if not len(features):
        raise ValueError("Feature selection matched zero atoms")

    fit = mda.Universe(str(rep.topology), str(rep.trajectory))
    reference = mda.Universe(str(rep.topology), str(rep.trajectory))
    reference.trajectory[0]
    fit_alignment = _install_alignment(
        fit,
        reference,
        alignment_selection
    ) if alignment_enabled else fit.atoms[[]]
    features_fit = fit.atoms[features.indices]
    reference_features = reference.atoms[features.indices].positions.copy()

    ligand = universe.select_atoms(selections["ligand"]) if selections.get("ligand") else None
    ligand_fit = fit.atoms[ligand.indices] if ligand is not None and len(ligand) else None
    reference_ligand = (
        reference.atoms[ligand.indices].positions.copy() if ligand_fit is not None else None
    )
    salt_state = initialize_salt_bridges(universe, cfg)
    ligand_state = initialize_ligand_contacts(universe, protein, ligand, cfg)
    atomic_parquet(
        raw / "protein_chains.parquet",
        pd.DataFrame({"chain": sorted({chain_label(atom) for atom in protein})})
    )

    rg_spec = cfg["analyses"].get("radius_of_gyration", {})
    rg_atoms = universe.select_atoms(str(rg_spec.get("selection") or selections["protein"]))
    if rg_spec.get("enabled", True) and not len(rg_atoms):
        raise ValueError("radius_of_gyration.selection matched zero atoms")

    densities = [DensityProfileAccumulator(
        universe,
        cfg,
        run.spec,
        run.run_id
    )for run in analysis_runs(cfg, "density_1d")]
    tr = cfg["trajectory"]
    start, stop, stride, chunk = (
        tr.get("start"), tr.get("stop"), int(tr["stride"]), int(tr["chunk_frames"])
    )
    if start is not None:
        universe.trajectory[int(start)]
    membrane = initialize_membrane(universe, cfg)
    rmsf_acc = Welford((len(features), 3))

    groups = {}
    for group_id, group_spec in analysis_groups(cfg).items():
        atoms = universe.select_atoms(group_spec["selection"])
        measured = features.intersection(atoms)
        if not len(atoms) or not len(measured):
            raise ValueError(f"Analysis group {group_id!r} resolved to zero atoms/features")
        align_selection = group_spec.get("alignment_selection") if alignment_enabled else None
        if align_selection:
            group_fit = mda.Universe(str(rep.topology), str(rep.trajectory))
            group_reference = mda.Universe(str(rep.topology), str(rep.trajectory))
            group_reference.trajectory[0]
            _install_alignment(group_fit, group_reference, align_selection)
            measure_fit = group_fit.atoms[measured.indices]
            measure_ref = group_reference.atoms[measured.indices].positions.copy()
        else:
            group_fit = None
            measure_fit = fit.atoms[measured.indices]
            measure_ref = reference.atoms[measured.indices].positions.copy()
        groups[group_id] = {
            "atoms": atoms, "measure": measured, "measure_fit": measure_fit, "fit": group_fit,
            "reference": measure_ref, "rmsf": Welford((
                len(measured),
                3
            )), "rg": rg_atoms.intersection(atoms),
        }

    processed = part = 0
    if checkpoint.exists():
        try:
            processed, part = _restore_checkpoint(
                checkpoint,
                cfg_hash,
                rmsf_acc,
                groups,
                densities,
                membrane
            )
        except Exception:
            for path in raw.glob("*-part-*.parquet"):
                path.unlink(missing_ok=True)
            for path in cache.glob("aligned-part-*.npz"):
                path.unlink(missing_ok=True)
            checkpoint.unlink(missing_ok=True)
            rmsf_acc = Welford((len(features), 3))
            for group in groups.values():
                group["rmsf"] = Welford((len(group["measure"]), 3))
            densities = [DensityProfileAccumulator(
                universe,
                cfg,
                run.spec,
                run.run_id
            )for run in analysis_runs(cfg, "density_1d")]
            for density in densities:
                density.clear_frame_parts(repout)
            if start is not None:
                universe.trajectory[int(start)]
            membrane = initialize_membrane(universe, cfg)
            processed = part = 0

    first = (0 if start is None else int(start)) + processed * stride
    rows, group_rows, salt_rows, ligand_rows, coords, coord_times = [], [], [], [], [], []
    for ts in universe.trajectory[slice(first, stop, stride)]:
        tps = time_ps(ts, processed, cfg)
        for density in densities:
            density.update(universe, tps)
        membrane_values = membrane_metrics(universe, membrane)
        fit.trajectory[ts.frame]
        rmsd_value = (
            float(
                rms.rmsd(
                    features_fit.positions, reference_features, center=False, superposition=False
                )
            ) if cfg["analyses"]["rmsd"].get("enabled", True) else np.nan
        )
        rg_value = float(rg_atoms.radius_of_gyration()) if rg_spec.get("enabled", True) else np.nan
        ligand_rmsd = float(rms.rmsd(
            ligand_fit.positions,
            reference_ligand,
            center=False,
            superposition=False
        )) if reference_ligand is not None else np.nan
        rmsf_acc.update(features_fit.positions.copy())

        for group_id, group in groups.items():
            if group["fit"] is not None:
                group["fit"].trajectory[ts.frame]
            positions = group["measure_fit"].positions.copy()
            group["rmsf"].update(positions)
            group_rows.append({
                "time_ps": tps, "group": group_id,
                "rmsd_A": float(rms.rmsd(
                    positions,
                    group["reference"],
                    center=False,
                    superposition=False
                )),
                "radius_of_gyration_A": float(
                    group["rg"].radius_of_gyration()
                ) if len(group["rg"]) else np.nan,
            })

        salt_rows.extend(salt_bridge_events(salt_state, tps, _valid_box(universe)))
        ligand_min, ligand_count, new_contacts = ligand_contact_events(
            ligand_state,
            tps,
            _valid_box(universe)
        )
        ligand_rows.extend(new_contacts)
        row = {"time_ps": tps, "rmsd_A": rmsd_value, "radius_of_gyration_A": rg_value}
        if ligand_fit is not None:
            row["ligand_rmsd_A"] = ligand_rmsd
        if cfg["analyses"]["ligand_contacts"].get("enabled", False):
            row["ligand_min_distance_A"] = ligand_min
            row["ligand_contact_residue_count"] = ligand_count
        row.update(membrane_values)
        rows.append(row)
        if cfg["cache"].get("aligned_coordinates", True):
            coords.append(features_fit.positions.astype(
                cfg["cache"].get("dtype", "float32"),
                copy=True
            ))
            coord_times.append(tps)
        processed += 1
        if len(rows) >= chunk:
            _flush(
                raw,
                cache,
                rows,
                group_rows,
                salt_rows,
                ligand_rows,
                coords,
                coord_times,
                part,
                cfg
            )
            for density in densities:
                density.flush_frames(repout, part)
            rows, group_rows, salt_rows, ligand_rows, coords, coord_times = [], [], [], [], [], []
            part += 1
            _checkpoint_payload(
                checkpoint,
                cfg_hash,
                processed,
                part,
                rmsf_acc,
                groups,
                densities,
                membrane
            )

    if rows or group_rows or salt_rows or ligand_rows or coords:
        _flush(
            raw,
            cache,
            rows,
            group_rows,
            salt_rows,
            ligand_rows,
            coords,
            coord_times,
            part,
            cfg
        )
        for density in densities:
            density.flush_frames(repout, part)
        part += 1
        _checkpoint_payload(
            checkpoint,
            cfg_hash,
            processed,
            part,
            rmsf_acc,
            groups,
            densities,
            membrane
        )

    if cfg["cache"].get("aligned_coordinates", True):
        atomic_parquet(
            cache / "feature_atoms.parquet",
            pd.DataFrame({
                "atom_index": features.indices,
                "residue": [residue_label(atom) for atom in features],
                "atom_name": features.names
            })
        )

    rmsf_values = rmsf_acc.rmsf()
    atomic_parquet(raw / "rmsf.parquet", pd.DataFrame({
        "atom_index": features.indices, "residue": [residue_label(atom) for atom in features],
        "atom_name": features.names,
        "rmsf_A": rmsf_values,
        "equivalent_bfactor_A2": (8.0 * np.pi ** 2 / 3.0) * rmsf_values ** 2,
    }))
    group_rmsf = []
    group_residues = []
    for group_id, group in groups.items():
        group_rmsf.extend(
            (group_id, int(atom.index), residue_label(atom), atom.name, float(value)) for atom,
            value in zip(group["measure"], group["rmsf"].rmsf())
        )
        group_residues.extend((
            group_id,
            int(res.ix),
            str(res.resname),
            int(res.resid)
        )for res in group["atoms"].residues)
    if group_rmsf:
        atomic_parquet(
            raw / "group_rmsf.parquet",
            pd.DataFrame(
                group_rmsf,
                columns=["group", "atom_index", "residue", "atom_name", "rmsf_A"]
            )
        )
        atomic_parquet(
            raw / "group_residues.parquet",
            pd.DataFrame(
                group_residues,
                columns=["group", "topology_residue_index", "resname", "resid"]
            ).drop_duplicates()
        )
    for density in densities:
        density.finalize(repout)

    try:
        bfactor = (8.0 * np.pi ** 2 / 3.0) * rmsf_values ** 2
        if not hasattr(universe.atoms, "tempfactors"):
            universe.add_TopologyAttr("tempfactors", np.zeros(len(universe.atoms)))
        universe.atoms.tempfactors = 0.0
        features.tempfactors = bfactor
        protein.write(str(raw / "rmsf_bfactor.pdb"))
    except Exception as exc:
        print(f"Warning: could not write RMSF B-factor PDB: {exc}")
    return{
        "frames": processed,
        "topology_atoms": len(universe.atoms),
        "protein_atoms": len(protein),
        "feature_atoms": len(features),
        "alignment_atoms": len(fit_alignment)
    }


def _flush(raw, cache, rows, group_rows, salt_rows, ligand_rows, coords, times, part, cfg):
    # Chunked writes cap RAM use and make completed trajectory segments
    # recoverable after wall-time/job failures.
    tag = f"part-{part:06d}"
    if rows:
        atomic_parquet(raw / f"timeseries-{tag}.parquet", pd.DataFrame(rows))
    if group_rows:
        atomic_parquet(raw / f"group_timeseries-{tag}.parquet", pd.DataFrame(group_rows))
    if salt_rows:
        columns = [
            "time_ps",
            "residue_1",
            "residue_2",
            "chain_1",
            "chain_2",
            "scope",
            "interaction_class",
            "min_distance_A",
        ]
        atomic_parquet(
            raw / f"salt_bridges-{tag}.parquet",
            pd.DataFrame(salt_rows, columns=columns),
        )
    if ligand_rows and cfg["analyses"]["ligand_contacts"].get("save_atom_events", True):
        columns = [
            "time_ps",
            "ligand_atom",
            "protein_atom",
            "protein_residue",
            "distance_A",
        ]
        atomic_parquet(
            raw / f"ligand_contacts-{tag}.parquet",
            pd.DataFrame(ligand_rows, columns=columns),
        )
    if coords:
        tmp = cache / f"aligned-{tag}.tmp.npz"
        # float32 aligned coordinates are accurate enough for RMSF while
        # roughly halving this cache's storage footprint.
        np.savez_compressed(tmp, xyz=np.stack(coords), time_ps=np.asarray(times, float))
        tmp.replace(cache / f"aligned-{tag}.npz")


def temporal_rmsf_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    """
    Load aligned coordinates needed for block-bootstrap RMSF profiles.
    """
    cache = Path(repout) / "cache"
    atom_file = cache / "feature_atoms.parquet"
    parts = sorted(cache.glob("aligned-part-*.npz"))
    if not atom_file.exists() or not parts:
        return {}
    atoms = pd.read_parquet(atom_file)
    loaded = [np.load(path) for path in parts]
    times = np.concatenate([part["time_ps"] for part in loaded]).astype(float)
    xyz = np.concatenate([part["xyz"] for part in loaded]).astype(np.float32)
    for part in loaded:
        part.close()
    out = {
        "rmsf": {
            "times": times,
            "coords": xyz,
            "atom_indices": atoms["atom_index"].to_numpy(int),
            "atom_residues": atoms["residue"].astype(str).to_numpy(),
            "reducer": "rmsf",
        }
    }
    group_file = Path(repout) / "raw" / "group_rmsf.parquet"
    if group_file.exists():
        group_atoms = pd.read_parquet(group_file)
        index_lookup = {int(atom): i for i, atom in enumerate(atoms["atom_index"].to_numpy(int))}
        for group_id, group in group_atoms.groupby("group", sort=False):
            columns =[
                index_lookup[
                    int(atom)
                ]for atom in group["atom_index"].to_numpy(int) if int(atom) in index_lookup
            ]
            if columns:
                out[f"rmsf__group_{group_id}"] = {
                    "times": times,
                    "coords": xyz[:, columns, :],
                    "atom_indices": atoms.iloc[columns]["atom_index"].to_numpy(int),
                    "atom_residues": atoms.iloc[columns]["residue"].astype(str).to_numpy(),
                    "reducer": "rmsf",
                }
    return out
