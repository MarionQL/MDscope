"""Frame-wise 1D mass/number density with profile cache and resumable accumulators."""
from __future__ import annotations
from dataclasses import dataclass,field
from pathlib import Path
import numpy as np
import pandas as pd
from .config import DEFAULT_INSTANCE_ID, analysis_groups, analysis_runs, atomic_parquet, ensure_dir, group_label, profile_dir
AXIS_INDEX={"x":0,"y":1,"z":2}
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


def _center_coordinate(atoms, axis_index: int) -> float:
    try:
        center = atoms.center_of_mass()
    except (AttributeError, ValueError, TypeError):
        center = atoms.center_of_geometry()
    if not np.all(np.isfinite(center)):
        center = atoms.center_of_geometry()
    return float(center[axis_index])


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
            _center_coordinate(self.center_atoms, self.ai)
            if self.center_atoms is not None
            else length / 2.0
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
