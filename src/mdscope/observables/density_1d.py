"""Frame-wise 1D mass/number density with profile cache and resumable accumulators."""
from __future__ import annotations
from dataclasses import dataclass,field
import json
from pathlib import Path
import numpy as np
import pandas as pd
from ..config import DEFAULT_INSTANCE_ID, analysis_groups, analysis_runs, atomic_json, atomic_parquet, ensure_dir, group_label, profile_dir
from ..core.lipid_geometry import bilayer_midplane, planar_box

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


def _membrane_reference(spec):
    return {
        "coordinate_reference": "phosphate_leaflet_midplane",
        "selection": spec["membrane_center_selection"],
        "axis": str(spec.get("axis", "z")).lower(),
        "bin_width_A": float(spec.get("bin_width_A", 1.0)),
        "half_width_A": float(spec["half_width_A"]),
    }


def missing_membrane_density_runs(repout, cfg):
    """Find density runs that need a trajectory pass without reusing box bins."""
    missing = []
    for run in analysis_runs(cfg, "density_1d"):
        if not run.spec.get("membrane_center_selection"):
            continue
        directory = profile_dir(repout, "density_1d", run.run_id)
        marker = directory / "reference.json"
        if not marker.is_file():
            missing.append(run)
            continue
        saved = json.loads(marker.read_text())
        if (saved["reference"] != _membrane_reference(run.spec)
                or saved["selections"] != run.spec["selections"]):
            raise ValueError(f"{directory}: density reference differs from configuration; use a new run ID")
        frame_directory = saved["frame_directory"]
        if frame_directory not in {"frames", "replayed_frames"}:
            raise ValueError(f"{directory}: invalid membrane density frame directory")
        if not saved["targets"] or any(
            not (directory / f"{target}.parquet").is_file()
            or not any((directory / frame_directory).glob(f"{target}__part-*.npz"))
            for target in saved["targets"]
        ):
            raise ValueError(f"{directory}: incomplete membrane-relative density output")
    return missing


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
        self.membrane_reference = None
        if spec.get("membrane_center_selection"):
            self.membrane_reference = _membrane_reference(spec)
            self.membrane_selections = spec["selections"]
            self.phosphates = universe.select_atoms(spec["membrane_center_selection"])
            if len(self.phosphates) < 4 or len(self.phosphates) != len(self.phosphates.residues):
                raise ValueError("Membrane density requires one anchor per lipid and at least four anchors")
        length = float(universe.dimensions[self.ai]) if universe.dimensions is not None else 0.0
        if length <= 0:
            raise ValueError("density_1d requires valid periodic box dimensions")
        if self.membrane_reference is None:
            nbins = max(1, int(np.ceil(length / self.bin_width)))
            self.edges = np.linspace(-length / 2.0, length / 2.0, nbins + 1)
        else:
            half_bins = round(self.membrane_reference["half_width_A"] / self.bin_width)
            self.edges = self.bin_width * np.arange(-half_bins, half_bins + 1)
            nbins = len(self.edges) - 1
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
        if self.membrane_reference is not None:
            planar_box(universe.dimensions, self.ai)
            if self.edges[-1] > length / 2:
                raise ValueError("Membrane density half_width_A exceeds half the current box length")
            origin = bilayer_midplane(self.phosphates.positions[:, self.ai], length)
        else:
            origin = (_center_coordinate(self.center_atoms, self.ai)
                      if self.center_atoms is not None else length / 2.0)
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
        if self.membrane_reference is not None:
            out[prefix + "reference"] = np.asarray(json.dumps(self.membrane_reference, sort_keys=True))
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
        if self.membrane_reference is not None:
            key = prefix + "reference"
            if key not in checkpoint or json.loads(str(checkpoint[key].item())) != self.membrane_reference:
                raise ValueError(f"density_1d[{self.run_id}] checkpoint reference differs")
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

    def finalize(self, repout: Path, *, replay=False):
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
        if self.membrane_reference is not None:
            frame_directory = "replayed_frames" if replay else "frames"
            if replay:
                (out / "frames").rename(out / frame_directory)
            atomic_json(out / "reference.json", {
                "frame_directory": frame_directory,
                "reference": self.membrane_reference,
                "selections": self.membrane_selections,
                "frames": self.targets[0].frames,
                "targets": [target.key for target in self.targets],
            })


def temporal_density_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    """
    Load chunked per-frame density profiles for hierarchical block
    bootstrap statistics.
    """
    out = {}
    missing = missing_membrane_density_runs(repout, cfg)
    if missing:
        raise ValueError(f"{repout}: missing membrane-relative density; run MDscope with this config")
    for run in analysis_runs(cfg, "density_1d"):
        directory = profile_dir(repout, "density_1d", run.run_id)
        relative = run.spec.get("membrane_center_selection")
        recorded = json.loads((directory / "reference.json").read_text()) if relative else None
        directory = directory / (recorded["frame_directory"] if recorded else "frames")
        if not directory.exists():
            continue
        keys = sorted({path.name.split(
            "__part-",
            1
        )[0] for path in directory.glob("*__part-*.npz")})
        if recorded is not None and set(keys) != set(recorded["targets"]):
            raise ValueError(f"{directory}: membrane density frame targets differ from saved reference")
        common_times = None
        for key in keys:
            chunks = sorted(directory.glob(f"{key}__part-*.npz"))
            time_chunks, density_chunks = [], []
            positions = None
            for path in chunks:
                with np.load(path) as part:
                    current = np.asarray(part["positions_A"], float)
                    if positions is None:
                        positions = current
                    elif not np.array_equal(positions, current):
                        raise ValueError(f"{path}: density grid changes between chunks")
                    time_chunks.append(part["time_ps"])
                    density_chunks.append(part["density"])
            times = np.concatenate(time_chunks).astype(float)
            values = np.concatenate(density_chunks).astype(float)
            if recorded is not None:
                if (len(times) != recorded["frames"] or not np.isfinite(times).all()
                        or np.any(np.diff(times) <= 0)
                        or (common_times is not None and not np.array_equal(common_times, times))):
                    raise ValueError(f"{directory}: incomplete or inconsistent membrane density frame grid")
                common_times = times
            out[f"density_1d__{run.run_id}__{key}"] = {
                "density_bin_width_A": float(run.spec.get("bin_width_A", 1.0)),
                "times": times,
                "values": values,
                "positions": positions,
                "labels": np.asarray([f"{x:.2f}" for x in positions], object),
                "reducer": "mean",
            }
    return out


def align_density_profiles(specs: dict, profile: str) -> dict:
    """Rebin offset uniform density slabs on a shared zero-anchored axis.

    Legacy frame caches contain slab centers. The accumulator uses uniform
    edges, so their physical widths can be recovered from the center spacing.
    Treat each saved density as constant within its slab and retain only
    target bins fully supported by every source grid.
    """
    axes = [np.asarray(spec["positions"], dtype=float) for spec in specs.values()]
    if len(axes) < 2 or all(np.array_equal(axis, axes[0]) for axis in axes[1:]):
        return specs

    source_edges = []
    for centers in axes:
        spacing = np.diff(centers)
        if not len(spacing) or not np.isfinite(centers).all():
            raise ValueError(f"{profile}: cannot recover density edges from cached centers")
        step = float(np.median(spacing))
        if step <= 0 or not np.allclose(spacing, step, rtol=1e-4, atol=2e-5):
            raise ValueError(f"{profile}: cached density grid must be increasing and uniform")
        source_edges.append(np.r_[
            centers[0] - step / 2, (centers[:-1] + centers[1:]) / 2,
            centers[-1] + step / 2,
        ])

    widths = {spec["density_bin_width_A"] for spec in specs.values()}
    width = next(iter(widths))
    if len(widths) != 1 or not np.isfinite(width) or width <= 0:
        raise ValueError(f"{profile}: density profiles require one positive target bin width")
    lower = max(edges[0] for edges in source_edges)
    upper = min(edges[-1] for edges in source_edges)
    first = int(np.ceil(lower / width))
    last = int(np.floor(upper / width))
    target_edges = width * np.arange(first, last + 1)
    target_edges = target_edges[(target_edges >= lower) & (target_edges <= upper)]
    if len(target_edges) < 2:
        raise ValueError(f"{profile}: density grids have no fully shared {width:g} Å bin")
    positions = (target_edges[:-1] + target_edges[1:]) / 2
    labels = positions.astype(str)

    aligned = {}
    for (key, spec), edges in zip(specs.items(), source_edges):
        values = np.asarray(spec["values"], dtype=float)
        if values.ndim != 2 or values.shape[1] != len(edges) - 1:
            raise ValueError(f"{profile}: frame density matrix does not match cached grid")
        rebinned = np.empty((len(values), len(positions)), dtype=float)
        for i, (left, right) in enumerate(zip(target_edges[:-1], target_edges[1:])):
            start = np.searchsorted(edges, left, side="right") - 1
            stop = np.searchsorted(edges, right, side="left")
            overlap = np.minimum(edges[start + 1:stop + 1], right) - np.maximum(
                edges[start:stop], left
            )
            rebinned[:, i] = np.sum(values[:, start:stop] * overlap, axis=1) / width
        aligned[key] = {**spec, "values": rebinned, "positions": positions, "labels": labels}
    return aligned
