"""FreeSASA residue solvent-accessible surface area."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from ..core.profiling import _residue_label
from ..pipeline.progress import track
from ..config import ensure_dir,atomic_parquet,time_ps
def _freesasa_structure(atoms):
    try:
        import freesasa
    except ImportError as exc:
        raise RuntimeError("SASA requires FreeSASA (`pip install freesasa`).") from exc
    structure = freesasa.Structure()
    for atom in atoms:
        x, y, z = map(float, atom.position)
        structure.addAtom(
            str(atom.name),
            str(atom.resname),
            int(atom.resid),
            str(atom.segid),
            x,
            y,
            z
        )
    return freesasa, structure


# Explicit van der Waals radii (Angstrom) for heterogeneous membrane atoms.
# This is a geometric approximation, not a force-field SASA parameterization.
ELEMENT_RADII_A = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
                   "P": 1.80, "S": 1.80, "F": 1.47, "CL": 1.75,
                   "BR": 1.85, "I": 1.98}


def sasa_instances(universe, cfg):
    """Build independent target/occluder SASA selections from one Universe."""
    from ..config import analysis_runs
    try:
        import freesasa
    except ImportError as exc:
        raise RuntimeError("SASA requires freesasa") from exc
    instances = {}
    for run in analysis_runs(cfg, "sasa"):
        spec = run.spec
        target_sel = str(spec.get("target_selection") or spec.get("selection")
                         or cfg["selections"]["protein"])
        occ_sel = str(spec.get("occluder_selection") or target_sel)
        target = universe.select_atoms(target_sel)
        occluders = universe.select_atoms(occ_sel)
        if not len(target) or not len(occluders):
            raise ValueError(f"SASA run {run.run_id}: target/occluder selection is empty")
        target_idx = np.flatnonzero(np.isin(occluders.indices, target.indices))
        if len(target_idx) != len(target):
            raise ValueError(f"SASA run {run.run_id}: occluders must include every target atom")
        residues = list(target.residues)
        residue_ids = np.asarray([r.ix for r in residues], dtype=int)
        mapping = np.searchsorted(residue_ids, target.resindices)
        # Indices from MDAnalysis residues are sorted by topology index.
        if not np.array_equal(residue_ids[mapping], target.resindices):
            raise ValueError("SASA target residue mapping failed")
        radii_spec = spec.get("element_radii_A", {})
        radii_model = {**ELEMENT_RADII_A, **{str(k).upper(): float(v)
                                         for k, v in radii_spec.items()}}
        elements = getattr(occluders, "elements", None)
        if elements is None:
            raise ValueError("Environment SASA requires topology atom elements")
        symbols = [str(e).upper() for e in elements]
        missing = sorted(set(symbols) - radii_model.keys())
        if missing:
            raise ValueError(f"SASA run {run.run_id}: no radii for elements {missing}")
        radii = [float(radii_model[e]) for e in symbols]
        params = freesasa.Parameters()
        algorithm = str(spec.get("algorithm", "shrake_rupley")).lower()
        if algorithm in {"shrake_rupley", "shrake-rupley", "sr"}:
            params.setAlgorithm(freesasa.ShrakeRupley)
            if spec.get("n_points") is not None:
                params.setNPoints(int(spec["n_points"]))
        elif algorithm in {"lee_richards", "lee-richards", "lr"}:
            params.setAlgorithm(freesasa.LeeRichards)
            if spec.get("n_slices") is not None:
                params.setNSlices(int(spec["n_slices"]))
        else:
            raise ValueError("sasa.algorithm must be shrake_rupley or lee_richards")
        if spec.get("probe_radius_A") is not None:
            params.setProbeRadius(float(spec["probe_radius_A"]))
        if spec.get("n_threads") is not None:
            params.setNThreads(int(spec["n_threads"]))
        instances[run.run_id] = {
            "occluders": occluders, "radii": radii, "target_idx": target_idx,
            "residues": residues, "mapping": mapping, "params": params,
            "library": freesasa,
        }
    return instances


def sasa_residue_rows(state, time):
    """Return per-residue protein SASA in nm², including occlusion by other atoms."""
    atoms = state["occluders"]
    result = state["library"].calcCoord(
        atoms.positions.astype(np.float64, copy=False).ravel().tolist(),
        state["radii"], state["params"],
    )
    areas = np.fromiter(
        (result.atomArea(int(i)) for i in state["target_idx"]),
        dtype=float, count=len(state["target_idx"]),
    )
    values = np.bincount(
        state["mapping"], weights=areas, minlength=len(state["residues"])
    ) / 100.0
    return [(time, int(res.ix), _residue_label(res), float(area))
            for res, area in zip(state["residues"], values)]


def _run_sasa(universe, cfg: dict, repout: Path):
    """Standalone SASA path using the same implementation as the fused observer."""
    states = sasa_instances(universe, cfg)
    tr = cfg["trajectory"]
    rows = {key: [] for key in states}
    frames = universe.trajectory[tr.get("start"):tr.get("stop"):int(tr["stride"])]
    for frame_index, ts in enumerate(track(
            frames, total=len(frames), desc="SASA",
            enabled=cfg.get("progress", {}).get("enabled", True))):
        t = time_ps(ts, frame_index, cfg)
        for key, state in states.items():
            rows[key].extend(sasa_residue_rows(state, t))
    for key, values in rows.items():
        output_dir = (Path(repout) / "profiles" / "sasa" /
                      key if key != "default" else Path(repout) / "profiles" / "sasa")
        atomic_parquet(ensure_dir(output_dir) / "residue_sasa.parquet",
                       pd.DataFrame(values, columns=[
                           "time_ps", "topology_residue_index", "residue", "sasa_nm2"
                       ]))
