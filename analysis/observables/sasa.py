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


def _sasa_radii(atoms):
    """Classify static atom identities once; preserve FreeSASA's exact radii."""
    freesasa, structure = _freesasa_structure(atoms)
    if structure.nAtoms() != len(atoms):
        raise ValueError("FreeSASA excluded atoms from the SASA selection")
    return freesasa, [structure.radius(i) for i in range(len(atoms))]


def _sasa_areas(freesasa, atoms, radii, params):
    """Calculate SASA using current coordinates and cached atomic radii."""
    coordinates = atoms.positions.astype(np.float64, copy=False).ravel().tolist()
    result = freesasa.calcCoord(coordinates, radii, params)
    return np.fromiter(
        (result.atomArea(i) for i in range(len(atoms))),
        dtype=float, count=len(atoms)
    )


def _run_sasa(universe, cfg: dict, repout: Path):
    spec, tr = cfg["analyses"]["sasa"], cfg["trajectory"]
    selection = str(spec.get("selection") or cfg["selections"]["protein"])
    atoms = universe.select_atoms(selection)
    if not len(atoms):
        raise ValueError(f"sasa.selection matched zero atoms: {selection}")
    residues = list(atoms.residues)
    residue_lookup = {int(res.ix): i for i, res in enumerate(residues)}
    atom_residue = np.fromiter(
        (residue_lookup[int(ix)] for ix in atoms.resindices),
        dtype=int,
        count=len(atoms)
    )
    freesasa, radii = _sasa_radii(atoms)
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
    rows = []
    frames=universe.trajectory[tr.get("start"):tr.get("stop"):int(tr["stride"])]
    for frame_index, ts in enumerate(track(frames,total=len(frames),desc="SASA",
                enabled=cfg.get("progress",{}).get("enabled",True))):
        areas = _sasa_areas(freesasa, atoms, radii, params)
        residue_areas = (
            np.bincount(
                atom_residue,
                weights=areas,
                minlength=len(residues),
            )
            / 100.0
        )
        rows.extend(
            (
                time_ps(ts, frame_index, cfg),
                int(res.ix),
                _residue_label(res),
                float(area),
            )
            for res, area in zip(residues, residue_areas)
        )
    out = ensure_dir(repout / "profiles" / "sasa")
    atomic_parquet(
        out / "residue_sasa.parquet",
        pd.DataFrame(
            rows,
            columns=[
                "time_ps",
                "topology_residue_index",
                "residue",
                "sasa_nm2",
            ],
        )
    )
