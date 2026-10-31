"""Residue secondary structure assignments from DSSP."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from MDAnalysis.analysis.dssp import DSSP
from ..profile_common import _residue_label
from ..config import ensure_dir,atomic_parquet
def _run_dssp(universe, cfg: dict, repout: Path):
    spec, tr = cfg["analyses"]["dssp"], cfg["trajectory"]
    selection = str(spec.get("selection") or cfg["selections"]["protein"])
    atoms = universe.select_atoms(selection)
    if not len(atoms):
        raise ValueError(f"dssp.selection matched zero atoms: {selection}")
    analysis = DSSP(
        atoms,
        guess_hydrogens=bool(spec.get("guess_hydrogens", True)),
        heavyatom_names=tuple(
            spec.get("heavyatom_names", ("N", "CA", "C", "O O1 OT1"))
        ),
        hydrogen_name=str(spec.get("hydrogen_name", "H HN HT1 HT2 HT3")),
    )
    kwargs = {
        "start": tr.get("start"),
        "stop": tr.get("stop"),
        "step": int(tr["stride"]),
    }
    backend = str(spec.get("backend", "serial"))
    if backend != "serial":
        kwargs["backend"] = backend
        if spec.get("n_workers") is not None:
            kwargs["n_workers"] = int(spec["n_workers"])
        if spec.get("n_parts") is not None:
            kwargs["n_parts"] = int(spec["n_parts"])
    analysis.run(verbose=bool(cfg.get("progress",{}).get("enabled",True)),**kwargs)
    residues = list(atoms.select_atoms("protein").residues)
    resids = np.asarray(analysis.results.resids)
    if len(residues) != len(resids):
        mapped, used = [], set()
        for resid in resids:
            candidates = [r for r in residues if r.resid == resid and r.ix not in used]
            if not candidates:
                raise RuntimeError(
                    "Could not map DSSP residues back to topology residues"
                )
            mapped.append(candidates[0])
            used.add(candidates[0].ix)
        residues = mapped
    codes = np.asarray(analysis.results.dssp).astype(str)
    codes[codes == "-"] = "C"
    times = np.asarray(analysis.times, float)
    if cfg.get("time", {}).get("dt_ps") is not None:
        # Keep DSSP on the same user-defined analyzed-frame clock as the
        # streamed analyses.
        times = float(cfg["time"].get(
            "origin_ps",
            0.0
        )) + np.arange(len(times)) * float(cfg["time"]["dt_ps"])
    nframes, nres = codes.shape
    out = ensure_dir(repout / "profiles" / "dssp")
    atomic_parquet(out / "assignments.parquet", pd.DataFrame({
        "time_ps": np.repeat(times, nres),
        "topology_residue_index": np.tile([int(r.ix) for r in residues], nframes),
        "residue": np.tile([_residue_label(r) for r in residues], nframes),
        "dssp": codes.ravel(),
    }))
