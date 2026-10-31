"""RMSF accumulator and profile reader. RMSF is deliberately separate from RMSD."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
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
        if self.n < 1:
            return np.full(self.mean.shape[0], np.nan)
        return np.sqrt((self.M2 / self.n).sum(axis=1))

def temporal_rmsf_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    """Load aligned coordinates needed for block-bootstrap RMSF profiles."""
    if not cfg["analyses"]["rmsf"].get("enabled", True):
        return {}
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
            "atom_residue_indices": (
                atoms["topology_residue_index"].to_numpy(int)
                if "topology_residue_index" in atoms
                else None
            ),
            "atom_residues": atoms["residue"].astype(str).to_numpy(),
            "reducer": "rmsf",
        }
    }
    group_file = Path(repout) / "raw" / "group_rmsf.parquet"
    if group_file.exists():
        group_atoms = pd.read_parquet(group_file)
        index_lookup = {int(atom): i for i, atom in enumerate(atoms["atom_index"].to_numpy(int))}
        for group_id, group in group_atoms.groupby("group", sort=False):
            columns = [
                index_lookup[
                    int(atom)
                ]for atom in group["atom_index"].to_numpy(int) if int(atom) in index_lookup
            ]
            if columns:
                out[f"rmsf__group_{group_id}"] = {
                    "times": times,
                    "coords": xyz[:, columns, :],
                    "atom_indices": atoms.iloc[columns]["atom_index"].to_numpy(int),
                    "atom_residue_indices": (
                        atoms.iloc[columns]["topology_residue_index"].to_numpy(int)
                        if "topology_residue_index" in atoms
                        else None
                    ),
                    "atom_residues": (
                        atoms.iloc[columns]["residue"].astype(str).to_numpy()
                    ),
                    "reducer": "rmsf",
                }
    return out
