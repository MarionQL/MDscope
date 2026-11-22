"""RMSF accumulator and profile reader. RMSF is deliberately separate from RMSD."""
from __future__ import annotations

from pathlib import Path
from zipfile import BadZipFile, ZipFile

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

def temporal_rmsf_inputs(repout: Path):
    """Locate the persisted inputs shared by inference and readiness inspection."""
    cache = Path(repout) / "cache"
    return cache / "feature_atoms.parquet", sorted(cache.glob("aligned-part-*.npz"))


def temporal_rmsf_readiness(repout: Path, expected_times):
    """Inspect saved coordinate capability without opening a trajectory."""
    atom_file, parts = temporal_rmsf_inputs(repout)
    if not atom_file.is_file() or not parts:
        return False, "Temporal RMSF requires saved feature_atoms metadata and aligned coordinate cache"
    try:
        atoms = pd.read_parquet(atom_file)
        if (atoms.empty or not {"atom_index", "residue"} <= set(atoms) or
                not pd.api.types.is_integer_dtype(atoms["atom_index"]) or
                atoms["atom_index"].isna().any() or (atoms["atom_index"] < 0).any() or
                not atoms["atom_index"].is_unique):
            raise ValueError("empty or ambiguous feature atom identity")
        static = pd.read_parquet(Path(repout) / "raw/rmsf.parquet")
        identity = ["atom_index", "residue"]
        if "topology_residue_index" in atoms or "topology_residue_index" in static:
            identity.append("topology_residue_index")
            if ("topology_residue_index" not in atoms or
                    not pd.api.types.is_integer_dtype(atoms["topology_residue_index"]) or
                    atoms["topology_residue_index"].isna().any()):
                raise ValueError("invalid topology residue identity")
        if (not set(identity) <= set(atoms) or not set(identity) <= set(static) or
                not atoms[identity].equals(static[identity])):
            raise ValueError("feature atom identity disagrees with static RMSF")
        offset = 0
        for path in parts:
            with ZipFile(path) as archive, archive.open("xyz.npy") as payload:
                version = np.lib.format.read_magic(payload)
                if version == (1, 0):
                    shape, _, dtype = np.lib.format.read_array_header_1_0(payload)
                elif version == (2, 0):
                    shape, _, dtype = np.lib.format.read_array_header_2_0(payload)
                else:
                    raise ValueError(f"unsupported aligned cache array version {version}")
                if dtype.kind not in "fiu":
                    raise ValueError("aligned cache coordinates must be numeric")
                if archive.getinfo("xyz.npy").file_size != payload.tell() + int(np.prod(shape)) * dtype.itemsize:
                    raise ValueError("aligned cache coordinate payload is truncated")
                if archive.testzip() is not None:
                    raise ValueError("aligned cache coordinate payload failed CRC verification")
            with np.load(path, allow_pickle=False) as part:
                times = part["time_ps"]
                if times.ndim != 1 or shape != (len(times), len(atoms), 3):
                    raise ValueError("aligned cache frame/atom shape mismatch")
                if not np.array_equal(times, expected_times[offset:offset + len(times)]):
                    raise ValueError("aligned cache time grid disagrees with summary")
                offset += len(times)
        if offset != len(expected_times):
            raise ValueError("aligned cache frame coverage disagrees with summary")
    except (OSError, ValueError, KeyError, TypeError, EOFError, BadZipFile) as exc:
        return False, f"Invalid temporal RMSF inputs: {exc}"
    return True, "Validated temporal RMSF coordinate cache and feature atom metadata"


def temporal_rmsf_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    """Load aligned coordinates needed for block-bootstrap RMSF profiles."""
    if not cfg.get("analyses", {}).get("rmsf", {}).get("enabled", True):
        return {}
    atom_file, parts = temporal_rmsf_inputs(repout)
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
