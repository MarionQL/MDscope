"""Discover contiguous residue changes from complete replica profiles."""
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import atomic_parquet


_COLUMNS = [
    "profile", "system_1", "system_2", "chain", "start_residue", "end_residue",
    "start_label", "end_label", "n_residues", "n_replicates_1", "n_replicates_2",
    "mean_1", "mean_2", "difference_system2_minus_system1",
    "min_replica_separation", "rank",
]


def scan_residue_regions(profiles, system_replicates):
    """Rank runs of at least three residues with nonoverlapping replica ranges.

    Use complete finite replica coverage at each residue. Effects and separation
    are in the input profile's units. The screen is descriptive.
    """
    selected = profiles.profile.str.match(
        r"^(?:dssp_(?:helix|sheet|coil)_occupancy$|rmsf$|sasa_mean(?:__|$)"
        r"|sasa_buried_nm2__|residue_occlusion__|residue_membrane_distance__(?!.*__group_))", na=False
    )
    profiles = profiles.loc[selected]
    rows = []
    for profile, group in profiles.groupby("profile", sort=True):
        values = group.pivot(
            index="position_label", columns=["system", "replicate"], values="value",
        )
        residues = values.index.to_series().str.extract(
            r"^(?P<chain>[^:]+):[A-Za-z]+(?P<residue>-?\d+)$",
        ).dropna()
        residues = residues[residues.chain.ne("?")].copy()
        residues["residue"] = residues.residue.astype(int)
        residues = residues[~residues.duplicated(["chain", "residue"], keep=False)]
        residues = residues.sort_values(["chain", "residue"])
        for left, right in combinations(system_replicates, 2):
            left_keys = [(left, rep) for rep in system_replicates[left]]
            right_keys = [(right, rep) for rep in system_replicates[right]]
            if min(len(left_keys), len(right_keys)) < 2:
                continue
            a = values.reindex(index=residues.index, columns=left_keys).to_numpy(float)
            b = values.reindex(index=residues.index, columns=right_keys).to_numpy(float)
            complete = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
            higher = b.min(axis=1) - a.max(axis=1)
            lower = a.min(axis=1) - b.max(axis=1)
            candidates = residues.assign(
                direction=np.where(complete & (higher > 0), 1,
                                   np.where(complete & (lower > 0), -1, 0)),
                mean_1=a.mean(axis=1),
                mean_2=b.mean(axis=1),
                separation=np.maximum(higher, lower),
            )
            boundaries = (
                candidates.chain.ne(candidates.chain.shift())
                | candidates.residue.diff().ne(1)
                | candidates.direction.ne(candidates.direction.shift())
            )
            regions = []
            for _, region in candidates.groupby(boundaries.cumsum(), sort=False):
                if len(region) < 3 or region.direction.iloc[0] == 0:
                    continue
                mean_1, mean_2 = region[["mean_1", "mean_2"]].mean()
                regions.append({
                    "profile": profile, "system_1": left, "system_2": right,
                    "chain": region.chain.iloc[0],
                    "start_residue": region.residue.iloc[0],
                    "end_residue": region.residue.iloc[-1],
                    "start_label": region.index[0], "end_label": region.index[-1],
                    "n_residues": len(region),
                    "n_replicates_1": len(left_keys), "n_replicates_2": len(right_keys),
                    "mean_1": mean_1, "mean_2": mean_2,
                    "difference_system2_minus_system1": mean_2 - mean_1,
                    "min_replica_separation": region.separation.min(),
                })
            regions.sort(key=lambda row: -abs(row["difference_system2_minus_system1"]))
            rows.extend({**row, "rank": rank} for rank, row in enumerate(regions, 1))
    return pd.DataFrame(rows, columns=_COLUMNS)


def scan_saved_residue_regions(statistics_dir: Path) -> Path:
    """Scan an existing statistics export without accessing replica directories."""
    manifest = json.loads((statistics_dir / "manifest.json").read_text())
    if manifest.get("status") != "complete" or manifest.get("families", {}).get("profiles") != "enabled":
        raise ValueError("Residue region scanning requires completed profile statistics")
    replicas = {}
    for entry in manifest["replicate_paths"]:
        path = Path(entry)
        replicas.setdefault(path.parent.name, []).append(path.name)
    source = statistics_dir / "profiles" / "replicate_profiles.parquet"
    output = source.with_name("residue_regions.parquet")
    atomic_parquet(output, scan_residue_regions(pd.read_parquet(source), replicas))
    return output
