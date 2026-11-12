"""Per-residue heavy-atom contact fractions, grouped by partner."""
from pathlib import Path

import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance, self_capped_distance

from ..interactions.interaction_tables import chain_label, residue_label

OCCLUSION_FILE = "residue_occlusion.parquet"
PARTNERS = ("protein", "hs", "lipid", "any")
FRACTION_COLUMNS = tuple(f"{partner}_fraction" for partner in PARTNERS)
IDENTITY_COLUMNS = ("time_ps", "topology_residue_index", "residue", "heavy_atom_count",
                    "contact_cutoff_A")


class ResidueOcclusion:
    """Count each contacted target atom once per partner and once in the union."""

    def __init__(self, protein, hs, lipids, cutoff):
        self.protein, self.hs, self.lipids = protein, hs, lipids
        self.cutoff = float(cutoff)
        if not len(protein) or not np.isfinite(self.cutoff) or self.cutoff <= 0:
            raise ValueError("Residue occlusion requires protein heavy atoms and a positive cutoff")
        indices = np.concatenate([protein.indices, hs.indices, lipids.indices])
        if len(np.unique(indices)) != len(indices):
            raise ValueError("Residue occlusion protein, HS and lipid selections must be disjoint")

        self.resindices, first, self.atom_to_residue, self.atom_counts = np.unique(
            protein.resindices, return_index=True, return_inverse=True, return_counts=True,
        )
        representatives = protein[first]
        self.labels = [residue_label(atom) for atom in representatives]
        if len(set(self.labels)) != len(self.labels):
            raise ValueError("Residue occlusion requires unique chain and residue labels")
        if any(chain_label(atom) == "?" for atom in representatives):
            raise ValueError("Residue occlusion requires chain IDs or segment IDs")
        self.chains = np.asarray([
            f"{atom.segindex}:{chain_label(atom)}" for atom in representatives
        ])[self.atom_to_residue]
        self.resids = representatives.resids[self.atom_to_residue]

    def frame(self, box):
        positions = self.protein.positions
        contacted = np.zeros((len(positions), len(PARTNERS)), dtype=bool)
        first, second = self_capped_distance(
            positions, self.cutoff, box=box, return_distances=False,
        ).T
        nonlocal_pairs = ((self.chains[first] != self.chains[second])
                          | (np.abs(self.resids[first] - self.resids[second]) > 2))
        contacted[first[nonlocal_pairs], 0] = True
        contacted[second[nonlocal_pairs], 0] = True
        for column, partner in enumerate((self.hs, self.lipids), 1):
            if len(partner):
                hits = capped_distance(positions, partner.positions, self.cutoff,
                                       box=box, return_distances=False)
                contacted[hits[:, 0], column] = True
        contacted[:, -1] = contacted[:, :-1].any(axis=1)
        counts = np.column_stack([
            np.bincount(self.atom_to_residue[mask], minlength=len(self.resindices))
            for mask in contacted.T
        ])
        return (counts / self.atom_counts[:, None]).astype(np.float32)

    def table(self, times, fractions):
        nframes = len(times)
        table = pd.DataFrame(
            np.concatenate(fractions, axis=0), columns=FRACTION_COLUMNS,
        )
        table.insert(0, "time_ps", np.repeat(np.asarray(times, dtype=float), len(self.resindices)))
        table.insert(1, "topology_residue_index", np.tile(self.resindices, nframes))
        table.insert(2, "residue", np.tile(self.labels, nframes))
        table.insert(3, "heavy_atom_count", np.tile(self.atom_counts, nframes))
        table.insert(4, "contact_cutoff_A", self.cutoff)
        return table


def temporal_residue_occlusion_profiles(repout: Path, cfg: dict) -> dict:
    """Read framewise partner fractions for the ordinary profile statistics."""
    spec = cfg.get("analyses", {}).get("sdc4_mechanisms", {})
    path = Path(repout) / "extensions/sdc4_mechanisms" / OCCLUSION_FILE
    if not spec.get("enabled") or not path.is_file():
        return {}
    data = pd.read_parquet(path, columns=[*IDENTITY_COLUMNS, *FRACTION_COLUMNS])
    if data.empty:
        raise ValueError(f"{path}: empty residue occlusion table")
    if not data["contact_cutoff_A"].eq(float(spec.get("contact_cutoff_A", 5.))).all():
        raise ValueError(f"{path}: residue occlusion cutoff differs from configuration")
    labels = data[["topology_residue_index", "residue"]].drop_duplicates().set_index(
        "topology_residue_index",
    )["residue"]
    pivot = data.pivot(
        index="time_ps", columns="topology_residue_index", values=list(FRACTION_COLUMNS),
    )
    return {
        f"residue_occlusion__{partner}": {
            "times": pivot.index.to_numpy(float),
            "values": pivot[column].to_numpy(float),
            "positions": np.arange(len(pivot[column].columns), dtype=float),
            "labels": labels.loc[pivot[column].columns].to_numpy(str),
            "reducer": "mean",
        }
        for partner, column in zip(PARTNERS, FRACTION_COLUMNS)
    }
