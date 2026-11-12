"""Per-residue heavy-atom contact fractions, grouped by partner."""
from pathlib import Path

import numpy as np
import pandas as pd
from pyarrow.parquet import read_schema
from MDAnalysis.lib.distances import capped_distance, self_capped_distance

from ..interactions.interaction_tables import chain_label, residue_label

OCCLUSION_FILE = "residue_occlusion.parquet"
PARTNERS = ("protein", "hs", "lipid", "any")
WATER_COLUMNS = ("water_count", "polar_atom_count", "water_cutoff_A")
FRACTION_COLUMNS = tuple(f"{partner}_fraction" for partner in PARTNERS)
IDENTITY_COLUMNS = ("time_ps", "topology_residue_index", "residue", "heavy_atom_count",
                    "contact_cutoff_A")


def residue_environment_pending(path, spec):
    if not path.is_file():
        return True
    return spec.get("residue_water_cutoff_A") is not None and not set(WATER_COLUMNS).issubset(
        read_schema(path).names
    )


class ResidueOcclusion:
    """Count each contacted target atom once per partner and once in the union."""

    def __init__(self, protein, hs, lipids, cutoff, water=None, water_cutoff=None):
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
        self.water, self.water_cutoff = water, water_cutoff
        if water_cutoff is not None:
            if not len(water) or len(np.unique(water.resindices)) != len(water):
                raise ValueError("Residue hydration requires one oxygen per water molecule")
            polar = np.char.startswith(protein.names.astype(str), "N") | np.char.startswith(
                protein.names.astype(str), "O")
            self.polar_atoms = protein[polar]
            self.polar_residues = self.atom_to_residue[polar]
            self.polar_counts = np.bincount(self.polar_residues, minlength=len(self.resindices))
            if not len(self.polar_atoms):
                raise ValueError("Residue hydration requires protein N or O atoms")

    def water_frame(self, box):
        hits = capped_distance(self.polar_atoms.positions, self.water.positions, self.water_cutoff,
                               box=box, method="nsgrid", return_distances=False)
        n_water = len(self.water)
        pairs = np.unique(self.polar_residues[hits[:, 0]] * n_water + hits[:, 1])
        return np.bincount(pairs // n_water, minlength=len(self.resindices)).astype(np.int32)


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
        fractions = (counts / self.atom_counts[:, None]).astype(np.float32)
        if self.water_cutoff is not None:
            return np.column_stack((fractions, self.water_frame(box)))
        return fractions

    def table(self, times, fractions):
        nframes = len(times)
        values = np.concatenate(fractions, axis=0)
        table = pd.DataFrame(values[:, :len(FRACTION_COLUMNS)].astype(np.float32), columns=FRACTION_COLUMNS)
        if self.water_cutoff is not None:
            self.add_water(table, values[:, -1].astype(np.int32))
        table.insert(0, "time_ps", np.repeat(np.asarray(times, dtype=float), len(self.resindices)))
        table.insert(1, "topology_residue_index", np.tile(self.resindices, nframes))
        table.insert(2, "residue", np.tile(self.labels, nframes))
        table.insert(3, "heavy_atom_count", np.tile(self.atom_counts, nframes))
        table.insert(4, "contact_cutoff_A", self.cutoff)
        return table

    def add_water(self, table, counts):
        table["water_count"] = counts
        table["polar_atom_count"] = np.tile(self.polar_counts, len(table) // len(self.resindices))
        table["water_cutoff_A"] = self.water_cutoff


def temporal_residue_occlusion_profiles(repout: Path, cfg: dict) -> dict:
    """Read framewise partner fractions for the ordinary profile statistics."""
    spec = cfg.get("analyses", {}).get("sdc4_mechanisms", {})
    path = Path(repout) / "extensions/sdc4_mechanisms" / OCCLUSION_FILE
    if not spec.get("enabled") or not spec.get("residue_occlusion") or not path.is_file():
        return {}
    hydration = spec.get("residue_water_cutoff_A") is not None
    extra = list(WATER_COLUMNS) if hydration else []
    data = pd.read_parquet(path, columns=[*IDENTITY_COLUMNS, *FRACTION_COLUMNS, *extra])
    if data.empty:
        raise ValueError(f"{path}: empty residue occlusion table")
    if not data["contact_cutoff_A"].eq(float(spec.get("contact_cutoff_A", 5.))).all():
        raise ValueError(f"{path}: residue occlusion cutoff differs from configuration")
    if hydration and not data["water_cutoff_A"].eq(float(spec["residue_water_cutoff_A"])).all():
        raise ValueError(f"{path}: residue water cutoff differs from configuration")
    labels = data[["topology_residue_index", "residue"]].drop_duplicates().set_index(
        "topology_residue_index",
    )["residue"]
    pivot = data.pivot(
        index="time_ps", columns="topology_residue_index",
        values=[*FRACTION_COLUMNS, *extra[:2]],
    )
    columns = pivot[FRACTION_COLUMNS[0]].columns
    common = {"times": pivot.index.to_numpy(float), "positions": np.arange(len(columns), dtype=float),
              "labels": labels.loc[columns].to_numpy(str), "reducer": "mean"}
    out = {
        f"residue_occlusion__{partner}": {**common, "values": pivot[column].to_numpy(float)}
        for partner, column in zip(PARTNERS, FRACTION_COLUMNS)
    }
    if hydration:
        polar = pivot["polar_atom_count"].iloc[0].to_numpy() > 0
        water = pivot["water_count"].to_numpy(float)
        # No polar sites means the coordination observable is undefined.
        selected = {**common, "positions": common["positions"][polar], "labels": common["labels"][polar]}
        water = water[:, polar]
        out["residue_water_coordination__all"] = {**selected, "values": water}
        for partner, column in zip(PARTNERS, FRACTION_COLUMNS):
            occupancy = pivot[column].to_numpy() > 0
            if partner == "any":
                partner, occupancy = "no_partner", ~occupancy
            state = occupancy[:, polar]
            out[f"residue_water_contact_occupancy__{partner}"] = {**selected, "values": state}
            out[f"residue_water_coordination__during_{partner}"] = {
                **selected, "values": np.where(state, water, np.nan), "conditional": True,
            }
    return out
