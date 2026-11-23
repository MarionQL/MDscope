"""Whole-protein interprotomer heavy-atom contact maps.

One row per observed residue pair/frame plus complete per-residue time series.
Only contacts across the two selected protein chains count: intra-chain packing,
HS and other nonprotein atoms are not silently included.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from ..core.persistence import atomic_parquet
from ..interactions.interaction_tables import chain_label, residue_label

PAIR_COLUMNS = [
    "time_ps", "topology_residue_index_A", "residue_A",
    "topology_residue_index_B", "residue_B", "atom_contact_pairs",
    "min_distance_A", "cutoff_A",
]
RESIDUE_COLUMNS = [
    "time_ps", "chain", "topology_residue_index", "residue",
    "protein_heavy_atom_count", "contacting_heavy_atom_fraction",
    "opposite_partner_residues", "atom_contact_pairs", "any_contact",
    "cutoff_A",
]


def _protein_chain(protein, tm):
    names = {chain_label(atom) for atom in tm}
    if len(names) != 1 or "?" in names:
        raise ValueError("TM selection must identify exactly one protein chain")
    name = next(iter(names))
    selection = protein[np.fromiter((chain_label(atom) == name for atom in protein),
                                    bool, count=len(protein))]
    if not len(selection):
        raise ValueError(f"No protein heavy atoms for TM chain {name}")
    return selection


class InterprotomerContacts:
    def __init__(self, protein, tm_chains, cutoff):
        if len(tm_chains) != 2:
            raise ValueError("Interprotomer contacts require exactly two TM chain selections")
        self.cutoff = float(cutoff)
        if not np.isfinite(self.cutoff) or self.cutoff <= 0:
            raise ValueError("Interprotomer contact cutoff must be positive")
        keys = tuple(tm_chains)
        self.a = _protein_chain(protein, tm_chains[keys[0]])
        self.b = _protein_chain(protein, tm_chains[keys[1]])
        if np.intersect1d(self.a.indices, self.b.indices).size:
            raise ValueError("Interprotomer protein groups overlap")
        self.keys = keys
        self.groups = {}
        for key, atoms in zip(keys, (self.a, self.b)):
            ids, first, inv, counts = np.unique(
                atoms.resindices, return_index=True, return_inverse=True,
                return_counts=True,
            )
            representatives = atoms[first]
            labels = [residue_label(atom) for atom in representatives]
            if len(set(labels)) != len(labels):
                raise ValueError("Ambiguous protein residue labels")
            self.groups[key] = {
                "indices": ids, "labels": np.asarray(labels, str),
                "inverse": inv, "counts": counts,
            }

    def frame(self, box, time_ps):
        """Return observed pairs and complete residue rows for one frame."""
        hits, distances = capped_distance(
            self.a.positions, self.b.positions, max_cutoff=self.cutoff,
            box=box, return_distances=True,
        )
        ga, gb = (self.groups[k] for k in self.keys)
        n_a, n_b = len(ga["indices"]), len(gb["indices"])
        count_a = np.zeros(n_a, np.int64)
        count_b = np.zeros(n_b, np.int64)
        touched_a = np.zeros(n_a, np.int64)
        touched_b = np.zeros(n_b, np.int64)
        partners_a = np.zeros(n_a, np.int64)
        partners_b = np.zeros(n_b, np.int64)
        pairs = pd.DataFrame(columns=PAIR_COLUMNS)
        if len(hits):
            ai = ga["inverse"][hits[:, 0]]
            bi = gb["inverse"][hits[:, 1]]
            count_a = np.bincount(ai, minlength=n_a)
            count_b = np.bincount(bi, minlength=n_b)
            touched_a = np.bincount(
                ga["inverse"][np.unique(hits[:, 0])], minlength=n_a)
            touched_b = np.bincount(
                gb["inverse"][np.unique(hits[:, 1])], minlength=n_b)
            grouped = pd.DataFrame({
                "ai": ai, "bi": bi, "distance": distances
            }).groupby(["ai", "bi"], sort=True)["distance"].agg(
                atom_contact_pairs="size", min_distance_A="min",
            ).reset_index()
            partners_a = np.bincount(grouped.ai.to_numpy(int), minlength=n_a)
            partners_b = np.bincount(grouped.bi.to_numpy(int), minlength=n_b)
            pairs = pd.DataFrame({
                "time_ps": time_ps,
                "topology_residue_index_A": ga["indices"][grouped.ai],
                "residue_A": ga["labels"][grouped.ai],
                "topology_residue_index_B": gb["indices"][grouped.bi],
                "residue_B": gb["labels"][grouped.bi],
                "atom_contact_pairs": grouped.atom_contact_pairs.to_numpy(int),
                "min_distance_A": grouped.min_distance_A.to_numpy(float),
                "cutoff_A": self.cutoff,
            }, columns=PAIR_COLUMNS)
        rows = []
        for key, g, counts, touched, partners in (
            (self.keys[0], ga, count_a, touched_a, partners_a),
            (self.keys[1], gb, count_b, touched_b, partners_b),
        ):
            rows.append(pd.DataFrame({
                "time_ps": time_ps, "chain": key,
                "topology_residue_index": g["indices"],
                "residue": g["labels"],
                "protein_heavy_atom_count": g["counts"],
                "contacting_heavy_atom_fraction": touched / g["counts"],
                "opposite_partner_residues": partners,
                "atom_contact_pairs": counts,
                "any_contact": (partners > 0).astype(np.int8),
                "cutoff_A": self.cutoff,
            }, columns=RESIDUE_COLUMNS))
        return pairs, pd.concat(rows, ignore_index=True)


def summarize_interprotomer_contacts(outdir: Path, pair_parts, residue_parts, times):
    """Validate frame completeness and aggregate bounded chunks before publishing."""
    from ..core.persistence import atomic_parquet_parts

    pair_parts, residue_parts = list(pair_parts), list(residue_parts)
    times = np.asarray(times, float)
    if (not pair_parts or len(pair_parts) != len(residue_parts)
            or not len(times) or not np.isfinite(times).all()
            or len(np.unique(times)) != len(times)):
        raise ValueError("Interprotomer summary requires matching chunks and unique finite times")

    expected_times = set(times)
    seen_times = set()
    identities = None
    residue_summaries, pair_summaries = [], []
    pair_keys = ["topology_residue_index_A", "residue_A",
                 "topology_residue_index_B", "residue_B"]
    residue_keys = ["chain", "topology_residue_index", "residue"]
    measures = {
        "contact_frames": ("any_contact", "sum"),
        "contact_atom_fraction": ("contacting_heavy_atom_fraction", "sum"),
        "partner_count": ("opposite_partner_residues", "sum"),
        "atom_pairs": ("atom_contact_pairs", "sum"),
    }
    for pair_part, residue_part in zip(pair_parts, residue_parts):
        residue = pd.read_parquet(residue_part)
        pairs = pd.read_parquet(pair_part)
        if not set(RESIDUE_COLUMNS).issubset(residue) or not set(PAIR_COLUMNS).issubset(pairs):
            raise ValueError("Interprotomer chunks lack required contact columns")
        chunk_times = set(residue.time_ps.unique())
        if not chunk_times or not chunk_times <= expected_times or chunk_times & seen_times:
            raise ValueError("Interprotomer chunk contains missing, unexpected or repeated times")
        seen_times.update(chunk_times)
        if (residue.duplicated(["time_ps", "topology_residue_index"]).any()
                or not pairs.time_ps.isin(chunk_times).all()
                or pairs.duplicated(["time_ps", *pair_keys[:1], pair_keys[2]]).any()):
            raise ValueError("Duplicate or out-of-grid interprotomer contact rows")
        signature = residue[[
            "chain", "topology_residue_index", "residue",
            "protein_heavy_atom_count", "cutoff_A",
        ]].drop_duplicates().sort_values("topology_residue_index").reset_index(drop=True)
        if signature.topology_residue_index.duplicated().any():
            raise ValueError("Changing interprotomer residue identity within a chunk")
        if identities is None:
            identities = signature
        elif not signature.equals(identities):
            raise ValueError("Inconsistent interprotomer residue identities between chunks")
        per_id = residue.groupby(residue_keys, sort=True, as_index=False).agg(
            frames=("time_ps", "size"), **measures,
        )
        if not per_id.frames.eq(len(chunk_times)).all():
            raise ValueError("Incomplete interprotomer residue frame grid")
        residue_summaries.append(per_id)
        if not pairs.empty:
            pair_summaries.append(pairs.groupby(pair_keys, sort=True, as_index=False).agg(
                frames_present=("time_ps", "size"),
                total_atom_pairs=("atom_contact_pairs", "sum"),
                min_distance_observed_A=("min_distance_A", "min"),
            ))
    if seen_times != expected_times:
        raise ValueError("Interprotomer frame grid differs from SDC4 timeseries")

    by_residue = pd.concat(residue_summaries).groupby(
        residue_keys, sort=True, as_index=False,
    ).agg({name: "sum" for name in ["frames", *measures]})
    if not by_residue.frames.eq(len(times)).all():
        raise ValueError("Incomplete interprotomer residue frame grid")
    for name, source in (
        ("contact_occupancy", "contact_frames"),
        ("mean_contacting_heavy_atom_fraction", "contact_atom_fraction"),
        ("mean_opposite_partner_residues", "partner_count"),
        ("mean_atom_contact_pairs", "atom_pairs"),
    ):
        by_residue[name] = by_residue[source] / by_residue.frames
    by_residue = by_residue[[*residue_keys, "frames", "contact_occupancy",
                             "mean_contacting_heavy_atom_fraction",
                             "mean_opposite_partner_residues", "mean_atom_contact_pairs"]]
    if pair_summaries:
        by_pair = pd.concat(pair_summaries).groupby(
            pair_keys, sort=True, as_index=False,
        ).agg(
            frames_present=("frames_present", "sum"),
            total_atom_pairs=("total_atom_pairs", "sum"),
            min_distance_observed_A=("min_distance_observed_A", "min"),
        )
        by_pair["contact_occupancy"] = by_pair.frames_present / len(times)
        by_pair["mean_atom_pairs_per_frame"] = by_pair.total_atom_pairs / len(times)
        by_pair = by_pair.drop(columns="total_atom_pairs")
    else:
        by_pair = pd.DataFrame(columns=[
            *pair_keys, "frames_present", "contact_occupancy",
            "mean_atom_pairs_per_frame", "min_distance_observed_A",
        ])

    outdir = Path(outdir)
    atomic_parquet_parts(outdir / "interprotomer_pair_frames.parquet", pair_parts)
    atomic_parquet_parts(outdir / "interprotomer_residue_frames.parquet", residue_parts)
    atomic_parquet(outdir / "interprotomer_residue_summary.parquet", by_residue)
    atomic_parquet(outdir / "interprotomer_pair_summary.parquet", by_pair)


def temporal_interprotomer_profiles(repout: Path, cfg: dict) -> dict:
    """Expose full-protomer residue response to the existing profile statistics."""
    spec = cfg.get("analyses", {}).get("sdc4_mechanisms", {})
    path = Path(repout) / "extensions/sdc4_mechanisms/interprotomer_residue_frames.parquet"
    if not spec.get("enabled") or not spec.get("interprotomer_contacts") or not path.is_file():
        return {}
    data = pd.read_parquet(path, columns=[
        "time_ps", "topology_residue_index", "residue", "any_contact",
        "opposite_partner_residues", "contacting_heavy_atom_fraction", "cutoff_A",
    ])
    if data.empty or not data.cutoff_A.eq(float(spec.get("contact_cutoff_A", 5.))).all():
        raise ValueError(f"{path}: empty or incompatible contact cutoff")
    if data.duplicated(["time_ps", "topology_residue_index"]).any():
        raise ValueError(f"{path}: duplicate residue frames")
    labels = data[["topology_residue_index", "residue"]].drop_duplicates()
    if labels.duplicated("topology_residue_index").any():
        raise ValueError(f"{path}: ambiguous residue labels")
    labels = labels.set_index("topology_residue_index")["residue"]
    outputs = {}
    for key, column in (
        ("interprotomer_contact_occupancy", "any_contact"),
        ("interprotomer_partner_count", "opposite_partner_residues"),
        ("interprotomer_atom_contact_fraction", "contacting_heavy_atom_fraction"),
    ):
        pivot = data.pivot(index="time_ps", columns="topology_residue_index", values=column)
        outputs[key] = {
            "times": pivot.index.to_numpy(float),
            "positions": np.arange(len(pivot.columns), dtype=float),
            "labels": labels.loc[pivot.columns].to_numpy(str),
            "values": pivot.to_numpy(float),
            "reducer": "mean",
        }
    return outputs
