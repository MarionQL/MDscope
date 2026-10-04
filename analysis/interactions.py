"""
Salt bridges, hydrogen bonds, ligand contacts, and shared interaction
event/lifetime handling.
"""

from __future__ import annotations

import glob
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from MDAnalysis.lib.distances import capped_distance, distance_array

from .config import atomic_parquet, ensure_dir


INTERACTION_KEY_COLUMNS = {
    "salt_bridges": ["residue_1", "residue_2", "interaction_class"],
    "ligand_residue": ["protein_residue"],
    "ligand_atom_residue": ["ligand_atom", "protein_residue"],
    "ligand_atom_atom": ["ligand_atom", "protein_atom"],
    "hydrogen_bonds": [
        "donor",
        "hydrogen",
        "acceptor",
        "interaction_class",
    ],
}


def interaction_key_columns(kind: str, table: pd.DataFrame | None = None) -> list[str]:
    """Return stable identity columns for a persisted interaction family."""
    if kind in INTERACTION_KEY_COLUMNS:
        return list(INTERACTION_KEY_COLUMNS[kind])
    if table is None or table.empty:
        return []
    excluded = {
        "time_ps",
        "distance_A",
        "min_distance_A",
        "angle_deg",
        "scope",
        "donor_chain",
        "acceptor_chain",
        "chain_1",
        "chain_2",
        "event_index",
        "start_ps",
        "end_ps",
        "lifetime_ps",
        "n_observations",
        "left_censored",
        "right_censored",
        "completed",
        "n_events",
        "n_event_starts",
        "n_completed_events",
        "censored_event_fraction",
        "event_rate_per_ns",
        "mean_lifetime_ps",
        "median_lifetime_ps",
        "max_lifetime_ps",
        "total_lifetime_ps",
        "frames_present",
        "occupancy",
    }
    return [column for column in table.columns if column not in excluded]


def chain_label(atom) -> str:
    return str(getattr(atom, "chainID", "") or getattr(atom, "segid", "") or "?")


def residue_label(atom) -> str:
    return f"{chain_label(atom)}:{atom.resname}{atom.resid}"


def atom_label(atom) -> str:
    return f"{residue_label(atom)}:{atom.name}"


def load_parts(pattern: Path) -> pd.DataFrame:
    files = sorted(glob.glob(str(pattern)))
    return pd.concat(
        [pd.read_parquet(path) for path in files],
        ignore_index=True
    ) if files else pd.DataFrame()


def sample_dt(times) -> float:
    times = np.sort(np.unique(np.asarray(times, float)))
    delta = np.diff(times)
    delta = delta[np.isfinite(delta) & (delta > 0)]
    return float(np.median(delta)) if len(delta) else 0.0


def event_lifetimes(
    events: pd.DataFrame,
    key_cols: list[str],
    all_times,
    gap_factor: float = 1.5
) -> pd.DataFrame:
    # Segment events from time gaps rather than Python frame-by-frame
    # state machines; this is both shorter and faster on large event
    # tables.
    cols = key_cols +[
        "event_index",
        "start_ps",
        "end_ps",
        "lifetime_ps",
        "n_observations",
        "left_censored",
        "right_censored",
        "completed"
    ]
    if events.empty:
        return pd.DataFrame(columns=cols)
    all_times = np.sort(np.unique(np.asarray(all_times, float)))
    dt = sample_dt(all_times)
    gap = gap_factor * dt if dt > 0 else np.inf
    first_t = float(all_times[0]) if len(all_times) else np.nan
    last_t = float(all_times[-1]) if len(all_times) else np.nan
    d = events[key_cols + ["time_ps"]].drop_duplicates().sort_values(key_cols + ["time_ps"]).copy()
    groupers = [d[c] for c in key_cols]
    d["event_index"] = d.groupby(
        key_cols,
        dropna=False
    )["time_ps"].diff().gt(gap).groupby(groupers, dropna=False).cumsum().fillna(0).astype(int) + 1
    life = d.groupby(
        key_cols + ["event_index"],
        dropna=False
    )["time_ps"].agg(start_ps="min", end_ps="max", n_observations="size").reset_index()
    life["lifetime_ps"] = life["end_ps"] - life["start_ps"] + dt
    life["left_censored"] = np.isclose(life["start_ps"], first_t)
    life["right_censored"] = np.isclose(life["end_ps"], last_t)
    life["completed"] = ~(life["left_censored"] | life["right_censored"])
    return life[cols]


def lifetime_summary(
    life: pd.DataFrame,
    key_cols: list[str],
    occupancy: pd.DataFrame,
    total_time_ps: float
) -> pd.DataFrame:
    cols = key_cols + [
        "n_events", "n_event_starts", "n_completed_events", "censored_event_fraction",
        "event_rate_per_ns",
        "mean_lifetime_ps", "median_lifetime_ps", "max_lifetime_ps", "total_lifetime_ps",
        "frames_present", "occupancy",
    ]
    if occupancy.empty:
        return pd.DataFrame(columns=cols)
    event_counts = life.groupby(
        key_cols,
        dropna=False
    ).size().rename("n_events") if not life.empty else pd.Series(dtype=float)
    observed_starts = life[~life["left_censored"]].groupby(
        key_cols,
        dropna=False
    ).size().rename("n_event_starts") if not life.empty else pd.Series(dtype=float)
    completed = life[life["completed"]] if not life.empty else life
    completed_stats = completed.groupby(key_cols, dropna=False)["lifetime_ps"].agg(
        n_completed_events="count", mean_lifetime_ps="mean", median_lifetime_ps="median",
        max_lifetime_ps="max"
    ) if not completed.empty else pd.DataFrame()
    out = occupancy.copy().set_index(key_cols)
    out["n_events"] = event_counts.reindex(out.index).fillna(0).astype(int)
    out["n_event_starts"] = observed_starts.reindex(out.index).fillna(0).astype(int)
    if not completed_stats.empty:
        out = out.join(completed_stats, how="left")
    else:
        out["n_completed_events"] = 0
        out["mean_lifetime_ps"] = np.nan
        out["median_lifetime_ps"] = np.nan
        out["max_lifetime_ps"] = np.nan
    out["n_completed_events"] = out["n_completed_events"].fillna(0).astype(int)
    out["censored_event_fraction"] = np.where(
        out["n_events"] > 0,
        1.0 - out["n_completed_events"] / out["n_events"],
        np.nan
    )
    out["event_rate_per_ns"] = (
        out["n_event_starts"] / (total_time_ps / 1000.0) if total_time_ps > 0 else np.nan
    )
    # Total observed interaction time is occupancy times the sampled
    # trajectory duration.
    out["total_lifetime_ps"] = out["occupancy"] * total_time_ps
    return out.reset_index()[cols]


def count_timeseries(
    events: pd.DataFrame,
    all_times,
    prefix: str,
    class_col: str | None = None,
    key_cols: list[str] | None = None
) -> pd.DataFrame:
    base = pd.DataFrame({"time_ps": np.asarray(all_times, float)})
    if events.empty:
        base[f"{prefix}_total_count"] = 0
        return base
    keys = key_cols or [col for col in events.columns if col not in {
        "time_ps",
        "distance_A",
        "min_distance_A",
        "angle_deg",
        "scope",
        "donor_chain",
        "acceptor_chain",
        "chain_1",
        "chain_2"
    }]
    present = events.drop_duplicates(["time_ps"] + keys)
    base = (
        base.merge(
            present.groupby("time_ps").size().rename(f"{prefix}_total_count"),
            left_on="time_ps",
            right_index=True,
            how="left",
        )
    )
    if class_col and class_col in present:
        for cls, group in present.groupby(class_col, sort=False):
            safe = str(cls).replace(" ", "_").replace("/", "_")
            base = (
                base.merge(
                    group.groupby("time_ps").size().rename(f"{prefix}_{safe}_count"),
                    left_on="time_ps",
                    right_index=True,
                    how="left",
                )
            )
    for col in base.columns[1:]:
        base[col] = pd.to_numeric(base[col], errors="coerce").fillna(0).astype(int)
    return base


def summarize_interaction_kind(
    events: pd.DataFrame,
    keys: list[str],
    all_times,
    gap_factor: float,
    outdir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    present = (
        events.drop_duplicates(
            ["time_ps"] + keys
        ) if not events.empty else pd.DataFrame(columns=["time_ps"] + keys)
    )
    times = np.sort(np.unique(np.asarray(all_times, float)))
    event_times = (
        np.sort(np.unique(events["time_ps"].to_numpy(float)))
        if not events.empty and "time_ps" in events
        else np.array([], dtype=float)
    )
    if not len(times):
        times = event_times
    elif len(event_times):
        indices = np.searchsorted(times, event_times)
        clipped = np.minimum(indices, len(times) - 1)
        matched = (indices < len(times)) & np.isclose(
            times[clipped],
            event_times,
        )
        if not np.all(matched):
            raise ValueError(
                "Interaction event times do not match the analysis time grid. "
                "Check time.dt_ps and trajectory slicing."
            )
    nframes = len(times)
    dt = sample_dt(times)
    total_time_ps = nframes * dt if dt > 0 else 0.0
    if present.empty or nframes == 0:
        occupancy = pd.DataFrame(columns=keys + ["frames_present", "occupancy"])
    else:
        occupancy = present.groupby(
            keys,
            dropna=False
        ).size().rename("frames_present").reset_index()
        occupancy["occupancy"] = occupancy["frames_present"] / nframes
    life = event_lifetimes(present, keys, times, gap_factor)
    summary = lifetime_summary(life, keys, occupancy, total_time_ps)
    ensure_dir(outdir)
    atomic_parquet(outdir / "presence.parquet", present)
    atomic_parquet(
        outdir / "occupancy.parquet",
        occupancy.sort_values("occupancy", ascending=False) if not occupancy.empty else occupancy
    )
    atomic_parquet(outdir / "lifetimes.parquet", life)
    atomic_parquet(outdir / "summary.parquet", summary)
    return occupancy, life, summary


def postprocess_interactions(
    repout: Path,
    cfg: dict,
    base_times: np.ndarray
) -> list[pd.DataFrame]:
    """
    Create the same presence/occupancy/lifetime/count products for every
    interaction family.
    """
    counts = []
    gap_factor = float(cfg.get("analyses", {}).get("lifetimes", {}).get("gap_factor", 1.5))

    if cfg["analyses"]["salt_bridges"].get("enabled", False):
        events = load_parts(repout / "raw" / "salt_bridges-part-*.parquet")
        summarize_interaction_kind(
            events,
            ["residue_1", "residue_2", "interaction_class"],
            base_times,
            gap_factor,
            repout / "interactions" / "salt_bridges"
        )
        counts.append(count_timeseries(
            events,
            base_times,
            "salt_bridge",
            "interaction_class",
            ["residue_1", "residue_2", "interaction_class"]
        ))

    if cfg["analyses"]["ligand_contacts"].get("enabled", False):
        events = load_parts(repout / "raw" / "ligand_contacts-part-*.parquet")
        for kind, keys in {
            "ligand_residue": ["protein_residue"],
            "ligand_atom_residue": ["ligand_atom", "protein_residue"],
            "ligand_atom_atom": ["ligand_atom", "protein_atom"],
        }.items():
            summarize_interaction_kind(
                events,
                keys,
                base_times,
                gap_factor,
                repout / "interactions" / kind
            )

    hb_dir = repout / "interactions" / "hydrogen_bonds"
    hb_file = hb_dir / "events.parquet"
    if cfg["analyses"]["hydrogen_bonds"].get("enabled", False) and hb_file.exists():
        hb = pd.read_parquet(hb_file)
        summarize_interaction_kind(
            hb,
            ["donor", "hydrogen", "acceptor", "interaction_class"],
            base_times,
            gap_factor,
            hb_dir
        )
        count_file = hb_dir / "count_timeseries.parquet"
        if count_file.exists():
            counts.append(pd.read_parquet(count_file))
    return counts


def interaction_tables(repout: Path) -> dict[str, dict[str, pd.DataFrame]]:
    out = {}
    base = Path(repout) / "interactions"
    if not base.exists():
        return out
    for kind_dir in sorted(path for path in base.iterdir() if path.is_dir()):
        tables = {}
        for name in ("presence", "occupancy", "lifetimes", "summary"):
            path = kind_dir / f"{name}.parquet"
            if path.exists():
                tables[name] = pd.read_parquet(path)
        if tables:
            out[kind_dir.name] = tables
    return out


# Interaction detectors
NEGATIVE_ATOMS = "(resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2)"
POSITIVE_ATOMS = "(resname LYS and name NZ) or (resname ARG and name NH1 NH2)"


def initialize_salt_bridges(universe, cfg: dict):
    spec = cfg["analyses"]["salt_bridges"]
    if not spec.get("enabled", False):
        return None
    return {
        "negative": universe.select_atoms(f"protein and ({NEGATIVE_ATOMS})"),
        "positive": universe.select_atoms(f"protein and ({POSITIVE_ATOMS})"),
        "cutoff_A": float(spec.get("cutoff_A", 4.0)),
    }


def salt_bridge_events(state, time_ps: float, box=None) -> list[tuple]:
    # capped_distance returns only pairs inside the cutoff, avoiding
    # construction of the full acidicxbasic distance matrix.
    if not state or not len(state["negative"]) or not len(state["positive"]):
        return []
    neg, pos = state["negative"], state["positive"]
    pairs, distances = capped_distance(
        neg.positions,
        pos.positions,
        max_cutoff=state["cutoff_A"],
        box=box,
        return_distances=True
    )
    best = {}
    for (i, j), distance in zip(pairs, distances):
        pair = tuple(sorted((residue_label(neg[int(i)]), residue_label(pos[int(j)]))))
        best[pair] = min(float(distance), best.get(pair, np.inf))
    return [
        (
            time_ps, r1, r2, r1.split(":", 1)[0], r2.split(":", 1)[0],
            "intrachain" if r1.split(":", 1)[0] == r2.split(":", 1)[0] else "interchain",
            (
                f"intra_{r1.split(':', 1)[0]}"
                if r1.split(":", 1)[0] == r2.split(":", 1)[0]
                else "interchain"
            ),
            distance,
        )
        for (r1, r2), distance in best.items()
    ]


def initialize_ligand_contacts(universe, protein, ligand, cfg: dict):
    spec = cfg["analyses"]["ligand_contacts"]
    if not spec.get("enabled", False) or ligand is None or not len(ligand):
        return None
    return {
        "protein": protein,
        "ligand": ligand,
        "cutoff_A": float(spec.get("cutoff_A", 4.0)),
        "save_atom_events": bool(spec.get("save_atom_events", True)),
        "ligand_labels": np.asarray(
            [atom_label(atom) for atom in ligand],
            dtype=object,
        ),
        "protein_atom_labels": np.asarray(
            [atom_label(atom) for atom in protein],
            dtype=object,
        ),
        "protein_residue_labels": np.asarray(
            [residue_label(atom) for atom in protein],
            dtype=object,
        ),
    }


def _minimum_distance(left, right, box=None, block_atoms: int = 4096) -> float:
    """Return the exact minimum distance without a full giant matrix."""
    if not len(left) or not len(right):
        return np.nan
    minimum = np.inf
    block_atoms = max(1, int(block_atoms))
    for start in range(0, len(right), block_atoms):
        stop = min(start + block_atoms, len(right))
        distances = distance_array(
            left.positions,
            right.positions[start:stop],
            box=box,
        )
        minimum = min(minimum, float(np.min(distances)))
    return minimum


def ligand_contact_events(state, time_ps: float, box=None):
    if not state:
        return np.nan, 0, []
    ligand, protein = state["ligand"], state["protein"]
    pairs, distances = capped_distance(
        ligand.positions,
        protein.positions,
        max_cutoff=state["cutoff_A"],
        box=box,
        return_distances=True,
    )
    if not len(pairs):
        return _minimum_distance(ligand, protein, box), 0, []

    li = pairs[:, 0].astype(int)
    pi = pairs[:, 1].astype(int)
    protein_residues = state["protein_residue_labels"]
    rows = []
    if state.get("save_atom_events", True):
        rows = list(
            zip(
                np.full(len(pairs), time_ps),
                state["ligand_labels"][li],
                state["protein_atom_labels"][pi],
                protein_residues[pi],
                distances.astype(float),
            )
        )
    return (
        float(np.min(distances)),
        int(len(np.unique(protein_residues[pi]))),
        rows,
    )


def run_hydrogen_bonds(rep, cfg: dict, repout: Path):
    spec = cfg["analyses"]["hydrogen_bonds"]
    if not spec.get("enabled", False):
        return
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    tr = cfg["trajectory"]
    between = spec.get("between")
    if between is None and spec.get("mode", "protein_protein") == "protein_ligand":
        ligand = cfg["selections"].get("ligand")
        if not ligand:
            raise ValueError("hydrogen_bonds.mode=protein_ligand requires selections.ligand")
        between = [cfg["selections"]["protein"], ligand]
    elif between is None:
        between = [cfg["selections"]["protein"], cfg["selections"]["protein"]]
    analysis = HydrogenBondAnalysis(
        universe,
        donors_sel=spec.get("donors_sel"),
        hydrogens_sel=spec.get("hydrogens_sel"),
        acceptors_sel=spec.get("acceptors_sel"),
        between=between,
        d_h_cutoff=float(spec.get("d_h_cutoff_A", 1.2)),
        d_a_cutoff=float(spec.get("d_a_cutoff_A", 3.0)),
        d_h_a_angle_cutoff=float(spec.get("angle_cutoff_deg", 150.0)),
        update_selections=bool(spec.get("update_selections", False)),
    )
    try:
        analysis.run(start=tr.get("start"), stop=tr.get("stop"), step=int(tr["stride"]))
    except Exception as exc:
        raise RuntimeError(
            "Hydrogen-bond analysis failed. Use a bond-aware topology or "
            "explicit donor/hydrogen/acceptor selections."
        ) from exc

    arr = np.asarray(analysis.results.hbonds)
    out = ensure_dir(repout / "interactions" / "hydrogen_bonds")
    analyzed_times = np.asarray(analysis.times, float)
    if cfg.get("time", {}).get("dt_ps") is not None:
        analyzed_times = (
            float(cfg["time"].get("origin_ps", 0.0))
            + np.arange(len(analyzed_times), dtype=float)
            * float(cfg["time"]["dt_ps"])
        )
    if arr.size == 0:
        atomic_parquet(
            out / "events.parquet",
            pd.DataFrame(
                columns=[
                    "time_ps",
                    "donor",
                    "hydrogen",
                    "acceptor",
                    "donor_chain",
                    "acceptor_chain",
                    "scope",
                    "interaction_class",
                    "distance_A",
                    "angle_deg",
                ]
            ),
        )
        atomic_parquet(
            out / "count_timeseries.parquet",
            pd.DataFrame({
                "time_ps": analyzed_times,
                "hydrogen_bond_total_count": 0,
            }),
        )
        return

    frame_to_time = dict(
        zip(np.asarray(analysis.frames, int), analyzed_times)
    )
    donor_idx, hydrogen_idx, acceptor_idx = arr[
        :,
        1
    ].astype(int), arr[:, 2].astype(int), arr[:, 3].astype(int)
    donors, hydrogens, acceptors = (
        universe.atoms[donor_idx], universe.atoms[hydrogen_idx], universe.atoms[acceptor_idx]
    )
    # Labels are materialized once and then NumPy-indexed for every
    # detected H-bond instead of rebuilding strings event by event.
    donor_labels = np.asarray([atom_label(a) for a in donors], object)
    hydrogen_labels = np.asarray([atom_label(a) for a in hydrogens], object)
    acceptor_labels = np.asarray([atom_label(a) for a in acceptors], object)
    donor_chains = np.asarray([chain_label(a) for a in donors], object)
    acceptor_chains = np.asarray([chain_label(a) for a in acceptors], object)
    same = donor_chains == acceptor_chains
    interaction_class = np.where(
        same,
        np.char.add("intra_", donor_chains.astype(str)),
        "interchain"
    )
    frame_numbers = arr[:, 0].astype(int)
    df = pd.DataFrame({
        "time_ps": np.fromiter(
            (frame_to_time.get(int(f), np.nan) for f in frame_numbers),
            float,
            count=len(frame_numbers)
        ),
        "donor": donor_labels,
        "hydrogen": hydrogen_labels,
        "acceptor": acceptor_labels,
        "donor_chain": donor_chains,
        "acceptor_chain": acceptor_chains,
        "scope": np.where(same, "intrachain", "interchain"),
        "interaction_class": interaction_class,
        "distance_A": arr[:, 4].astype(float),
        "angle_deg": arr[:, 5].astype(float),
    })
    atomic_parquet(out / "events.parquet", df)
    atomic_parquet(
        out / "count_timeseries.parquet",
        count_timeseries(
            df,
            analyzed_times,
            "hydrogen_bond",
            "interaction_class",
            ["donor", "hydrogen", "acceptor", "interaction_class"]
        )
    )
