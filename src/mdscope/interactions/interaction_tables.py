"""Canonical interaction event identities, occupancy and lifetime tables."""
from __future__ import annotations
import glob
from pathlib import Path
import numpy as np
import pandas as pd
from ..config import atomic_parquet,ensure_dir
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
    "lipid_hydrogen_bonds": [
        "donor_molecule_id", "acceptor_molecule_id",
        "donor", "hydrogen", "acceptor", "interaction_class",
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
    lipid_hb_file=repout/"interactions"/"lipid_hydrogen_bonds"/"events.parquet"
    if cfg.get("analyses",{}).get("lipid_hydrogen_bonds",{}).get("enabled",False) and lipid_hb_file.exists():
        events=pd.read_parquet(lipid_hb_file)
        hb_dir=repout/"interactions"/"lipid_hydrogen_bonds"
        summarize_interaction_kind(events,
            ["donor_molecule_id","acceptor_molecule_id",
             "donor","hydrogen","acceptor","interaction_class"],
            base_times,gap_factor,hb_dir)
        counts_path=hb_dir/"count_timeseries.parquet"
        if counts_path.exists():
            counts.append(pd.read_parquet(counts_path))
    return counts


def interaction_tables(
    repout: Path,
    exclude_kinds: set[str] | None = None,
    table_names: tuple[str, ...] = ("presence", "occupancy", "lifetimes", "summary"),
) -> dict[str, dict[str, pd.DataFrame]]:
    out = {}
    excluded = exclude_kinds or set()
    base = Path(repout) / "interactions"
    if not base.exists():
        return out
    for kind_dir in sorted(path for path in base.iterdir() if path.is_dir()):
        if kind_dir.name in excluded:
            continue
        tables = {}
        for name in table_names:
            path = kind_dir / f"{name}.parquet"
            if path.exists():
                tables[name] = pd.read_parquet(path)
        if tables:
            out[kind_dir.name] = tables
    return out
