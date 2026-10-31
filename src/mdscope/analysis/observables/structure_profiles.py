"""Structural profile dispatcher and shared time-series/profile tables."""
from __future__ import annotations
from pathlib import Path
import MDAnalysis as mda
import numpy as np
import pandas as pd
from ..config import analysis_groups,atomic_parquet,ensure_dir
from .dssp import _run_dssp
from .sasa import _run_sasa
def run_structure_profiles(rep, cfg: dict, repout: Path):
    if not cfg["analyses"]["dssp"].get(
        "enabled",
        False
    ) and not cfg["analyses"]["sasa"].get("enabled", False):
        return
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    groups = analysis_groups(cfg)
    if groups:
        rows = []
        for group_id, group_spec in groups.items():
            atoms = universe.select_atoms(str(group_spec["selection"]))
            rows.extend(
                {
                    "group": group_id,
                    "topology_residue_index": int(residue.ix),
                }
                for residue in atoms.residues
            )
        if rows:
            atomic_parquet(
                ensure_dir(repout / "profiles") / "group_residues.parquet",
                pd.DataFrame(rows).drop_duplicates(),
            )
    if cfg["analyses"]["dssp"].get("enabled", False):
        _run_dssp(universe, cfg, repout)
    if cfg["analyses"]["sasa"].get("enabled", False):
        universe.trajectory[0]
        _run_sasa(universe, cfg, repout)

def postprocess_structure_profiles(
    repout: Path,
    cfg: dict
) -> tuple[list[pd.DataFrame], list[pd.DataFrame]]:
    """
    Return additional whole-system and group time-series produced by
    profile analyses.
    """
    whole, grouped = [], []
    raw_group = repout / "raw" / "group_residues.parquet"
    profile_group = repout / "profiles" / "group_residues.parquet"
    group_file = raw_group if raw_group.exists() else profile_group
    group_map = (
        pd.read_parquet(group_file)[
            ["group", "topology_residue_index"]
        ].drop_duplicates()
        if group_file.exists()
        else None
    )
    dssp_file = repout / "profiles" / "dssp" / "assignments.parquet"
    if dssp_file.exists():
        data = pd.read_parquet(dssp_file)
        if not data.empty:
            fractions = data.groupby(["time_ps", "dssp"]).size().unstack(fill_value=0)
            fractions = fractions.div(
                fractions.sum(axis=1),
                axis=0
            ).reset_index().rename(columns={
                "H": "helix_fraction",
                "E": "sheet_fraction",
                "C": "coil_fraction"
            })
            for col in ("helix_fraction", "sheet_fraction", "coil_fraction"):
                if col not in fractions:
                    fractions[col] = 0.0
            fractions = fractions[
                [
                    "time_ps",
                    "helix_fraction",
                    "sheet_fraction",
                    "coil_fraction",
                ]
            ]
            atomic_parquet(
                repout / "profiles" / "dssp" / "fraction_timeseries.parquet",
                fractions,
            )
            whole.append(fractions)
            occupancy = data.groupby([
                "topology_residue_index",
                "residue",
                "dssp"
            ]).size().rename("frames").reset_index()
            occupancy["occupancy"] = occupancy["frames"] / occupancy.groupby([
                "topology_residue_index",
                "residue"
            ])["frames"].transform("sum")
            atomic_parquet(
                repout / "profiles" / "dssp" / "residue_occupancy.parquet",
                occupancy,
            )
            if group_map is not None:
                gd = data.merge(group_map, on="topology_residue_index", how="inner")
                if not gd.empty:
                    gfrac = (
                        gd.groupby(["group", "time_ps", "dssp"])
                        .size()
                        .unstack(fill_value=0)
                    )
                    gfrac = gfrac.div(
                        gfrac.sum(axis=1),
                        axis=0
                    ).reset_index().rename(columns={
                        "H": "helix_fraction",
                        "E": "sheet_fraction",
                        "C": "coil_fraction"
                    })
                    for col in ("helix_fraction", "sheet_fraction", "coil_fraction"):
                        if col not in gfrac:
                            gfrac[col] = 0.0
                    gfrac = gfrac[[
                        "group",
                        "time_ps",
                        "helix_fraction",
                        "sheet_fraction",
                        "coil_fraction"
                    ]]
                    atomic_parquet(
                        repout
                        / "profiles"
                        / "dssp"
                        / "group_fraction_timeseries.parquet",
                        gfrac
                    )
                    grouped.append(gfrac)
    sasa_file = repout / "profiles" / "sasa" / "residue_sasa.parquet"
    if sasa_file.exists():
        data = pd.read_parquet(sasa_file)
        if not data.empty:
            total = data.groupby(
                "time_ps",
                as_index=False
            )["sasa_nm2"].sum().rename(columns={"sasa_nm2": "total_sasa_nm2"})
            atomic_parquet(
                repout / "profiles" / "sasa" / "total_timeseries.parquet",
                total,
            )
            whole.append(total)
            summary = data.groupby(
                ["topology_residue_index", "residue"]
            )["sasa_nm2"].agg(
                mean="mean",
                sd="std",
                median="median",
                min="min",
                max="max"
            ).reset_index()
            atomic_parquet(
                repout / "profiles" / "sasa" / "residue_summary.parquet",
                summary,
            )
            if group_map is not None:
                gd = data.merge(group_map, on="topology_residue_index", how="inner")
                if not gd.empty:
                    gtotal = gd.groupby(
                        ["group", "time_ps"],
                        as_index=False
                    )["sasa_nm2"].sum().rename(columns={"sasa_nm2": "total_sasa_nm2"})
                    atomic_parquet(
                        repout / "profiles" / "sasa" / "group_total_timeseries.parquet",
                        gtotal
                    )
                    grouped.append(gtotal)
    return whole, grouped


def temporal_structure_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    """
    Return per-frame residue matrices for SASA and DSSP occupancy
    profiles.
    """
    out = {}
    sasa_file = Path(repout) / "profiles" / "sasa" / "residue_sasa.parquet"
    if sasa_file.exists():
        data = pd.read_parquet(sasa_file)
        if not data.empty:
            residues = data[[
                "topology_residue_index",
                "residue"
            ]].drop_duplicates().sort_values("topology_residue_index")
            pivot = data.pivot(
                index="time_ps",
                columns="topology_residue_index",
                values="sasa_nm2"
            ).reindex(columns=residues.topology_residue_index)
            out["sasa_mean"] = {
                "times": pivot.index.to_numpy(float), "values": pivot.to_numpy(float),
                "positions": np.arange(
                    len(residues),
                    dtype=float
                ), "labels": residues.residue.to_numpy(str), "reducer": "mean",
            }
    dssp_file = Path(repout) / "profiles" / "dssp" / "assignments.parquet"
    if dssp_file.exists():
        data = pd.read_parquet(dssp_file)
        if not data.empty:
            residues = data[[
                "topology_residue_index",
                "residue"
            ]].drop_duplicates().sort_values("topology_residue_index")
            for code, label in (("H", "helix"), ("E", "sheet"), ("C", "coil")):
                indicator = data.assign(value=(data["dssp"] == code).astype(float))
                pivot = indicator.pivot(
                    index="time_ps",
                    columns="topology_residue_index",
                    values="value"
                ).reindex(columns=residues.topology_residue_index).fillna(0.0)
                out[f"dssp_{label}_occupancy"] = {
                    "times": pivot.index.to_numpy(float),
                    "values": pivot.to_numpy(float),
                    "positions": np.arange(
                        len(residues),
                        dtype=float
                    ), "labels": residues.residue.to_numpy(str), "reducer": "mean",
                }
    return out
