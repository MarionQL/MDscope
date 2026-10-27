"""Hydrogen-bond event detection; shares normal interaction occupancy/lifetime processing."""
from __future__ import annotations
from pathlib import Path
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from .config import atomic_parquet,ensure_dir
from .interaction_tables import atom_label,chain_label,count_timeseries
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
