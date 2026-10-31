"""Diagnostic plots for SDC4 lipid maps and shared conformational states."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


def plot_sdc4_mechanisms(cfg, root, system_repouts):
    if not cfg.get("plots", {}).get("enabled", True):
        return
    if not cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("enabled", False):
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for system, paths in system_repouts.items():
        for repout in paths:
            path = Path(repout) / "extensions/sdc4_mechanisms/lipid_density_maps.npz"
            if not path.exists():
                continue
            with np.load(path, allow_pickle=False) as data:
                species = data["species"].astype(str)
                leaf = data["leaflet"].astype(str)
                fold = data["enrichment_fold"]
                density = data["density_per_nm2"]
                x, y = data["x_edges_A"], data["y_edges_A"]
            for field, values, label in (
                ("density", density, "lipids / nm²"),
                ("enrichment", fold, "local / bulk lipid density"),
            ):
                n = len(species)
                fig, axes = plt.subplots(
                    n, 1, figsize=(6, max(3*n, 3)), squeeze=False,
                    constrained_layout=True,
                )
                for i, ax in enumerate(axes[:, 0]):
                    img = ax.pcolormesh(
                        x, y, values[i].T, shading="auto",
                        cmap="viridis", vmin=0,
                    )
                    ax.set_aspect("equal")
                    ax.set_xlabel("TM-centered x (Å), A → B")
                    ax.set_ylabel("TM-centered y (Å)")
                    ax.set_title(f"{species[i]} · {leaf[i]}")
                    fig.colorbar(img, ax=ax, label=label)
                fig.suptitle(f"{system}/{Path(repout).name}: protein-centered {field}")
                fig.savefig(path.parent / f"lipid_{field}_maps.png", dpi=180)
                plt.close(fig)

    state_dir = Path(root) / "project_summary/shared_states"
    pop_path = state_dir / "replicate_populations.parquet"
    if not pop_path.exists():
        return
    pop = pd.read_parquet(pop_path)
    if pop.empty:
        return
    systems = list(dict.fromkeys(pop["system"]))
    states = sorted(pop["state"].unique())
    fig, ax = plt.subplots(figsize=(max(7, len(states)*2.2), 4.5))
    width = .8 / max(len(systems), 1)
    for j, system in enumerate(systems):
        sub = pop[pop.system == system]
        means = [sub.loc[sub.state == state, "fraction"].mean() for state in states]
        xpos = np.arange(len(states)) + (j - (len(systems)-1)/2)*width
        ax.bar(xpos, means, width=width*.9, alpha=.6, label=system)
        for i, state in enumerate(states):
            values = sub.loc[sub.state == state, "fraction"].to_numpy(float)
            offsets = np.linspace(-width*.22, width*.22, len(values))
            ax.scatter(np.full(len(values), xpos[i]) + offsets, values,
                       color="black", s=13, zorder=3)
    ax.set_xticks(np.arange(len(states)), [f"State {s}" for s in states])
    ax.set_ylabel("Fraction of analyzed frames")
    ax.set_ylim(0, 1)
    ax.set_title("Shared-state populations (points = independent replicates)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(state_dir / "state_populations.png", dpi=200)
    plt.close(fig)
    diag = pd.read_parquet(state_dir / "sampling_diagnostics.parquet")
    fig, ax = plt.subplots(figsize=(max(6, len(diag)*1.1), 4))
    labels = [f"{row.system}/{row.replicate}" for row in diag.itertuples()]
    ax.bar(np.arange(len(diag)), diag.first_second_half_js_divergence)
    ax.set_xticks(np.arange(len(diag)), labels, rotation=45, ha="right")
    ax.set_ylabel("Jensen–Shannon divergence (bits)")
    ax.set_title("First- vs second-half state populations")
    fig.tight_layout()
    fig.savefig(state_dir / "state_sampling_diagnostics.png", dpi=200)
    plt.close(fig)
