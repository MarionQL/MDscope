"""Diagnostic plots for SDC4 lipid maps and shared conformational states."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from ..config import presentation_enabled


def _plot_sdc4_replica_figures(cfg, root, system_repouts):
    if not cfg.get("plots", {}).get("enabled", True):
        return
    if not cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("enabled", False):
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for system, paths in system_repouts.items():
        for repout in paths:
            specificity_file = (Path(repout) / "extensions/sdc4_mechanisms"
                                / "lipid_region_specificity.parquet")
            if specificity_file.exists():
                all_contacts = pd.read_parquet(specificity_file)
                chemistry = all_contacts.loc[
                    all_contacts.moiety.isin(
                        ("glycan", "ceramide_head", "ceramide_tail")
                    )
                ]
                if not chemistry.empty:
                    chem = chemistry.groupby(
                        ["region", "species", "leaflet", "moiety"], as_index=False
                    ).agg(mean_fraction=("contact_fraction", "mean"))
                    groups = ("glycan", "ceramide_head", "ceramide_tail")
                    fig, axs = plt.subplots(
                        len(groups), 1,
                        figsize=(max(8, 1.2 * chem.region.nunique()), 9),
                        constrained_layout=True,
                    )
                    for ax, moiety in zip(axs, groups):
                        panel = chem.loc[chem.moiety.eq(moiety)]
                        for (species, leaflet), group in panel.groupby(
                            ["species", "leaflet"], sort=True
                        ):
                            ax.scatter(
                                group.region, group.mean_fraction,
                                s=25, label=f"{species} ({leaflet})"
                            )
                        ax.set_ylim(-.03, 1.03)
                        ax.set_ylabel(f"{moiety} contact fraction")
                        ax.tick_params(axis="x", rotation=60)
                        if not panel.empty:
                            ax.legend(fontsize="x-small", ncol=3)
                    fig.savefig(
                        specificity_file.with_name(
                            "lipid_region_moiety_contact_fraction.png"
                        ), dpi=180,
                    )
                    plt.close(fig)
                contacts = all_contacts.loc[all_contacts.moiety.eq("any")]
                if not contacts.empty:
                    # Old artifacts lack matched counts: leave the corrected
                    # metric missing; never substitute legacy enrichment.
                    contacts = contacts.copy()
                    for column in ("matched_local_contact_fraction",
                                   "matched_local_composition_enrichment"):
                        if column not in contacts:
                            contacts[column] = np.nan
                    summary = contacts.groupby(
                        ["region", "species", "leaflet"], as_index=False
                    ).agg(
                        mean_contact_fraction=("contact_fraction", "mean"),
                        mean_local_enrichment=("local_composition_enrichment", "mean"),
                        mean_matched_local_contact_fraction=("matched_local_contact_fraction", "mean"),
                        mean_matched_local_enrichment=("matched_local_composition_enrichment", "mean"),
                        n_matched_local_frames=("matched_local_contact_fraction", "count"),
                        n_matched_enrichment_frames=("matched_local_composition_enrichment", "count"),
                        n_frames_with_species=("contact_fraction", "count"),
                    )
                    summary.to_csv(
                        specificity_file.with_name("lipid_region_replicate_summary.csv"),
                        index=False,
                    )
                    for metric, ylabel in (
                        ("mean_contact_fraction", "Fraction of available molecules contacting"),
                        ("mean_local_enrichment", "Legacy unmatched local enrichment (comparison only)"),
                        ("mean_matched_local_contact_fraction", "Fraction contacting within the local anchor population"),
                        ("mean_matched_local_enrichment", "Matched local contact share / local composition share"),
                    ):
                        if summary[metric].notna().sum() == 0:
                            continue
                        fig, ax = plt.subplots(figsize=(max(7, len(summary.region.unique()) * 1.3), 4))
                        for i, (species, leaflet) in enumerate(
                            summary[["species", "leaflet"]].drop_duplicates().itertuples(index=False, name=None)
                        ):
                            group = summary.loc[
                                summary.species.eq(species) & summary.leaflet.eq(leaflet)
                            ]
                            ax.scatter(
                                group.region, group[metric],
                                marker=["o", "s", "^", "D", "v"][i % 5],
                                label=f"{species} ({leaflet})", s=25,
                            )
                        if metric.endswith("enrichment"):
                            ax.axhline(1., linestyle="--", color="black", alpha=.4)
                        ax.set_ylabel(ylabel)
                        ax.tick_params(axis="x", rotation=65)
                        ax.legend(fontsize="x-small", ncol=2)
                        fig.tight_layout()
                        fig.savefig(
                            specificity_file.with_name(f"lipid_region_{metric}.png"),
                            dpi=180,
                        )
                        plt.close(fig)
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

    # Exploratory within-replica NxIP coupling diagnostic. Avoid pooled
    # frame-level correlations: trajectories, not frames, are replicates.
    for system, paths in system_repouts.items():
        for repout in paths:
            folder = Path(repout) / "extensions/sdc4_mechanisms"
            series_path = folder / "timeseries.parquet"
            if not series_path.exists():
                continue
            table = pd.read_parquet(series_path)
            fig, axes = plt.subplots(1, 2, figsize=(9, 3.6),
                                     constrained_layout=True)
            n_panels = 0
            for ax, chain in zip(axes, ("A", "B")):
                xkey = f"sdc4_ectodomain_{chain}_height_A"
                occlusion = [
                    c for c in table.columns
                    if c.startswith(f"sdc4_motif_NxIP_PRO{chain}_")
                    and c.endswith("_steric_occlusion_fraction")
                ]
                ax.set_xlabel(f"Chain {chain} ectodomain height (Å)")
                ax.set_ylabel("NxIP steric occlusion fraction")
                if xkey not in table or len(occlusion) != 1:
                    ax.text(.5, .5, "NxIP mapping unavailable",
                            ha="center", va="center", transform=ax.transAxes)
                    continue
                x = table[xkey].to_numpy(float)
                y = table[occlusion[0]].to_numpy(float)
                valid = np.isfinite(x) & np.isfinite(y)
                if not valid.any():
                    continue
                ax.scatter(x[valid], y[valid], s=6, alpha=.18, rasterized=True)
                ax.set_ylim(-.03, 1.03)
                n_panels += 1
            if n_panels:
                fig.suptitle(f"{system}/{Path(repout).name}: exploratory NxIP–conformation coupling")
                fig.savefig(folder / "nxip_height_occlusion_joint.png", dpi=190)
            plt.close(fig)

            # Per-proline Ramachandran distributions, not averages of wrapped
            # angles. Keep individual replicas separate.
            proline_prefixes = sorted({
                name[:-len("_phi_deg")] for name in table.columns
                if name.startswith("sdc4_proline_") and name.endswith("_phi_deg")
                and name[:-len("_phi_deg")] + "_psi_deg" in table.columns
            })
            if proline_prefixes:
                fig, axs = plt.subplots(
                    1, len(proline_prefixes), squeeze=False,
                    figsize=(3.2 * len(proline_prefixes), 3.6),
                    constrained_layout=True,
                )
                for ax, label in zip(axs[0], proline_prefixes):
                    phi = table[label + "_phi_deg"].to_numpy(float)
                    psi = table[label + "_psi_deg"].to_numpy(float)
                    keep = np.isfinite(phi) & np.isfinite(psi)
                    ax.scatter(phi[keep], psi[keep], s=5, alpha=.15,
                               rasterized=True)
                    ax.set(xlim=(-180, 180), ylim=(-180, 180),
                           xlabel="φ (degrees)", ylabel="ψ (degrees)",
                           title=label.removeprefix("sdc4_proline_"))
                    ax.set_aspect("equal", adjustable="box")
                fig.savefig(folder / "vvhplvpl_proline_ramachandran.png", dpi=180)
                plt.close(fig)



def _plot_shared_state_figures(cfg, root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
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



def plot_sdc4_mechanisms(cfg, root, system_repouts):
    if not cfg.get('plots',{}).get('enabled',True) or not cfg.get('analyses',{}).get('sdc4_mechanisms',{}).get('enabled',False):
        return
    if presentation_enabled(cfg,'plots','sdc4'):
        _plot_sdc4_replica_figures(cfg,root,system_repouts)
    if presentation_enabled(cfg,'plots','shared_states'):
        _plot_shared_state_figures(cfg,root)
