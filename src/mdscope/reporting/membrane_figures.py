"""Concise manuscript-facing plots constructed ONLY from the standard statistics tables.

Do not independently fit uncertainty or pool correlated trajectory frames.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import ensure_dir, figure_enabled, plot_enabled
from .figure_inventory import current_statistics
from .plot_ids import figure_dpi


def plot_lipid_enrichment_matrix(cfg,root):
    if not any(figure_enabled(cfg, "lipid_enrichment:" + name)
               for name in ("lipid_enrichment_heatmap", "lipid_area_density_heatmap")):
        return
    if not current_statistics(cfg, root, 'time_series'):
        return
    source=Path(root)/"project_summary"/"statistics"/"time_series"/"replicate_statistics.parquet"
    if not source.exists() or not plot_enabled(cfg, "lipid_enrichment"):return
    frame=pd.read_parquet(source)
    if not {"system","metric","mean"}.issubset(frame):return
    for suffix, title, filename in (
        ("_fold", "Composition enrichment: local / leaflet mole fraction", "lipid_enrichment_heatmap.png"),
        ("_area_density_ratio", "Density ratio: protein-excluded planar area model", "lipid_area_density_heatmap.png"),
    ):
        if not figure_enabled(cfg, f"lipid_enrichment:{Path(filename).stem}"):
            continue
        selected=frame[frame["metric"].str.startswith("lipid_footprint_")
                       &frame["metric"].str.endswith(suffix)].copy()
        if selected.empty:continue
        selected[["species", "leaflet", "band"]] = selected["metric"].str.extract(
            rf"^lipid_footprint_(.+)_(upper|lower)_(.+A){suffix}$"
        )
        selected["species"] += " | " + selected["band"].str.replace("_", "–", regex=False)
        # Replicate means are averaged equally.
        selected["column"] = selected["system"].astype(str) + " | " + selected["leaflet"]
        grouped=selected.groupby(["species","column"],dropna=False)["mean"].mean().reset_index()
        table=grouped.pivot(index="species",columns="column",values="mean")
        if table.empty:continue
        values=table.to_numpy(float)
        if not np.isfinite(values).any():
            continue
        signed=np.full_like(values,np.nan,dtype=float)
        good=np.isfinite(values)&(values>=0)
        # Plot 0 observed enrichment on a labeled finite clipping floor.
        signed[good]=np.log2(np.maximum(values[good],1e-3))
        longest = max(map(len, table.index.astype(str)), default=0)
        fig,ax=plt.subplots(figsize=(max(9,len(table.columns)*1.7) + min(4, longest*.07),
                                     max(4,len(table.index)*.54+2)))
        cmap=plt.get_cmap("RdBu_r").copy();cmap.set_bad("lightgrey")
        image=ax.imshow(signed,aspect="auto",vmin=-2.5,vmax=2.5,cmap=cmap)
        ax.set_xticks(range(len(table.columns)),table.columns,rotation=45,ha="right")
        ax.set_yticks(range(len(table.index)),table.index)
        ax.set_title(title+" (equal replica means)")
        ax.set_xlabel("Condition and leaflet")
        ax.set_ylabel("Lipid species")
        if values.shape[0] * values.shape[1] <= 80:
            for i in range(values.shape[0]):
                for j in range(values.shape[1]):
                    if np.isfinite(values[i,j]):
                        ax.text(j,i,f"{values[i,j]:.2f}×",ha="center",va="center",fontsize=9)
        cb=fig.colorbar(image,ax=ax,pad=.025);cb.set_label("log₂ ratio (display clipped)")
        fig.subplots_adjust(left=min(.42, max(.19, .09 + .007 * longest)),
                            right=.88, bottom=.20, top=.90)
        dest=ensure_dir(Path(root)/"project_summary"/"figures"/"membrane")
        fig.savefig(dest/filename,dpi=figure_dpi(cfg, f"lipid_enrichment:{Path(filename).stem}",
                                                  int(cfg.get("plots",{}).get("dpi",250))))
        plt.close(fig)



def plot_lipid_contacts(cfg, root, system_repouts):
    """Offline equal-replica mean neighbor matrices, distinct from H-bonds."""
    if (not plot_enabled(cfg, "lipid_contacts")
            or not cfg.get('analyses', {}).get('lipid_contacts', {}).get('enabled', False)):
        return
    from .lipid_contact_bootstrap import plot_lipid_contact_bootstrap
    if plot_lipid_contact_bootstrap(cfg, Path(root)):
        return
    # Legacy saved-data products have no bootstrap CI: render descriptive
    # species panels only rather than implying inferential support.
    indices = ['leaflet', 'species_a', 'species_b']
    metrics = ['neighbors_per_molecule', 'enrichment']
    tables = []
    for system, paths in system_repouts.items():
        for repout in paths:
            path = Path(repout) / 'extensions/lipid_contacts/matrix.parquet'
            if path.exists():
                data = pd.read_parquet(path)
                columns = metrics + ['pair_count_undirected', 'molecules_a', 'molecules_b']
                tables.append(data.groupby(indices, as_index=False)[columns].mean()
                              .assign(system=system, replicate=Path(repout).name))
    if not tables:
        return
    replicas = pd.concat(tables, ignore_index=True)
    dest = ensure_dir(Path(root) / 'project_summary/figures/lipid_contacts')
    if cfg.get("plots", {}).get("include") is None:
        replicas.to_csv(dest / 'replicate_matrices.csv', index=False)

    def draw(values, title, metric, filename, *, difference=False):
        if not figure_enabled(cfg, f"lipid_contacts:{Path(filename).stem}"):
            return
        if values.empty or not np.isfinite(values.to_numpy(float)).any():
            return
        longest = max((len(str(v)) for v in (*values.index, *values.columns)), default=0)
        size = max(7, min(16, 3 + .7 * max(values.shape)))
        fig, ax = plt.subplots(figsize=(size + min(4, .08 * longest), size))
        cmap = plt.get_cmap('RdBu_r' if difference else 'viridis').with_extremes(bad='lightgrey')
        scale = {}
        if difference:
            finite = np.abs(values.to_numpy()[np.isfinite(values.to_numpy())])
            limit = float(finite.max()) if len(finite) and finite.max() > 0 else 1.
            scale = dict(vmin=-limit, vmax=limit)
        image = ax.imshow(values, cmap=cmap, **scale)
        ax.set_xticks(range(len(values.columns)), values.columns, rotation=55, ha="right")
        ax.set_yticks(range(len(values)), values.index)
        ax.tick_params(axis="both", labelsize=max(7, min(11, 13 - .18 * max(values.shape))))
        ax.set_xlabel('Neighbor species B')
        ax.set_ylabel('Reference species A')
        ax.set_title(f'{title}: {metric.replace("_", " ")}')
        units = 'Neighbors per molecule' if metric == 'neighbors_per_molecule' else 'Contact share / possible-pair share'
        fig.colorbar(image, ax=ax, label=f'Difference in {units.lower()}' if difference else units)
        fig.subplots_adjust(left=min(.4, max(.17, .09 + .008 * longest)),
                            right=.88, bottom=min(.38, max(.17, .10 + .008 * longest)),
                            top=.9)
        fig.savefig(dest / filename, dpi=figure_dpi(cfg, f"lipid_contacts:{Path(filename).stem}",
                                                    int(cfg.get('plots', {}).get('dpi', 180))))
        plt.close(fig)

    for (system, leaflet), group in replicas.groupby(['system', 'leaflet']):
        for metric in metrics:
            values = group.groupby(['species_a', 'species_b'])[metric].mean().unstack()
            draw(values, f'{system} / {leaflet}', metric,
                 f'lipid_contact_matrix_{system}_{leaflet}_{metric}.png')

    # Report a contrast only when every configured replica contributes a
    # finite estimate and contains both species in the corresponding leaflet.
    eligible = replicas.loc[(replicas.molecules_a > 0) & (replicas.molecules_b > 0)]
    grouped = eligible.groupby(['system'] + indices)[metrics]
    means, counts = grouped.mean(), grouped.count()
    condition = {
        system: means.loc[system].where(counts.loc[system] == len(paths))
        for system, paths in system_repouts.items() if system in means.index.get_level_values('system')
    }
    for first, second in combinations(condition, 2):
        delta = condition[first] - condition[second]
        for leaflet, group in delta.groupby(level='leaflet'):
            for metric in metrics:
                values = group.droplevel('leaflet')[metric].unstack()
                if values.notna().any().any():
                    draw(values, f'{first} minus {second} / {leaflet}', metric,
                         f'lipid_contact_difference_{first}__{second}_{leaflet}_{metric}.png',
                         difference=True)

def plot_project_diffusion(cfg, root):
    """Render independent replica estimates from the saved numerical summary."""
    if not figure_enabled(cfg, "diffusion:diffusion_coefficients"):
        return
    out = Path(root) / "project_summary" / "diffusion"
    source = out / "replicate_coefficients.parquet"
    if not source.is_file():
        return
    data = pd.read_parquet(source)
    data = data.loc[np.isfinite(
        pd.to_numeric(data["diffusion_nm2_ns"], errors="coerce"))]
    if data.empty:
        return
    groups = data.groupby(["system", "metric"], sort=True)
    fig, ax = plt.subplots(figsize=(max(8, 1.5 * groups.ngroups), 5))
    labels, ticks = [], []
    for i, ((system, metric), rows) in enumerate(groups):
        values = rows["diffusion_nm2_ns"].dropna().to_numpy(float)
        if len(values):
            ax.scatter(np.full(len(values), i), values, s=45, alpha=.85)
            ax.plot([i - .22, i + .22], [values.mean()] * 2, lw=2.5, color="black")
            labels.append(f"{system}\n{metric}")
            ticks.append(i)
    ax.set_xticks(ticks, labels, rotation=25, ha="right")
    ax.set_ylabel("Lateral diffusion coefficient (nm²/ns)")
    ax.set_title("Single-origin MSD replicate estimates")
    fig.tight_layout()
    fig.savefig(out / "diffusion_coefficients.png", dpi=figure_dpi(
        cfg, "diffusion:diffusion_coefficients", int(cfg.get("plots", {}).get("dpi", 250))))
    plt.close(fig)
