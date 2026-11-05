"""Concise manuscript-facing plots constructed ONLY from the standard statistics tables.

Do not independently fit uncertainty or pool correlated trajectory frames.
"""
from __future__ import annotations
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ..config import ensure_dir, presentation_enabled
from .figure_inventory import current_statistics

def plot_lipid_enrichment_matrix(cfg,root):
    if not current_statistics(cfg, root, 'time_series'):
        return
    source=Path(root)/"project_summary"/"statistics"/"time_series"/"replicate_statistics.parquet"
    if not source.exists() or not presentation_enabled(cfg,"plots","lipid_enrichment"):return
    frame=pd.read_parquet(source)
    if not {"system","metric","mean"}.issubset(frame):return
    selected=frame[frame["metric"].str.startswith("lipid_enrichment_")
                   &frame["metric"].str.endswith("_fold")].copy()
    if selected.empty:return
    selected["species"]=selected["metric"].str.replace("^lipid_enrichment_","",regex=True).str.replace("_(upper|lower)_fold$","",regex=True)
    selected["leaflet"]=selected["metric"].str.extract("_(upper|lower)_fold$")[0]
    # Replicate means are averaged equally; not a pooled per-frame average.
    selected["column"]=selected["system"].astype(str)+" | "+selected["leaflet"].astype(str)
    grouped=selected.groupby(["species","column"],dropna=False)["mean"].mean().reset_index()
    table=grouped.pivot(index="species",columns="column",values="mean")
    if table.empty:return
    values=table.to_numpy(float)
    signed=np.full_like(values,np.nan,dtype=float)
    good=np.isfinite(values)&(values>=0)
    # Plot 0 observed enrichment on a labeled finite clipping floor.
    signed[good]=np.log2(np.maximum(values[good],1e-3))
    fig,ax=plt.subplots(figsize=(max(7,len(table.columns)*1.5),max(3.5,len(table.index)*.58+2)))
    cmap=plt.get_cmap("RdBu_r").copy();cmap.set_bad("lightgrey")
    image=ax.imshow(signed,aspect="auto",vmin=-2.5,vmax=2.5,cmap=cmap)
    ax.set_xticks(range(len(table.columns)),table.columns,rotation=45,ha="right")
    ax.set_yticks(range(len(table.index)),table.index)
    ax.set_title("Lipid enrichment/depletion near protein (replicate means)")
    ax.set_xlabel("Condition and leaflet")
    ax.set_ylabel("Lipid species")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if np.isfinite(values[i,j]):
                ax.text(j,i,f"{values[i,j]:.2f}×",ha="center",va="center",fontsize=9)
    cb=fig.colorbar(image,ax=ax,pad=.025);cb.set_label("log₂ enrichment (display clipped)")
    fig.tight_layout()
    dest=ensure_dir(Path(root)/"project_summary"/"figures"/"membrane")
    fig.savefig(dest/"lipid_enrichment_heatmap.png",dpi=int(cfg.get("plots",{}).get("dpi",250)))
    plt.close(fig)



def plot_lipid_contacts(cfg, root, system_repouts):
    """Offline equal-replica mean neighbor matrices, distinct from H-bonds."""
    if (not cfg.get('plots', {}).get('enabled', True)
            or not cfg.get('plots', {}).get('families', {}).get('lipid_contacts', True)
            or not cfg.get('analyses', {}).get('lipid_contacts', {}).get('enabled', False)):
        return
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
    replicas.to_csv(dest / 'replicate_matrices.csv', index=False)

    def draw(values, title, metric, filename, *, difference=False):
        fig, ax = plt.subplots(figsize=(7, 6))
        cmap = plt.get_cmap('RdBu_r' if difference else 'viridis').with_extremes(bad='lightgrey')
        scale = {}
        if difference:
            finite = np.abs(values.to_numpy()[np.isfinite(values.to_numpy())])
            limit = float(finite.max()) if len(finite) and finite.max() > 0 else 1.
            scale = dict(vmin=-limit, vmax=limit)
        image = ax.imshow(values, cmap=cmap, **scale)
        ax.set_xticks(range(len(values.columns)), values.columns, rotation=60)
        ax.set_yticks(range(len(values)), values.index)
        ax.set_xlabel('Neighbor species B')
        ax.set_ylabel('Reference species A')
        ax.set_title(f'{title}: {metric.replace("_", " ")}')
        units = 'Neighbors per molecule' if metric == 'neighbors_per_molecule' else 'Contact share / possible-pair share'
        fig.colorbar(image, ax=ax, label=f'Difference in {units.lower()}' if difference else units)
        fig.tight_layout()
        fig.savefig(dest / filename, dpi=int(cfg.get('plots', {}).get('dpi', 180)))
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
