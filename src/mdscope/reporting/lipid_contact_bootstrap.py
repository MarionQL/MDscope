"""Bootstrapped lipid-species and lipid-class contact heatmaps."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import ensure_dir, figure_enabled
from .figure_inventory import current_statistics
from .plot_ids import figure_dpi


_UNITS = {
    "pair_count_undirected": "Unique contacting molecule pairs / frame",
    "neighbors_per_molecule": "Neighbor molecules / reference molecule",
    "enrichment": "Contact share / random-mixing pair share",
}


def _draw(cfg, dest: Path, table: pd.DataFrame, *, ordered: list[str],
          value: str, supported: bool, title: str, stem: str,
          metric: str, limit: float, comparison: bool):
    plot_id = f"lipid_contacts:{stem}"
    if not figure_enabled(cfg, plot_id):
        return
    values = table.pivot(index="species_a", columns="species_b", values=value)
    values = values.reindex(index=ordered, columns=ordered)
    matrix = values.to_numpy(float)
    if not np.isfinite(matrix).any():
        return
    size = max(7, min(15, 3 + .74 * len(ordered)))
    fig, ax = plt.subplots(figsize=(size, size - .2))
    cmap = plt.get_cmap("RdBu_r" if comparison else "viridis").copy()
    cmap.set_bad("#e0e0e0")
    limits = dict(vmin=-limit, vmax=limit) if comparison else dict(vmin=0, vmax=limit)
    image = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, **limits)
    ax.set_xticks(range(len(ordered)), ordered, rotation=55, ha="right")
    ax.set_yticks(range(len(ordered)), ordered)
    ax.set_xlabel("Lipid " + ("class" if "class_" in stem else "species") + " B")
    ax.set_ylabel("Reference lipid " + ("class" if "class_" in stem else "species") + " A")
    ax.set_title(title + "\n" + metric.replace("_", " "))
    if comparison and supported:
        flags = table.pivot(index="species_a", columns="species_b",
                            values="ci_excludes_zero").reindex(
                                index=ordered, columns=ordered,
                            ).fillna(False).to_numpy(bool)
        for i, j in zip(*np.where(flags & np.isfinite(matrix))):
            ax.text(j, i, "•", ha="center", va="center", color="black", fontsize=17)
    fig.colorbar(image, ax=ax, label=("System 2 − system 1: " if comparison else "")
                 + _UNITS[metric])
    fig.text(.5, .01, ("Dots: pointwise bootstrap CI excludes zero; "
                      "not corrected for multiple testing" if comparison
                      else "Equal-replica mean; grey cells have unavailable population or support"),
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .025, 1, 1))
    fig.savefig(dest / f"{stem}.png",
                dpi=figure_dpi(cfg, plot_id, int(cfg.get("plots", {}).get("dpi", 200))),
                bbox_inches="tight")
    plt.close(fig)


def plot_lipid_contact_bootstrap(cfg: dict, root: Path) -> bool:
    """Draw saved per-leaflet bootstrap matrices; no trajectory or inferential work.

    Returns False if no completed bootstrap output is available so the caller
    may fall back to descriptive legacy species plots.
    """
    if not current_statistics(cfg, root, "lipid_contacts"):
        return False
    folder = root / "project_summary/statistics/lipid_contacts"
    systems_path, differences_path = (folder / "system_statistics.parquet",
                                      folder / "system_comparisons.parquet")
    if not systems_path.is_file() or not differences_path.is_file():
        return False
    systems = pd.read_parquet(systems_path)
    comparisons = pd.read_parquet(differences_path)
    if systems.empty:
        return False
    dest = ensure_dir(root / "project_summary/figures/lipid_contacts")
    specs = cfg.get("analyses", {}).get("lipid_contacts", {}).get("species_selections", {})
    classes = cfg.get("analyses", {}).get("lipid_contacts", {}).get("classes", {})
    ordering = {"species": list(specs), "class": list(classes)}
    for level in ("species", "class"):
        prefix = "lipid_contact_" if level == "species" else "lipid_contact_class_"
        for metric in _UNITS:
            s = systems[(systems.level == level) & (systems.metric == metric)]
            c = comparisons[(comparisons.level == level) & (comparisons.metric == metric)]
            finite_s = s["mean"].to_numpy(float)
            finite_s = finite_s[np.isfinite(finite_s)]
            finite_c = c["difference_system2_minus_system1"].to_numpy(float)
            finite_c = np.abs(finite_c[np.isfinite(finite_c)])
            vmax = max(1e-12, float(finite_s.max())) if len(finite_s) else 1.
            diff_max = max(1e-12, float(finite_c.max())) if len(finite_c) else 1.
            for (system, leaflet), group in s.groupby(["system", "leaflet"], sort=True):
                species = set(group.species_a) | set(group.species_b)
                order = [v for v in ordering[level] if v in species]
                _draw(cfg, dest, group, ordered=order, value="mean", supported=False,
                      title=f"{system} · {leaflet} leaflet · {level}",
                      stem=f"{prefix}matrix_{system}_{leaflet}_{metric}",
                      metric=metric, limit=vmax, comparison=False)
            for (first, second, leaflet), group in c.groupby(
                    ["system_1", "system_2", "leaflet"], sort=True):
                species = set(group.species_a) | set(group.species_b)
                order = [v for v in ordering[level] if v in species]
                _draw(cfg, dest, group, ordered=order,
                      value="difference_system2_minus_system1",
                      supported=True,
                      title=f"{second} − {first} · {leaflet} leaflet · {level}",
                      stem=f"{prefix}difference_{first}__{second}_{leaflet}_{metric}",
                      metric=metric, limit=diff_max, comparison=True)
    return True
