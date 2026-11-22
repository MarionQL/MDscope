"""Bootstrapped lipid-species and lipid-class contact heatmaps."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import figure_enabled
from ..core.persistence import ensure_dir
from .figure_inventory import current_statistics
from .plot_ids import figure_dpi

# Short display names only. Statistical identities remain the original
# class or lipid-species keys (e.g., POPC, DOPC, GM1a).
_CLASS_LABELS = {
    "phosphatidylcholine": "PC",
    "phosphatidylethanolamine": "POPE",
    "anionic_phospholipids": "Anionic PL",
    "sphingomyelin": "PSM",
    "cholesterol": "CHL1",
    "glycosphingolipids": "GSL",
    "ceramide": "CER160",
}


def lipid_display_label(name: str, *, level: str) -> str:
    return _CLASS_LABELS.get(name, name) if level == "class" else name


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
    # Preserve room for long class names rather than enlarging data cells alone.
    level = "class" if stem.startswith("lipid_contact_class_") else "species"
    labels = [lipid_display_label(name, level=level) for name in ordered]
    longest = max(map(len, labels), default=0)
    size = max(7, min(17, 3 + .68 * len(ordered)))
    fig, ax = plt.subplots(figsize=(size + min(5, .085 * longest),
                                    size + min(2, .035 * longest)))
    cmap = plt.get_cmap("RdBu_r" if comparison else "viridis").copy()
    cmap.set_bad("#e0e0e0")
    limits = dict(vmin=-limit, vmax=limit) if comparison else dict(vmin=0, vmax=limit)
    image = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, **limits)
    ax.set_xticks(range(len(ordered)), labels, rotation=55, ha="right")
    ax.set_yticks(range(len(ordered)), labels)
    ax.tick_params(axis="both", labelsize=max(7, min(11, 13 - .18 * len(ordered))))
    ax.set_xlabel("Lipid " + ("class" if stem.startswith("lipid_contact_class_") else "species") + " B")
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
    fig.subplots_adjust(left=min(.38, max(.18, .085 + .009 * longest)),
                        bottom=min(.35, max(.17, .10 + .008 * longest)),
                        right=.88, top=.89)
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
