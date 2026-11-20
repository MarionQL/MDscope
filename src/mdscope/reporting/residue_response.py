"""Cross-metric bootstrap contrast heatmaps from saved statistics."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

from ..config import ensure_dir, figure_enabled
from .figure_inventory import current_statistics
from .plot_ids import safe_plot_stem


def plot_residue_response(cfg: dict, root: Path) -> None:
    """Render condition/chain panels; hollow dots mark pointwise CI exclusions.

    The highlighted regions were learned from the same data; they are not
    independently validated functional motifs or region-level significance.
    """
    if not figure_enabled(cfg, "profiles:residue_response:heatmap"):
        return
    if not current_statistics(cfg, root, "profiles"):
        return
    stats = root / "project_summary/statistics/profiles"
    source = stats / "residue_response_matrix.parquet"
    if not source.exists():
        return
    data = pd.read_parquet(source)
    if data.empty:
        return
    dest = ensure_dir(root / "project_summary/figures/profiles/residue_response")
    for (left, right, chain), group in data.groupby(
            ["system_1", "system_2", "chain"], sort=True):
        values = group.pivot(index="profile", columns="resid",
                             values="standardized_effect").sort_index()
        values = values.reindex(sorted(values.columns), axis=1)
        if values.empty or values.shape[1] < 2:
            continue
        supported = group.pivot(index="profile", columns="resid",
                                values="ci_excludes_zero").reindex(
                                    index=values.index, columns=values.columns)
        figure_width = max(10., min(26., 5. + .055 * len(values.columns)))
        figure_height = max(4., min(18., 1. + .30 * len(values.index)))
        fig, ax = plt.subplots(figsize=(figure_width, figure_height))
        matrix = values.to_numpy(float)
        image = ax.imshow(np.ma.masked_invalid(np.clip(matrix, -3, 3)),
                          interpolation="nearest", aspect="auto", cmap="RdBu_r",
                          norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3))
        ax.set_facecolor("#ececec")
        x = np.arange(len(values.columns))
        y = np.arange(len(values.index))
        # A small marker indicates that the existing pointwise bootstrap
        # interval excluded zero; it is not a post-clustering interval.
        points_y, points_x = np.where(supported.fillna(False).to_numpy(bool) &
                                      np.isfinite(matrix))
        ax.scatter(points_x, points_y, marker=".", s=5, color="black", alpha=.75)
        ticks = np.unique(np.r_[0, np.arange(0, len(x), max(1, len(x)//15)),
                                len(x)-1])
        ax.set_xticks(ticks, [str(values.columns[i]) for i in ticks], rotation=45)
        ax.set_yticks(y, values.index)
        ax.set_xlabel("Residue number (within chain)")
        ax.set_ylabel("Residue-resolved observable")
        ax.set_title(f"{chain}: {right} minus {left}\n"
                     "Scaled by within-residue replica SD; dots: pointwise 95% CI excludes zero")
        region_ids = group.loc[group.region_id.notna(), "region_id"].unique()
        for region_id in region_ids:
            res = group.loc[group.region_id == region_id, "resid"].unique()
            where = np.flatnonzero(np.isin(values.columns, res))
            if len(where):
                ax.add_patch(Rectangle((where.min() - .5, -.5), len(where),
                                       len(values.index), fill=False,
                                       edgecolor="black", linewidth=1.2))
        fig.colorbar(image, ax=ax, label="Standardized condition difference (clipped ±3)")
        fig.text(.5, .005, "Exploratory, data-selected regions; CI dots are pointwise, not region significance.",
                 ha="center", fontsize=8)
        fig.tight_layout(rect=(0, .025, 1, 1))
        stem = safe_plot_stem(f"{left}__{right}__{chain}")
        fig.savefig(dest / f"{stem}__heatmap.png",
                    dpi=int(cfg.get("plots", {}).get("dpi", 250)), bbox_inches="tight")
        plt.close(fig)
