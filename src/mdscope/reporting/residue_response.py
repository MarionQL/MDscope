"""Cross-metric bootstrap contrast heatmaps from saved statistics."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from ..config import figure_enabled
from ..contracts.plot_validation import safe_plot_stem
from ..core.persistence import ensure_dir
from .figure_inventory import current_statistics


def plot_residue_response(cfg: dict, root: Path) -> None:
    """Render condition/chain panels; hollow dots mark pointwise CI exclusions.

    The highlighted regions were learned from the same data; they are not
    independently validated functional motifs or region-level significance.
    """
    if not figure_enabled(cfg, "profiles:residue_response:heatmap"):
        return
    if not current_statistics(cfg, root, "profiles"):
        return
    from .figure_inventory import residue_response_status
    if residue_response_status(cfg, root) not in {"valid", "legacy_unverified"}:
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
        matrix = values.to_numpy(float)
        if not np.isfinite(matrix).any():
            continue
        supported = group.pivot(index="profile", columns="resid",
                                values="ci_excludes_zero").reindex(
                                    index=values.index, columns=values.columns)
        from urllib.parse import quote

        from .plot_ids import figure_path, save_figure
        from .templates import heatmap_pages
        stem = safe_plot_stem(f"{left}__{right}__{chain}")
        for fig, ax, rs, cs, page in heatmap_pages(np.clip(matrix, -3, 3), list(values.index), list(values.columns),
                title=f"{chain}: {right} minus {left}\nScaled by within-residue replica SD; dots: saved pointwise CI excludes zero",
                units="Standardized condition difference (clipped ±3)", vmin=-3, vmax=3, cmap="RdBu_r"):
            supported_page = supported.fillna(False).to_numpy(bool)[rs, cs]
            points_y, points_x = np.where(supported_page & np.isfinite(matrix[rs, cs]))
            ax.scatter(points_x, points_y, marker=".", s=5, color="black", alpha=.75)
            ax.set_xlabel("Residue number (within chain)")
            ax.set_ylabel("Residue-resolved observable")
            for region_id in group.loc[group.region_id.notna(), "region_id"].unique():
                residues = group.loc[group.region_id == region_id, "resid"].unique()
                where = np.flatnonzero(np.isin(values.columns[cs], residues))
                if len(where):
                    ax.add_patch(Rectangle((where.min()-.5, -.5), len(where), len(values.index[rs]),
                                           fill=False, edgecolor="black", linewidth=1.2))
            fig.text(.5, .005, "Exploratory, data-selected regions; CI dots are pointwise, not region significance. Raw-unit profiles remain separate.", ha="center", fontsize=8)
            canonical = figure_path(cfg, dest / f"{stem}__heatmap.png", identifier="profiles:residue_response:heatmap")
            canonical = canonical.with_stem(canonical.stem + "@page_" + quote(stem, safe="-_.") + f"_{page:03d}")
            fig.tight_layout(rect=(0, .04, 1, 1))
            save_figure(fig, cfg, canonical, dpi=int(cfg.get("plots", {}).get("dpi", 250)), bbox_inches="tight")
            plt.close(fig)
