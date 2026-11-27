"""
Shared plotting for time series, profiles, interactions, and
bootstrap/statistics tables.
"""

from __future__ import annotations

import colorsys
import math
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import (
    figure_enabled,
    plot_enabled,
    replicate_label,
    replicate_specs,
    system_label,
)
from ..contracts.plot_validation import safe_plot_stem as _safe
from ..core.persistence import archive_output, ensure_dir
from ..interactions import interaction_tables
from ..statistics.common import finite
from ..statistics.time_series import collect_time_series
from .figure_inventory import current_statistics
from .labels import _metric_meta
from .plot_ids import figure_dpi, figure_path, id_from_figure_path, save_figure
from .shared_scales import synchronize_y_axes
from .table_reports import add_table_page, table_page_enabled, table_report_session


def configure_plot_style(cfg: dict):
    fonts = cfg.get("plots", {}).get("font_sizes", {})
    plt.rcParams.update({
        "font.size": fonts.get("base", 13),
        "axes.titlesize": fonts.get("title", 16),
        "figure.titlesize": fonts.get("suptitle", 18),
        "legend.fontsize": fonts.get("legend", 10),
        "xtick.labelsize": fonts.get("ticks", 11),
        "ytick.labelsize": fonts.get("ticks", 11),
    })


def _save(fig, path: Path, cfg: dict):
    if path.parent.name not in {"primary", "secondary", "tertiary"}:
        path = figure_path(cfg, path)
    ensure_dir(path.parent)
    fig.tight_layout()
    dpi = figure_dpi(cfg, id_from_figure_path(path), cfg.get("plots", {}).get("dpi", 250))
    save_figure(fig, cfg, path, dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)


def _save_selected(fig, path: Path, cfg: dict, plot_id: str):
    if figure_enabled(cfg, plot_id):
        _save(fig, path, cfg)
    else:
        plt.close(fig)


def _system_colors(cfg: dict) -> dict[str, tuple]:
    cmap = plt.get_cmap("tab10")
    return {
        system: mcolors.to_rgb(cfg["systems"][system].get("color", cmap(i % 10))) for i,
        system in enumerate(cfg["systems"])
    }


def _replicate_colors(cfg: dict, system: str) -> dict[str, tuple]:
    reps = replicate_specs(cfg, system)
    base = _system_colors(cfg)[system]
    h, light, sat = colorsys.rgb_to_hls(*base)
    bounds = cfg.get("plots", {}).get("replicate_lightness_range", [0.28, 0.72])
    low, high = float(bounds[0]), float(bounds[1])
    levels = [light] if len(reps) == 1 else np.linspace(low, high, len(reps))
    return {
        rep["id"]: colorsys.hls_to_rgb(h, float(level), max(0.6, sat)) for rep,
        level in zip(reps, levels)
    }




def _profile_meta(profile_id: str) -> tuple[str, str, str]:
    if profile_id.startswith("lipid_order__"):
        parts = profile_id.split("__")
        title = "Lipid order: " + " ".join(parts[1:])
        return "Acyl-chain carbon position", "-S_CH", title
    return "Profile position", "Value", profile_id.replace("__", " * ").replace("_", " ")


def _subplot_grid(n: int):
    cols = min(3, max(1, n))
    rows = int(math.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4.2 * rows), squeeze=False)
    return fig, axes.ravel()


def _table_figure(title: str, dataframe: pd.DataFrame, path: Path, cfg: dict,
                  max_rows: int = 16, plot_id: str | None = None):
    identifier = plot_id or id_from_figure_path(path)
    if (cfg.get("plots", {}).get("include") is not None and identifier is None):
        return
    if identifier is not None and not table_page_enabled(cfg, identifier):
        return
    # Statistical tables are figures on purpose: they travel with the
    # plot set and are easy to drop directly into slides/manuscripts.
    if dataframe.empty:
        return
    if len(dataframe) > max_rows:
        for page, start in enumerate(range(0, len(dataframe), max_rows), 1):
            canonical = figure_path(cfg, path)
            paged = canonical.with_stem(canonical.stem + f"@page_{page:03d}")
            _table_figure(title + f" · page {page}", dataframe.iloc[start:start+max_rows],
                          paged, cfg, max_rows=max_rows, plot_id=identifier)
        return
    shown = dataframe.copy()
    from numbers import Real
    def display(value):
        if pd.isna(value):
            return ''
        return f'{value:.3g}' if isinstance(value, Real) and not isinstance(value, bool) else str(value)
    for col in shown.columns:
        shown[col] = shown[col].map(display)
    fig_height = max(2.8, 0.38 * (len(shown) + 2))
    fig, ax = plt.subplots(figsize=(max(8, 1.45 * len(shown.columns)), fig_height))
    ax.axis("off")
    ax.set_title(title, pad=12)
    plt.sca(ax)
    table = plt.table(
        cellText=shown.values,
        colLabels=shown.columns,
        cellLoc="center",
        loc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.25)
    table.auto_set_column_width(list(range(len(shown.columns))))
    if add_table_page(fig, path, identifier=identifier):
        plt.close(fig)
    else:
        _save(fig, path, cfg)



def rolling_time_mean(time_ps, values, window_ns: float):
    """Trailing time-weighted-by-observation mean; do not alter saved values.

    A sample contributes if its timestamp is within (t-window, t].
    Nonfinite measurements are excluded, not interpreted as zero.
    """
    times = np.asarray(time_ps, dtype=float)
    y = np.asarray(values, dtype=float)
    if times.ndim != 1 or y.shape != times.shape:
        raise ValueError("Rolling time series requires matching one-dimensional times and values")
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("Rolling time series requires finite, strictly increasing timestamps")
    if not len(times):
        return np.empty(0, dtype=float)
    lower = np.searchsorted(times, times - window_ns * 1000., side="right")
    valid = np.isfinite(y)
    counts = np.r_[0, np.cumsum(valid, dtype=np.int64)]
    totals = np.r_[0., np.cumsum(np.where(valid, y, 0.), dtype=float)]
    n = counts[1:] - counts[lower]
    result = np.full(len(times), np.nan)
    np.divide(totals[1:] - totals[lower], n, out=result, where=n > 0)
    return result

def _all_zero_lipid_footprint(metric, tables):
    """Suppress descriptive footprint plots only if every finite sample is zero.

    Missing values never count as measured zeros. Keep statistical tables and
    saved Parquet inputs even when descriptive plots are not informative.
    """
    if not str(metric).startswith("lipid_footprint_"):
        return False
    measured = [finite(table["value"]) for table in tables if table is not None]
    measured = [values for values in measured if len(values)]
    return bool(measured) and all(np.all(values == 0) for values in measured)


def _archive_suppressed_footprint_plots(root, folder, stem):
    """Preserve old zero-only PNGs outside the active figure inventory."""
    archive_root = root / ".archives" / "zero_only_footprint_figures"
    for variant in ("timeseries", "distribution", "autocorrelation"):
        path = folder / f"{stem}__{variant}.png"
        if path.is_file():
            archive_output(path, archive_root, path.name)


def plot_time_series(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]):
    if not plot_enabled(cfg, 'time_series'):
        return
    if not cfg.get("plots", {}).get("enabled", True):
        return
    stats_allowed = current_statistics(cfg, root, 'time_series')
    rolling_ns = cfg.get("plots", {}).get("time_series_rolling_window_ns")
    out = ensure_dir(root / "project_summary" / "figures" / "time_series")
    systems = list(system_repouts)
    keys = set()
    cache = {}
    include = cfg.get("plots", {}).get("include")
    sources = None if include is None else {"summary"}
    if include is not None:
        for identifier in include:
            if not identifier.startswith("time_series:"):
                continue
            stem = identifier.split(":", 2)[1]
            if "__group_lipid_contacts_" in stem:
                sources.add("lipid_contacts")
            if "__group_solvation_bridges_" in stem:
                sources.add("solvation_bridges")
            elif "__group_solvation_" in stem:
                sources.add("solvation")
            if "__group_hydration_" in stem:
                sources.add("hydration")
    for system, reps in system_repouts.items():
        for repout in reps:
            cache[(system, repout.name)] = collect_time_series(repout, sources=sources)
            keys.update(cache[(system, repout.name)].keys())
    rep_stats_path = (
        root / "project_summary" / "statistics" / "time_series" / "replicate_statistics.parquet"
    )
    sys_stats_path = (
        root / "project_summary" / "statistics" / "time_series" / "system_statistics.parquet"
    )
    needs_replicate_stats = (include is None or any(
        item.startswith("time_series:") and item.rsplit(":", 1)[-1] in
        {"autocorrelation", "bootstrap_table"} for item in include))
    needs_system_stats = (include is None or any(
        item.startswith("time_series:") and item.endswith(":bootstrap_table")
        for item in include))
    needs_comparisons = (include is None or any(
        item.startswith("time_series:") and item.endswith(":comparison_table")
        for item in include))
    rep_stats = (pd.read_parquet(rep_stats_path) if rep_stats_path.exists()
                 and stats_allowed and needs_replicate_stats else pd.DataFrame())
    sys_stats = (pd.read_parquet(sys_stats_path) if sys_stats_path.exists()
                 and stats_allowed and needs_system_stats else pd.DataFrame())
    comp_stats_path = (
        root / "project_summary" / "statistics" / "time_series" / "system_comparisons.parquet"
    )
    comp_stats = (pd.read_parquet(comp_stats_path) if comp_stats_path.exists()
                  and stats_allowed and needs_comparisons else pd.DataFrame())

    stems = {}
    for group, metric in keys:
        stem = _safe(metric + (f"__group_{group}" if group is not None else ""))
        if stem in stems and stems[stem] != (group, metric):
            raise ValueError(f"Plot ID collision for {stem}: {stems[stem]!r} and {(group, metric)!r}; rename the groups")
        stems[stem] = (group, metric)
    for group_id, metric in sorted(keys, key=lambda x: (str(x[0]), x[1])):
        title, ylabel = _metric_meta(metric)
        suffix = f"__group_{group_id}" if group_id is not None else ""
        stem = _safe(metric + suffix)
        selected = cfg.get("plots", {}).get("include")
        if selected is not None and not any(
                identifier.startswith(f"time_series:{stem}:") for identifier in selected):
            continue
        # A registered metric may have no finite observations in any replica.
        if not any(
            len(finite(table["value"]))
            for system in systems for repout in system_repouts[system]
            if (table := cache[(system, repout.name)].get((group_id, metric))) is not None
        ):
            continue
        suppress_zero = _all_zero_lipid_footprint(
            metric, (
                cache[(system, repout.name)].get((group_id, metric))
                for system in systems for repout in system_repouts[system]
            ),
        )
        if suppress_zero:
            _archive_suppressed_footprint_plots(root, out, stem)
        else:
            fig, axes = _subplot_grid(len(systems))
            distribution_values, distribution_labels, distribution_colors = [], [], []
            for ax, system in zip(axes, systems):
                colors = _replicate_colors(cfg, system)
                for repout in system_repouts[system]:
                    table = cache[(system, repout.name)].get((group_id, metric))
                    if table is None:
                        continue
                    values = finite(table["value"])
                    if not len(values):
                        continue
                    plotted = (rolling_time_mean(table["time_ps"], table["value"], rolling_ns)
                               if rolling_ns is not None else table["value"])
                    if np.isfinite(np.asarray(plotted, dtype=float)).any():
                        ax.plot(
                            table["time_ps"] / 1000.0,
                            plotted,
                            lw=1.1,
                            color=colors.get(repout.name),
                            label=replicate_label(cfg, system, repout.name)
                        )
                    distribution_values.append(values)
                    distribution_labels.append(
                        f"{system_label(cfg, system)}\n"
                        f"{replicate_label(cfg, system, repout.name)}"
                    )
                    distribution_colors.append(colors.get(repout.name))
                ax.set_title(system_label(cfg, system))
                ax.set_xlabel("Time (ns)")
                ax.set_ylabel(ylabel)
                if ax.has_data():
                    ax.legend(frameon=False)
                else:
                    ax.text(.5, .5, "No finite data", ha="center", va="center",
                            transform=ax.transAxes, color="0.45")
            synchronize_y_axes(axes[:len(systems)])
            for ax in axes[len(systems):]:
                ax.remove()
            fig.suptitle(title + (f" - {group_id}" if group_id is not None else "")
                         + (f" · trailing {rolling_ns:g} ns mean" if rolling_ns is not None else ""))
            _save_selected(fig, out / f"{stem}__timeseries.png", cfg, f"time_series:{stem}:timeseries")

            if cfg.get("plots", {}).get("distributions", True) and distribution_values:
                fig, (
                    ax_hist,
                    ax_violin,
                ) = plt.subplots(
                    2,
                    1,
                    figsize=(max(8, 1.1 * len(distribution_values)), 7),
                    gridspec_kw={"height_ratios": [2.2, 1]}
                )
                pooled = np.concatenate(distribution_values)
                bins = np.histogram_bin_edges(pooled, bins="fd") if len(pooled) > 2 else 10
                for values, label, color in zip(
                    distribution_values,
                    distribution_labels,
                    distribution_colors
                ):
                    ax_hist.hist(
                        values,
                        bins=bins,
                        density=True,
                        histtype="step",
                        linewidth=1.3,
                        color=color,
                        label=label.replace("\n", " / ")
                    )
                ax_hist.set_ylabel("Density")
                ax_hist.legend(frameon=False, ncol=2)
                parts = ax_violin.violinplot(distribution_values, showmeans=False, showmedians=True)
                for body, color in zip(parts["bodies"], distribution_colors):
                    body.set_facecolor(color)
                    body.set_alpha(0.6)
                ax_violin.set_xticks(
                    np.arange(1, len(distribution_labels) + 1),
                    distribution_labels,
                    rotation=25,
                    ha="right"
                )
                ax_violin.set_ylabel(ylabel)
                fig.suptitle(f"{title} distribution")
                _save_selected(fig, out / f"{stem}__distribution.png", cfg, f"time_series:{stem}:distribution")

            if cfg.get("plots", {}).get("autocorrelation", True) and not rep_stats.empty:
                acf_path = (
                    root / "project_summary" / "statistics" / "time_series" / "autocorrelation.parquet"
                )
                if acf_path.exists() and stats_allowed:
                    acf = pd.read_parquet(acf_path)
                    mask = (acf.metric == metric) & (acf.group.fillna("") == (group_id or ""))
                    subset = acf[mask]
                    if not subset.empty and np.isfinite(
                            subset["autocorrelation"].to_numpy(float)).any():
                        fig, axes = _subplot_grid(len(systems))
                        for ax, system in zip(axes, systems):
                            colors = _replicate_colors(cfg, system)
                            for rep_id, data in subset[subset.system == system].groupby(
                                "replicate",
                                sort=False
                            ):
                                ax.plot(
                                    data.lag_ps / 1000.0,
                                    data.autocorrelation,
                                    color=colors.get(str(rep_id)),
                                    lw=1.2,
                                    label=replicate_label(cfg, system, str(rep_id))
                                )
                            ax.axhline(0, lw=0.8, color="0.5")
                            ax.set_title(system_label(cfg, system))
                            ax.set_xlabel("Lag (ns)")
                            ax.set_ylabel("Autocorrelation")
                            if ax.lines and any(
                                    np.isfinite(line.get_ydata()).any() for line in ax.lines):
                                ax.legend(frameon=False)
                        synchronize_y_axes(axes[:len(systems)], symmetric=True)
                        for ax in axes[len(systems):]:
                            ax.remove()
                        fig.suptitle(f"{title} autocorrelation")
                        _save_selected(fig, out / f"{stem}__autocorrelation.png", cfg, f"time_series:{stem}:autocorrelation")

        if not cfg.get("plots", {}).get("table_reports_pdf", False):
            if cfg.get("plots", {}).get("bootstrap_tables", True) and not rep_stats.empty:
                r = (
                    rep_stats[
                        (rep_stats.metric == metric) & (rep_stats.group.fillna("") == (group_id or ""))
                    ].copy()
                )
                s = (
                    sys_stats[
                        (sys_stats.metric == metric) & (sys_stats.group.fillna("") == (group_id or ""))
                    ].copy() if not sys_stats.empty else pd.DataFrame()
                )
                rows = []
                for _, row in r.iterrows():
                    rows.append({
                        "level": f"{row.system}/{row.replicate}", "mean": row["mean"],
                        "mean CI": f"{row.mean_ci_low:.3g}-{row.mean_ci_high:.3g}",
                        "median": row["median"],
                        "tau(fr)": row.tau_int_frames, "tau(ps)": row.tau_int_ps,
                        "N_eff": row.effective_n,
                        "block(fr)": int(row.block_length_frames),
                        "blocks/pseudo": int(row.blocks_per_pseudostudy),
                    })
                for _, row in s.iterrows():
                    rows.append({
                        "level": f"{row.system} (system)", "mean": row["mean"],
                        "mean CI": f"{row.mean_ci_low:.3g}-{row.mean_ci_high:.3g}",
                        "median": row["median"],
                        "tau(fr)": (
                            f"{row.tau_int_frames_median:.3g} "
                            f"[{row.tau_int_frames_min:.3g},"
                            f"{row.tau_int_frames_max:.3g}]"
                        ),
                        "tau(ps)": "replicate-specific", "N_eff": np.nan,
                        "block(fr)": "replicate-specific",
                        "blocks/pseudo": (
                            f"{row.blocks_per_pseudostudy_median:.3g} "
                            f"[{int(row.blocks_per_pseudostudy_min)},"
                            f"{int(row.blocks_per_pseudostudy_max)}]"
                        ),
                    })
                _table_figure(
                    f"{title}: hierarchical bootstrap statistics",
                    pd.DataFrame(rows),
                    out / f"{stem}__bootstrap_table.png",
                    cfg,
                )

            if not comp_stats.empty:
                c = (
                    comp_stats[
    (
                            comp_stats.metric == metric
                        ) & (comp_stats.group.fillna("") == (group_id or ""))
                    ].copy()
                )
                if not c.empty:
                    shown = pd.DataFrame({
                        "comparison": c.apply(
                            lambda r: (
                                f"{system_label(cfg, r.system_2)} - "
                                f"{system_label(cfg, r.system_1)}"
                            ),
                            axis=1,
                        ),
                        "mean delta": c.mean_difference_system2_minus_system1,
                        "95% CI": c.apply(
                            lambda r: (
                                f"{r.mean_difference_ci_low:.3g}-"
                                f"{r.mean_difference_ci_high:.3g}"
                            ),
                            axis=1,
                        ),
                        "P(delta>0)": c.mean_probability_difference_gt_0,
                        "p": c.mean_p_value,
                        "q": c["mean_q_value"] if "mean_q_value" in c else np.nan,
                    })
                    _table_figure(
                        f"{title}: system comparison statistics",
                        shown,
                        out / f"{stem}__comparison_table.png",
                        cfg,
                    )


def plot_profiles(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]):
    if not plot_enabled(cfg, 'profiles'):
        return
    selected = cfg.get("plots", {}).get("include")
    if selected is not None and not any(
            identifier.startswith("profiles:") and not identifier.endswith(":distance_heatmap")
            for identifier in selected):
        return
    stats_allowed = current_statistics(cfg, root, 'profiles')
    if not stats_allowed:
        return
    out = ensure_dir(root / "project_summary" / "figures" / "profiles")
    stats_dir = root / "project_summary" / "statistics" / "profiles"
    replicate_path, system_path = (
        stats_dir / "replicate_profiles.parquet", stats_dir / "system_profiles.parquet"
    )
    if not replicate_path.exists() or not system_path.exists():
        return
    reps, systems_df = pd.read_parquet(replicate_path), pd.read_parquet(system_path)
    comparisons_path = stats_dir / "system_comparisons.parquet"
    needs_contrasts = (selected is None or any(
        item.startswith("profiles:") and item.rsplit(":", 1)[-1] in
        {"system_difference", "comparison_table"} for item in selected))
    needs_bootstrap_tables = (selected is None or any(
        item.startswith("profiles:") and item.endswith(":bootstrap_table")
        for item in selected))
    comparisons = (pd.read_parquet(comparisons_path)
                   if needs_contrasts and comparisons_path.exists() else pd.DataFrame())
    global_path = stats_dir / "system_profile_statistics.parquet"
    global_stats = (pd.read_parquet(global_path) if needs_bootstrap_tables
                    and global_path.exists() else pd.DataFrame())
    rep_summary_path = stats_dir / "replicate_profile_statistics.parquet"
    rep_summary = (pd.read_parquet(rep_summary_path) if needs_bootstrap_tables
                   and rep_summary_path.exists() else pd.DataFrame())
    systems = list(system_repouts)

    for profile_id in reps["profile"].drop_duplicates():
        stem = _safe(profile_id)
        if selected is not None and not any(
                identifier.startswith(f"profiles:{stem}:") for identifier in selected):
            continue
        data = reps[reps.profile == profile_id]
        sysdata = systems_df[systems_df.profile == profile_id]
        if not (np.isfinite(pd.to_numeric(data["value"], errors="coerce")).any()
                or (not sysdata.empty and np.isfinite(
                    pd.to_numeric(sysdata["mean"], errors="coerce")).any())):
            continue
        xlabel, ylabel, profile_title = _profile_meta(str(profile_id))
        stem = _safe(profile_id)
        fig, axes = _subplot_grid(len(systems))
        distributions, labels, colors_list = [], [], []
        for ax, system in zip(axes, systems):
            colors = _replicate_colors(cfg, system)
            for rep_id, subset in data[data.system == system].groupby("replicate", sort=False):
                subset = subset.sort_values("position")
                if not np.isfinite(pd.to_numeric(subset.value, errors="coerce")).any():
                    continue
                ax.plot(
                    subset.position,
                    subset.value,
                    lw=1.0,
                    alpha=0.8,
                    color=colors.get(str(rep_id)),
                    label=replicate_label(cfg, system, str(rep_id))
                )
                values = finite(subset.value)
                if len(values):
                    distributions.append(values)
                    labels.append(
                        f"{system_label(cfg, system)}\n"
                        f"{replicate_label(cfg, system, str(rep_id))}"
                    )
                    colors_list.append(colors.get(str(rep_id)))
            mean = sysdata[sysdata.system == system].sort_values("position")
            if not mean.empty and np.isfinite(
                    pd.to_numeric(mean["mean"], errors="coerce")).any():
                base = _system_colors(cfg)[system]
                ax.plot(mean.position, mean["mean"], lw=2.0, color=base, label="System mean")
                ax.fill_between(
                    mean.position.to_numpy(float),
                    mean.mean_ci_low.to_numpy(float),
                    mean.mean_ci_high.to_numpy(float),
                    color=base,
                    alpha=0.18
                )
            ax.set_title(system_label(cfg, system))
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if ax.has_data():
                ax.legend(frameon=False)
            else:
                ax.text(.5, .5, "No finite data", ha="center", va="center",
                        transform=ax.transAxes, color="0.45")
        synchronize_y_axes(axes[:len(systems)])
        for ax in axes[len(systems):]:
            ax.remove()
        fig.suptitle(profile_title)
        _save_selected(fig, out / f"{stem}__profile.png", cfg, f"profiles:{stem}:profile")

        if distributions:
            fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(distributions)), 5))
            parts = ax.violinplot(distributions, showmedians=True)
            for body, color in zip(parts["bodies"], colors_list):
                body.set_facecolor(color)
                body.set_alpha(0.6)
            ax.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=25, ha="right")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{profile_title}: value distribution")
            _save_selected(fig, out / f"{stem}__distribution.png", cfg, f"profiles:{stem}:distribution")

        comp = (
            comparisons[
                comparisons.profile == profile_id
            ] if not comparisons.empty else pd.DataFrame()
        )
        if not comp.empty and np.isfinite(pd.to_numeric(
                comp["difference_system2_minus_system1"], errors="coerce")).any():
            fig, ax = plt.subplots(figsize=(10, 5))
            for (s1, s2), subset in comp.groupby(["system_1", "system_2"], sort=False):
                subset = subset.sort_values("position")
                if not np.isfinite(pd.to_numeric(
                        subset.difference_system2_minus_system1,
                        errors="coerce")).any():
                    continue
                ax.plot(
                    subset.position,
                    subset.difference_system2_minus_system1,
                    label=(
                        f"{system_label(cfg, s2)} - "
                        f"{system_label(cfg, s1)}"
                    ),
                )
                ax.fill_between(
                    subset.position.to_numpy(float),
                    subset.difference_ci_low.to_numpy(float),
                    subset.difference_ci_high.to_numpy(float),
                    alpha=0.18
                )
            ax.axhline(0, color="0.4", lw=0.8)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(f"Difference in {ylabel}")
            ax.set_title(f"{profile_title}: system differences")
            ax.legend(frameon=False)
            _save_selected(fig, out / f"{stem}__system_difference.png", cfg, f"profiles:{stem}:system_difference")
            if not cfg.get("plots", {}).get("table_reports_pdf", False):
                rank_col = "q_value" if "q_value" in comp else "p_value"
                ranked = comp.sort_values(rank_col, na_position="last")
                shown = pd.DataFrame({
                    "comparison": ranked.apply(
                        lambda r: (
                            f"{system_label(cfg, r.system_2)} - "
                            f"{system_label(cfg, r.system_1)}"
                        ),
                        axis=1,
                    ),
                    "position": ranked.position_label,
                    "delta": ranked.difference_system2_minus_system1,
                    "CI low": ranked.difference_ci_low, "CI high": ranked.difference_ci_high,
                    "P(delta>0)": ranked.probability_difference_gt_0, "p": ranked.p_value,
                    "q": ranked["q_value"] if "q_value" in ranked else np.nan,
                })
                _table_figure(
                    f"{profile_title}: pointwise comparison statistics",
                    shown,
                    out / f"{stem}__comparison_table.png",
                    cfg,
                    max_rows=16,
                )

        if not cfg.get("plots", {}).get("table_reports_pdf", False):
            table_rows = []
            rep_subset = (
                rep_summary[
                    rep_summary.profile == profile_id
                ] if not rep_summary.empty else pd.DataFrame()
            )
            for _, row in rep_subset.iterrows():
                table_rows.append({
                    "level": (
                        f"{system_label(cfg, row.system)}/"
                        f"{replicate_label(cfg, row.system, row.replicate)}"
                    ),
                    "mean profile": row["mean"],
                    "profile CI": f"{row.mean_profile_ci_low:.3g}-{row.mean_profile_ci_high:.3g}",
                    "tau(fr)": row.tau_int_frames, "tau(ps)": row.tau_int_ps, "N_eff": row.effective_n,
                    "block(fr)": int(row.block_length_frames),
                    "blocks/pseudo": int(row.blocks_per_pseudostudy),
                })
            global_subset = (
                global_stats[
                    global_stats.profile == profile_id
                ] if not global_stats.empty else pd.DataFrame()
            )
            for _, row in global_subset.iterrows():
                table_rows.append({
                    "level": f"{system_label(cfg, row.system)} (system)",
                    "mean profile": row.mean_profile_value,
                    "profile CI": f"{row.mean_profile_ci_low:.3g}-{row.mean_profile_ci_high:.3g}",
                    "tau(fr)": (
                        f"{row.tau_int_frames_median:.3g} "
                        f"[{row.tau_int_frames_min:.3g},"
                        f"{row.tau_int_frames_max:.3g}]"
                    ),
                    "tau(ps)": "replicate-specific", "N_eff": np.nan,
                    "block(fr)": "replicate-specific",
                    "blocks/pseudo": (
                        f"{row.blocks_per_pseudostudy_median:.3g} "
                        f"[{int(row.blocks_per_pseudostudy_min)},"
                        f"{int(row.blocks_per_pseudostudy_max)}]"
                    ),
                })
            _table_figure(
                f"{profile_title}: hierarchical bootstrap statistics",
                pd.DataFrame(table_rows),
                out / f"{stem}__bootstrap_table.png",
                cfg,
            )


def _plot_residue_distance_heatmaps(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]):
    from ..config import analysis_runs, profile_dir
    for run in analysis_runs(cfg, "residue_membrane_distance"):
        run_out = (root / "project_summary" / "figures" / "profiles"
                   / _safe(f"residue_membrane_distance__{run.run_id}"))
        for system, repouts in system_repouts.items():
            plot_id = f"profiles:{run_out.name}:{_safe(system)}:distance_heatmap"
            if not figure_enabled(cfg, plot_id):
                continue
            available = []
            for repout in repouts:
                directory = profile_dir(repout, "residue_membrane_distance", run.run_id)
                matrix_path, time_path, residue_path = (
                    directory / "distance_matrix_A.npy", directory / "time_ps.npy",
                    directory / "residues.parquet"
                )
                if matrix_path.exists() and time_path.exists() and residue_path.exists():
                    available.append((repout, matrix_path, time_path, residue_path))
            if not available:
                continue
            fig, axes = plt.subplots(
                len(available),
                1,
                figsize=(12, max(4, 3.6 * len(available))),
                squeeze=False
            )
            # Each replicate measures the same membrane distance. Use one
            # color scale and one y/time range for interpretable comparisons.
            configured_vmax = run.spec.get("heatmap_max_A")
            if configured_vmax is None:
                maxima = []
                for _, matrix_path, _, _ in available:
                    matrix_data = np.load(matrix_path, mmap_mode="r")
                    finite = np.isfinite(matrix_data)
                    if finite.any():
                        maxima.append(float(np.nanmax(matrix_data)))
                shared_vmax = max(maxima) if maxima else None
            else:
                shared_vmax = float(configured_vmax)
            image = None
            for ax, (repout, matrix_path, time_path, residue_path) in zip(axes.ravel(), available):
                matrix = np.load(matrix_path, mmap_mode="r")
                times = np.load(time_path) / 1000.0
                residues = pd.read_parquet(residue_path)
                image = ax.imshow(
                    matrix,
                    aspect="auto",
                    origin="lower",
                    extent=[0, matrix.shape[1], float(times[0]), float(times[-1])],
                    vmin=0,
                    vmax=shared_vmax,
                    cmap="viridis"
                )
                step = max(1, matrix.shape[1] // 15)
                ticks = np.arange(0, matrix.shape[1], step)
                labels = residues.iloc[ticks]["residue"].astype(str).tolist()
                ax.set_xticks(ticks + 0.5, labels, rotation=45, ha="right")
                ax.set_ylabel("Time (ns)")
                ax.set_title(replicate_label(cfg, system, repout.name))
            synchronize_y_axes(axes)
            axes.ravel()[-1].set_xlabel("Residue")
            if image is not None:
                fig.colorbar(
                    image,
                    ax=axes.ravel().tolist(),
                    label="Minimum membrane distance (A)",
                    pad=0.02
                )
            fig.suptitle(f"{system_label(cfg, system)} - residue/membrane distance ({run.label})")
            _save(fig, run_out / f"{_safe(system)}__distance_heatmap.png", cfg)


def _interaction_key(data: pd.DataFrame) -> pd.Series:
    metric_cols = {
        "n_events",
        "mean_lifetime_ps",
        "median_lifetime_ps",
        "max_lifetime_ps",
        "total_lifetime_ps",
        "frames_present",
        "occupancy"
    }
    keys = [col for col in data.columns if col not in metric_cols]
    return data[keys].astype(str).agg(" | ".join, axis=1)


def _lipid_hbond_pair_matrix(tables):
    """Mean simultaneous atom-level bonds by species pair, independently per replica.

    Summing individual bond occupancies is the mean count per sampled frame,
    not a probability. A species is known present if it participates in any
    saved bond in that replica. Without a population inventory, absence and
    unobserved species cannot be distinguished and remain masked.
    """
    pairs = sorted({
        str(pair) for table in tables if table is not None
        for pair in table["interaction_class"].dropna().unique()
    })
    matrix = np.full((len(pairs), len(tables)), np.nan)
    for j, table in enumerate(tables):
        if table is None or table.empty:
            continue
        counts = table.groupby("interaction_class")["occupancy"].sum()
        species = {s for pair in counts.index for s in str(pair).split("__")}
        for i, pair in enumerate(pairs):
            if set(pair.split("__")).issubset(species):
                matrix[i, j] = counts.get(pair, 0.0)
    return pairs, matrix


def _plot_lipid_hbond_pairs(cfg, kind_out, columns, tables):
    pairs, matrix = _lipid_hbond_pair_matrix(tables)
    if not pairs:
        return
    fig, ax = plt.subplots(figsize=(
        max(8, 1.2 * len(columns)), max(4, 0.3 * len(pairs) + 2)
    ))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("lightgray")
    image = ax.imshow(np.ma.masked_invalid(matrix), aspect="auto", vmin=0, cmap=cmap)
    ax.set_xticks(np.arange(len(columns)), columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(pairs)), [pair.replace("__", " – ") for pair in pairs])
    ax.set_title("Lipid hydrogen bonds by species pair")
    ax.set_xlabel("Gray: species not observed or data unavailable")
    fig.colorbar(image, ax=ax, label="Mean H-bonds per sampled frame")
    _save_selected(fig, kind_out / "occupancy_heatmap.png", cfg,
                   "interactions:lipid_hydrogen_bonds:occupancy_heatmap")


def plot_interactions(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]):
    if not plot_enabled(cfg, 'interactions'):
        return
    stats_allowed = current_statistics(cfg, root, 'interactions')
    out = ensure_dir(root / "project_summary" / "figures" / "interactions")
    kinds = set()
    tables = {}
    selected = cfg.get("plots", {}).get("include")
    selected_kinds = (None if selected is None else {
        identifier.split(":")[1] for identifier in selected
        if identifier.startswith("interactions:")
    })
    if selected is None:
        table_names = ("presence", "occupancy", "lifetimes", "summary")
    else:
        needed = set()
        for item in selected:
            if not item.startswith("interactions:"):
                continue
            variant = item.rsplit(":", 1)[-1]
            if variant == "occupancy_heatmap":
                needed.add("occupancy")
            elif variant == "lifetime_distribution":
                needed.add("lifetimes")
        table_names = tuple(sorted(needed))
        kinds.update((selected_kinds or set()) - set(
            cfg.get("plots", {}).get("interaction_exclude_kinds", [])))
    for system, reps in system_repouts.items():
        for repout in reps:
            tables[(system, repout.name)] = (
                interaction_tables(repout,
                    exclude_kinds=set(cfg.get('plots',{}).get('interaction_exclude_kinds',[])),
                    include_kinds=selected_kinds, table_names=table_names)
                if table_names else {})
            kinds.update(tables[(system, repout.name)])
    stats_path = (
        root
        / "project_summary"
        / "statistics"
        / "interactions"
        / "system_interaction_metrics.parquet"
    )
    needs_stats = selected is None or any(
        item.startswith("interactions:") and item.rsplit(":", 1)[-1] in
        {"bootstrap_table", "comparison_table"} for item in selected)
    needs_comp = selected is None or any(
        item.startswith("interactions:") and item.endswith(":comparison_table")
        for item in selected)
    stats = (pd.read_parquet(stats_path) if stats_path.exists()
             and stats_allowed and needs_stats else pd.DataFrame())
    comp_path = (
        root
        / "project_summary"
        / "statistics"
        / "interactions"
        / "system_interaction_comparisons.parquet"
    )
    comp_stats = (pd.read_parquet(comp_path) if comp_path.exists()
                  and stats_allowed and needs_comp else pd.DataFrame())

    for kind in sorted(kinds):
        if selected is not None and not any(
                identifier.startswith(f"interactions:{_safe(kind)}:") for identifier in selected):
            continue
        kind_out = ensure_dir(out / _safe(kind))
        columns, occupancy_series, lifetime_sets, lifetime_labels = [], [], [], []
        lipid_pair_tables = []
        for system, reps in system_repouts.items():
            for repout in reps:
                kind_tables = tables[(system, repout.name)].get(kind, {})
                occ = kind_tables.get("occupancy", pd.DataFrame())
                columns.append(
                    f"{system_label(cfg, system)}\n"
                    f"{replicate_label(cfg, system, repout.name)}"
                )
                if kind == "lipid_hydrogen_bonds":
                    lipid_pair_tables.append(kind_tables.get("occupancy"))
                elif occ.empty:
                    occupancy_series.append(pd.Series(dtype=float))
                else:
                    occ = occ.copy()
                    occ["interaction"] = _interaction_key(occ)
                    occupancy_series.append(occ.set_index("interaction")["occupancy"])
                life = kind_tables.get("lifetimes", pd.DataFrame())
                if not life.empty and "lifetime_ps" in life:
                    if "completed" in life:
                        life = life[life["completed"]]
                    values = finite(life["lifetime_ps"]) / 1000.0
                    if len(values):
                        lifetime_sets.append(values)
                        lifetime_labels.append(columns[-1])

        if kind == "lipid_hydrogen_bonds" and lipid_pair_tables:
            _plot_lipid_hbond_pairs(cfg, kind_out, columns, lipid_pair_tables)

        union = (
            sorted(
                set().union(*(set(series.index) for series in occupancy_series))
            ) if occupancy_series else []
        )
        if union:
            matrix = np.column_stack([series.reindex(
                union,
                fill_value=0
            ).to_numpy(float) for series in occupancy_series])
            keep = np.max(matrix, axis=1) >= float(cfg.get("plots", {}).get("min_occupancy", 0.05))
            matrix, labels = matrix[keep], np.asarray(union, object)[keep]
            if len(labels):
                order = np.argsort(-np.max(matrix, axis=1))
                matrix, labels = matrix[order], labels[order]
                from .templates import heatmap_pages
                identifier = f"interactions:{_safe(kind)}:occupancy_heatmap"
                if figure_enabled(cfg, identifier):
                    for fig, ax, rs, cs, page in heatmap_pages(matrix, labels, columns,
                            title=f"{kind.replace('_', ' ').title()} occupancy",
                            units="Occupancy fraction", vmin=0, vmax=1):
                        path = figure_path(cfg, kind_out / "occupancy_heatmap.png")
                        if len(labels) > 36 or len(columns) > 24:
                            path = path.with_stem(path.stem + f"@page_{page:03d}")
                        _save(fig, path, cfg)

        if lifetime_sets:
            fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(lifetime_sets)), 5))
            ax.violinplot(lifetime_sets, showmedians=True)
            ax.set_xticks(
                np.arange(1, len(lifetime_labels) + 1),
                lifetime_labels,
                rotation=25,
                ha="right"
            )
            ax.set_ylabel("Lifetime (ns)")
            ax.set_title(f"{kind.replace('_', ' ').title()} lifetime distribution")
            _save_selected(fig, kind_out / "lifetime_distribution.png", cfg,
                           f"interactions:{_safe(kind)}:lifetime_distribution")

        if not cfg.get("plots", {}).get("table_reports_pdf", False) and not stats.empty:
            subset = stats[stats.interaction_type == kind]
            for metric in subset.metric.drop_duplicates():
                metric_data = subset[subset.metric == metric].sort_values(
                    "mean",
                    ascending=False
                )
                if metric_data.empty:
                    continue
                shown = pd.DataFrame({
                    "system": metric_data.system.map(lambda x: system_label(cfg, x)),
                    "interaction": metric_data.interaction.map(
                        lambda x: x if len(x) < 55 else x[:52] + "..."
                    ),
                    "mean": metric_data["mean"], "CI low": metric_data.mean_ci_low,
                    "CI high": metric_data.mean_ci_high,
                    "n reps": metric_data.n_replicates_total,
                    "tau(fr)": metric_data.apply(
                        lambda r: (
                            f"{r.tau_int_frames_median:.3g} "
                            f"[{r.tau_int_frames_min:.3g},"
                            f"{r.tau_int_frames_max:.3g}]"
                        ),
                        axis=1,
                    ),
                    "tau(ps)": metric_data.apply(
                        lambda r: (
                            f"{r.tau_int_ps_median:.3g} "
                            f"[{r.tau_int_ps_min:.3g},{r.tau_int_ps_max:.3g}]"
                        ),
                        axis=1,
                    ),
                    "block(fr)": metric_data.apply(
                        lambda r: (
                            f"{r.block_length_frames_median:.3g} "
                            f"[{int(r.block_length_frames_min)},"
                            f"{int(r.block_length_frames_max)}]"
                        ),
                        axis=1,
                    ),
                    "blocks/pseudo": metric_data.apply(
                        lambda r: (
                            f"{r.blocks_per_pseudostudy_median:.3g} "
                            f"[{int(r.blocks_per_pseudostudy_min)},"
                            f"{int(r.blocks_per_pseudostudy_max)}]"
                        ),
                        axis=1,
                    ),
                })
                _table_figure(
                    f"{kind}: {metric} bootstrap statistics",
                    shown,
                    kind_out / f"{_safe(metric)}__bootstrap_table.png",
                    cfg,
                    max_rows=12,
                )
                if not comp_stats.empty:
                    c = (
                        comp_stats[
                            (comp_stats.interaction_type == kind) & (comp_stats.metric == metric)
                        ].copy()
                    )
                    if not c.empty:
                        rank_col = "q_value" if "q_value" in c else "p_value"
                        c = c.sort_values(rank_col, na_position="last")
                        comparison = pd.DataFrame({
                            "comparison": c.apply(
                        lambda r: (
                            f"{system_label(cfg, r.system_2)} - "
                            f"{system_label(cfg, r.system_1)}"
                        ),
                        axis=1,
                    ),
                            "interaction": c.interaction.map(
                                lambda x: x if len(x) < 45 else x[:42] + "..."
                            ),
                            "delta": c.difference_system2_minus_system1,
                            "CI low": c.difference_ci_low, "CI high": c.difference_ci_high,
                            "P(delta>0)": c.probability_difference_gt_0, "p": c.p_value,
                            "q": c["q_value"] if "q_value" in c else np.nan,
                        })
                        _table_figure(
                            f"{kind}: {metric} comparison statistics",
                            comparison,
                            kind_out / f"{_safe(metric)}__comparison_table.png",
                            cfg,
                            max_rows=12,
                        )


def plot_all(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]):
    if cfg.get("plots", {}).get("include") == []:
        return
    if not cfg.get("plots", {}).get("enabled", True):
        return
    configure_plot_style(cfg)
    from .overview import plot_overviews
    plot_overviews(cfg, root, system_repouts)
    selected = cfg.get('plots', {}).get('include')
    tables_selected = selected is None or any(identifier.endswith((':bootstrap_table', ':comparison_table')) for identifier in selected)
    with table_report_session(root, enabled=cfg.get("plots", {}).get("table_reports_pdf", False) and tables_selected, cfg=cfg):
        if plot_enabled(cfg, "time_series"):
            plot_time_series(cfg, root, system_repouts)
        if plot_enabled(cfg, "profiles"):
            plot_profiles(cfg, root, system_repouts)
            _plot_residue_distance_heatmaps(cfg, root, system_repouts)
            from .residue_response import plot_residue_response
            plot_residue_response(cfg, root)
        if plot_enabled(cfg, "interactions"):
            plot_interactions(cfg, root, system_repouts)

        if cfg.get("plots", {}).get("table_reports_pdf", False):
            from .statistical_tables import render_statistical_reports
            render_statistical_reports(cfg, root)
