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

from .config import ensure_dir, replicate_label, replicate_specs, system_label
from .interactions import interaction_tables
from .statistics import collect_time_series, finite


METRIC_LABELS = {
    "rmsd_A": ("RMSD", "A"),
    "radius_of_gyration_A": ("Radius of gyration", "A"),
    "ligand_rmsd_A": ("Ligand RMSD", "A"),
    "ligand_min_distance_A": ("Ligand minimum distance", "A"),
    "ligand_contact_residue_count": ("Ligand-contacting residues", "Count"),
    "membrane_thickness_A": ("Membrane thickness", "A"),
    "area_per_lipid_A2": ("Area per lipid", "A^2"),
    "area_per_lipid_upper_A2": ("Upper-leaflet area per lipid", "A^2"),
    "area_per_lipid_lower_A2": ("Lower-leaflet area per lipid", "A^2"),
    "membrane_msd_lateral_A2": ("Membrane lateral MSD", "A^2"),
    "total_sasa_nm2": ("Total SASA", "nm^2"),
    "helix_fraction": ("Helix fraction", "Fraction"),
    "sheet_fraction": ("Sheet fraction", "Fraction"),
    "coil_fraction": ("Coil fraction", "Fraction"),
}


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


def _safe(text) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in str(text)).strip("_")


def _save(fig, path: Path, cfg: dict):
    ensure_dir(path.parent)
    fig.tight_layout()
    fig.savefig(path, dpi=int(cfg.get("plots", {}).get("dpi", 250)), bbox_inches="tight")
    plt.close(fig)


def _system_colors(cfg: dict) -> dict[str, tuple]:
    cmap = plt.get_cmap("tab10")
    return{
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
    return{
        rep["id"]: colorsys.hls_to_rgb(h, float(level), max(0.6, sat)) for rep,
        level in zip(reps, levels)
    }


def _metric_meta(metric: str):
    if metric in METRIC_LABELS:
        return METRIC_LABELS[metric]
    if metric == "lipid_order_mean":
        return "Mean lipid acyl-chain order", "-S_CH"
    if metric.startswith("lipid_order_") and metric.endswith("_mean"):
        label = metric[len("lipid_order_"):-len("_mean")].replace("_", " ")
        return f"Mean lipid order: {label}", "-S_CH"
    if metric.endswith("_count"):
        return metric.replace("_", " ").title(), "Count"
    if metric.endswith("_A"):
        return metric[:-2].replace("_", " ").title(), "A"
    if metric.endswith("_A2"):
        return metric[:-3].replace("_", " ").title(), "A^2"
    if metric.endswith("_fraction"):
        return metric.replace("_", " ").title(), "Fraction"
    return metric.replace("_", " ").title(), metric


def _profile_meta(profile_id: str) -> tuple[str, str, str]:
    if profile_id.startswith("lipid_order__"):
        parts = profile_id.split("__")
        if len(parts) == 2:
            title = f"Lipid order: {parts[1]}"
        else:
            title = f"Lipid order: {parts[1]} {parts[2]}"
        return "Acyl-chain carbon position", "-S_CH", title
    return "Profile position", "Value", profile_id.replace("__", " * ").replace("_", " ")


def _subplot_grid(n: int):
    cols = min(3, max(1, n))
    rows = int(math.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4.2 * rows), squeeze=False)
    return fig, axes.ravel()


def _table_figure(title: str, dataframe: pd.DataFrame, path: Path, cfg: dict, max_rows: int = 16):
    # Statistical tables are figures on purpose: they travel with the
    # plot set and are easy to drop directly into slides/manuscripts.
    if dataframe.empty:
        return
    shown = dataframe.head(max_rows).copy()
    for col in shown.columns:
        if pd.api.types.is_float_dtype(shown[col]):
            shown[col] = shown[col].map(lambda x: "" if pd.isna(x) else f"{x:.3g}")
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
    _save(fig, path, cfg)


def plot_time_series(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]):
    if not cfg.get("plots", {}).get("enabled", True):
        return
    out = ensure_dir(root / "project_summary" / "figures" / "time_series")
    systems = list(system_repouts)
    keys = set()
    cache = {}
    for system, reps in system_repouts.items():
        for repout in reps:
            cache[(system, repout.name)] = collect_time_series(repout)
            keys.update(cache[(system, repout.name)].keys())
    rep_stats_path = (
        root / "project_summary" / "statistics" / "time_series" / "replicate_statistics.parquet"
    )
    sys_stats_path = (
        root / "project_summary" / "statistics" / "time_series" / "system_statistics.parquet"
    )
    rep_stats = pd.read_parquet(rep_stats_path) if rep_stats_path.exists() else pd.DataFrame()
    sys_stats = pd.read_parquet(sys_stats_path) if sys_stats_path.exists() else pd.DataFrame()
    comp_stats_path = (
        root / "project_summary" / "statistics" / "time_series" / "system_comparisons.parquet"
    )
    comp_stats = pd.read_parquet(comp_stats_path) if comp_stats_path.exists() else pd.DataFrame()

    for group_id, metric in sorted(keys, key=lambda x: (str(x[0]), x[1])):
        title, ylabel = _metric_meta(metric)
        suffix = f"__group_{group_id}" if group_id is not None else ""
        stem = _safe(metric + suffix)
        fig, axes = _subplot_grid(len(systems))
        distribution_values, distribution_labels, distribution_colors = [], [], []
        for ax, system in zip(axes, systems):
            colors = _replicate_colors(cfg, system)
            for repout in system_repouts[system]:
                table = cache[(system, repout.name)].get((group_id, metric))
                if table is None:
                    continue
                ax.plot(
                    table["time_ps"] / 1000.0,
                    table["value"],
                    lw=1.1,
                    color=colors.get(repout.name),
                    label=replicate_label(cfg, system, repout.name)
                )
                values = finite(table["value"])
                if len(values):
                    distribution_values.append(values)
                    distribution_labels.append(
                        f"{system_label(cfg, system)}\n"
                        f"{replicate_label(cfg, system, repout.name)}"
                    )
                    distribution_colors.append(colors.get(repout.name))
            ax.set_title(system_label(cfg, system))
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel(ylabel)
            ax.legend(frameon=False)
        for ax in axes[len(systems):]:
            ax.remove()
        fig.suptitle(title + (f" - {group_id}" if group_id is not None else ""))
        _save(fig, out / f"{stem}__timeseries.png", cfg)

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
            _save(fig, out / f"{stem}__distribution.png", cfg)

        if cfg.get("plots", {}).get("autocorrelation", True) and not rep_stats.empty:
            acf_path = (
                root / "project_summary" / "statistics" / "time_series" / "autocorrelation.parquet"
            )
            if acf_path.exists():
                acf = pd.read_parquet(acf_path)
                mask = (acf.metric == metric) & (acf.group.fillna("") == (group_id or ""))
                subset = acf[mask]
                if not subset.empty:
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
                        ax.legend(frameon=False)
                    for ax in axes[len(systems):]:
                        ax.remove()
                    fig.suptitle(f"{title} autocorrelation")
                    _save(fig, out / f"{stem}__autocorrelation.png", cfg)

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
    out = ensure_dir(root / "project_summary" / "figures" / "profiles")
    stats_dir = root / "project_summary" / "statistics" / "profiles"
    replicate_path, system_path = (
        stats_dir / "replicate_profiles.parquet", stats_dir / "system_profiles.parquet"
    )
    if not replicate_path.exists() or not system_path.exists():
        return
    reps, systems_df = pd.read_parquet(replicate_path), pd.read_parquet(system_path)
    comparisons_path = stats_dir / "system_comparisons.parquet"
    comparisons = (
        pd.read_parquet(comparisons_path) if comparisons_path.exists() else pd.DataFrame()
    )
    global_path = stats_dir / "system_profile_statistics.parquet"
    global_stats = pd.read_parquet(global_path) if global_path.exists() else pd.DataFrame()
    rep_summary_path = stats_dir / "replicate_profile_statistics.parquet"
    rep_summary = (
        pd.read_parquet(rep_summary_path) if rep_summary_path.exists() else pd.DataFrame()
    )
    systems = list(system_repouts)

    for profile_id in reps["profile"].drop_duplicates():
        data = reps[reps.profile == profile_id]
        sysdata = systems_df[systems_df.profile == profile_id]
        xlabel, ylabel, profile_title = _profile_meta(str(profile_id))
        stem = _safe(profile_id)
        fig, axes = _subplot_grid(len(systems))
        distributions, labels, colors_list = [], [], []
        for ax, system in zip(axes, systems):
            colors = _replicate_colors(cfg, system)
            for rep_id, subset in data[data.system == system].groupby("replicate", sort=False):
                subset = subset.sort_values("position")
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
            if not mean.empty:
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
            ax.legend(frameon=False)
        for ax in axes[len(systems):]:
            ax.remove()
        fig.suptitle(profile_title)
        _save(fig, out / f"{stem}__profile.png", cfg)

        if distributions:
            fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(distributions)), 5))
            parts = ax.violinplot(distributions, showmedians=True)
            for body, color in zip(parts["bodies"], colors_list):
                body.set_facecolor(color)
                body.set_alpha(0.6)
            ax.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=25, ha="right")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{profile_title}: value distribution")
            _save(fig, out / f"{stem}__distribution.png", cfg)

        comp = (
            comparisons[
                comparisons.profile == profile_id
            ] if not comparisons.empty else pd.DataFrame()
        )
        if not comp.empty:
            fig, ax = plt.subplots(figsize=(10, 5))
            for (s1, s2), subset in comp.groupby(["system_1", "system_2"], sort=False):
                subset = subset.sort_values("position")
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
            _save(fig, out / f"{stem}__system_difference.png", cfg)
            rank_col = "q_value" if "q_value" in comp else "p_value"
            ranked = comp.sort_values(rank_col, na_position="last").head(16)
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
    from .config import analysis_runs, profile_dir
    for run in analysis_runs(cfg, "residue_membrane_distance"):
        run_out = (
            ensure_dir(
                root
                / "project_summary"
                / "figures"
                / "profiles"
                / _safe(f"residue_membrane_distance__{run.run_id}")
            )
        )
        for system, repouts in system_repouts.items():
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
            image = None
            for ax, (repout, matrix_path, time_path, residue_path) in zip(axes.ravel(), available):
                matrix = np.load(matrix_path, mmap_mode="r")
                times = np.load(time_path) / 1000.0
                residues = pd.read_parquet(residue_path)
                vmax = run.spec.get("heatmap_max_A")
                image = ax.imshow(
                    matrix,
                    aspect="auto",
                    origin="lower",
                    extent=[0, matrix.shape[1], float(times[0]), float(times[-1])],
                    vmin=0,
                    vmax=float(vmax) if vmax is not None else None,
                    cmap="viridis"
                )
                step = max(1, matrix.shape[1] // 15)
                ticks = np.arange(0, matrix.shape[1], step)
                labels = residues.iloc[ticks]["residue"].astype(str).tolist()
                ax.set_xticks(ticks + 0.5, labels, rotation=45, ha="right")
                ax.set_ylabel("Time (ns)")
                ax.set_title(replicate_label(cfg, system, repout.name))
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
    metric_cols ={
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


def plot_interactions(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]):
    out = ensure_dir(root / "project_summary" / "figures" / "interactions")
    kinds = set()
    tables = {}
    for system, reps in system_repouts.items():
        for repout in reps:
            tables[(system, repout.name)] = interaction_tables(repout)
            kinds.update(tables[(system, repout.name)])
    stats_path = (
        root
        / "project_summary"
        / "statistics"
        / "interactions"
        / "system_interaction_metrics.parquet"
    )
    stats = pd.read_parquet(stats_path) if stats_path.exists() else pd.DataFrame()
    comp_path = (
        root
        / "project_summary"
        / "statistics"
        / "interactions"
        / "system_interaction_comparisons.parquet"
    )
    comp_stats = pd.read_parquet(comp_path) if comp_path.exists() else pd.DataFrame()

    for kind in sorted(kinds):
        kind_out = ensure_dir(out / _safe(kind))
        columns, occupancy_series, lifetime_sets, lifetime_labels = [], [], [], []
        for system, reps in system_repouts.items():
            for repout in reps:
                kind_tables = tables[(system, repout.name)].get(kind, {})
                occ = kind_tables.get("occupancy", pd.DataFrame())
                columns.append(
                    f"{system_label(cfg, system)}\n"
                    f"{replicate_label(cfg, system, repout.name)}"
                )
                if occ.empty:
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
                max_rows = min(40, len(labels))
                matrix, labels = matrix[:max_rows], labels[:max_rows]
                fig, ax = plt.subplots(figsize=(
                    max(8, 1.2 * len(columns)),
                    max(4, 0.3 * len(labels) + 2)
                ))
                image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap="viridis")
                ax.set_xticks(np.arange(len(columns)), columns, rotation=30, ha="right")
                ax.set_yticks(np.arange(len(labels)), labels)
                ax.set_title(f"{kind.replace('_', ' ').title()} occupancy")
                fig.colorbar(image, ax=ax, label="Occupancy")
                _save(fig, kind_out / "occupancy_heatmap.png", cfg)

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
            _save(fig, kind_out / "lifetime_distribution.png", cfg)

        if not stats.empty:
            subset = stats[stats.interaction_type == kind]
            for metric in subset.metric.drop_duplicates():
                metric_data = subset[subset.metric == metric].sort_values(
                    "mean",
                    ascending=False
                ).head(12)
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
                        c = c.sort_values(rank_col, na_position="last").head(12)
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
    if not cfg.get("plots", {}).get("enabled", True):
        return
    configure_plot_style(cfg)
    if cfg.get("plots", {}).get("time_series", True):
        plot_time_series(cfg, root, system_repouts)
    if cfg.get("plots", {}).get("profiles", True):
        plot_profiles(cfg, root, system_repouts)
        _plot_residue_distance_heatmaps(cfg, root, system_repouts)
    if cfg.get("plots", {}).get("interactions", True):
        plot_interactions(cfg, root, system_repouts)
