"""Condition-centered SDC4 diagnostics using independently preserved replicas.

All inputs are saved replica numerical products. A subplot/trace represents
one replicate; no frame-level pooling, cross-replica correlation, or invented
bootstrap uncertainty occurs here. Missing replicas remain visibly missing.
"""
from __future__ import annotations

from pathlib import Path
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..config import figure_enabled, system_label
from ..core.persistence import archive_output
from .plot_ids import figure_dpi
from .shared_scales import synchronize_y_axes
from .sdc4_exports import export_region_replica_summaries


REGION_METRICS = {
    "lipid_region_mean_contact_fraction": ("contact_fraction",
        "Fraction of available molecules contacting"),
    "lipid_region_mean_local_enrichment": ("local_composition_enrichment",
        "Legacy unmatched local enrichment (comparison only)"),
    "lipid_region_mean_matched_local_contact_fraction": (
        "matched_local_contact_fraction",
        "Fraction contacting within local anchor population"),
    "lipid_region_mean_matched_local_enrichment": (
        "matched_local_composition_enrichment",
        "Matched local contact share / local composition share"),
}
MOIETY_STEM = "lipid_region_moiety_contact_fraction"
MAP_STEMS = ("lipid_density_maps", "lipid_enrichment_maps")
SERIES_STEMS = ("nxip_height_occlusion_joint", "vvhplvpl_proline_ramachandran")
CONDITION_STEMS = frozenset((*REGION_METRICS, MOIETY_STEM, *MAP_STEMS, *SERIES_STEMS))


def _save(fig, cfg, out, system, stem, *, dpi=180):
    destination = out / f"{stem}.png"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=figure_dpi(cfg, f"sdc4:{system}:{stem}", dpi),
                bbox_inches="tight")
    plt.close(fig)


def _tables(reps, filename):
    tables = {}
    for rep in reps:
        path = Path(rep) / "extensions/sdc4_mechanisms" / filename
        if path.is_file():
            tables[Path(rep).name] = pd.read_parquet(path)
    return tables


def _region_label(region: str) -> str:
    """Short display names; saved region identifiers are never altered."""
    if region in {"TM_A", "TM_B"}:
        return region.replace("_", " ")
    motif = re.fullmatch(
        r"motif_(NxIP|acidic_patch|proline_hydrophobic)_(PRO[AB])_(\d+)",
        str(region),
    )
    if motif:
        name, chain, occurrence = motif.groups()
        title = {"NxIP": "NxIP", "acidic_patch": "Acidic",
                 "proline_hydrophobic": "Proline"}[name]
        return f"{title} {chain[-1]} · {occurrence}"
    return "\n".join(textwrap.wrap(str(region).replace("_", " "), 17))


def _archive_uninformative_region_figure(out: Path, stem: str) -> None:
    """Keep old PNGs outside the active project-summary figure inventory."""
    previous = out / f"{stem}.png"
    if previous.is_file():
        archive_output(previous, out.parents[3] / ".archives" /
                       "uninformative_sdc4_region_figures", f"{out.name}_{stem}.png")


def _region_figures(cfg, system, reps, out):
    stems = [name for name in (*REGION_METRICS, MOIETY_STEM)
             if figure_enabled(cfg, f"sdc4:{system}:{name}")]
    if not stems:
        return
    tables = _tables(reps, "lipid_region_specificity.parquet")
    if not tables:
        return
    export_region_replica_summaries(out.parents[3], system, reps, tables)
    names = [Path(rep).name for rep in reps]
    for stem in stems:
        moieties = ("glycan", "ceramide_head", "ceramide_tail") if stem == MOIETY_STEM else ("any",)
        column, ylabel = (("contact_fraction", "Contact fraction")
                          if stem == MOIETY_STEM else REGION_METRICS[stem])

        # Calculate displayed replicate means once; never pool frames/replicas.
        panels = {}
        regions = set()
        active_series = set()
        omitted_zero_series = 0
        for repname in names:
            tab = tables.get(repname)
            if tab is None or column not in tab or not {
                    "region", "species", "leaflet", "moiety"}.issubset(tab):
                continue
            for moiety in moieties:
                panel = tab.loc[tab.moiety.eq(moiety)].copy()
                panel[column] = pd.to_numeric(panel[column], errors="coerce")
                panel = panel.loc[np.isfinite(panel[column])]
                if panel.empty:
                    continue
                grouped = panel.groupby(
                    ["region", "species", "leaflet"], as_index=False
                ).agg(value=(column, "mean"))
                panels[(repname, moiety)] = grouped
                for (species, leaflet), series in grouped.groupby(
                        ["species", "leaflet"], sort=True):
                    if np.any(series.value.to_numpy(float) != 0):
                        active_series.add((str(species), str(leaflet)))
                        regions.update(series.region.astype(str))
                    else:
                        omitted_zero_series += 1

        # The all-zero region example is uninformative visually but remains
        # fully available in persisted Parquet and per-replica summary CSV.
        if not active_series:
            _archive_uninformative_region_figure(out, stem)
            continue

        region_order = sorted(regions)
        positions = {region: i for i, region in enumerate(region_order)}
        series_order = sorted(active_series)
        palette = plt.get_cmap("tab20")
        colors = {key: palette(i % 20) for i, key in enumerate(series_order)}
        # One readable legend for the entire condition; none inside panels.
        # Separate label space from data space instead of tight_layout's
        # repeated legend/tick compression.
        nrows, ncols = len(names), len(moieties)
        width = 16.5 if ncols == 1 else max(19., 5.0 * ncols + 5.)
        height = max(5.5, 3.2 * nrows + 1.6)
        fig, axes = plt.subplots(nrows, ncols, figsize=(width, height),
                                 squeeze=False, sharex=True, sharey=True)
        for i, repname in enumerate(names):
            for j, moiety in enumerate(moieties):
                ax = axes[i, j]
                ax.set_title(f"{repname} / {moiety}", fontsize=11, pad=8)
                grouped = panels.get((repname, moiety))
                if grouped is None:
                    ax.text(.5, .5, "Not measured", transform=ax.transAxes,
                            ha="center", va="center", color="0.45")
                    continue
                shown = False
                for (species, leaflet), series in grouped.groupby(
                        ["species", "leaflet"], sort=True):
                    key = (str(species), str(leaflet))
                    if key not in active_series or not np.any(
                            series.value.to_numpy(float) != 0):
                        continue
                    x = [positions[str(region)] for region in series.region]
                    ax.plot(x, series.value.to_numpy(float), ".-", ms=4,
                            lw=1.25, color=colors[key])
                    shown = True
                if not shown:
                    ax.text(.5, .5, "Only zero-valued series", transform=ax.transAxes,
                            ha="center", va="center", color="0.45")
                if "enrichment" in column:
                    ax.axhline(1, color="0.5", lw=.8, ls="--")
                ax.grid(axis="y", alpha=.2, lw=.6)

        for ax in axes.flat:
            ax.set_xlim(-.35, len(region_order) - .65 if len(region_order) > 1 else .35)
            ax.set_xticks(range(len(region_order)))
        for ax in axes[-1]:
            ax.set_xticklabels([_region_label(region) for region in region_order],
                               fontsize=9, rotation=0)
            ax.set_xlabel("Protein region", labelpad=10)
        if column.endswith("fraction"):
            for ax in axes.flat:
                ax.set_ylim(-.03, 1.03)
        else:
            synchronize_y_axes(axes)

        fig.supylabel(ylabel, x=.015, fontsize=12)
        found = sum(name in tables for name in names)
        title = stem.removeprefix("lipid_region_").replace("_", " ")
        fig.suptitle(f"{system_label(cfg, system)} · lipid region {title}"
                     f" ({found}/{nrows} replicas)", fontsize=16, y=.97)
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], color=colors[key], marker=".", lw=1.3,
                          label=f"{key[0]} / {key[1]}") for key in series_order]
        fig.legend(handles=handles, loc="center left",
                   bbox_to_anchor=(.78, .53), fontsize=9, ncol=2,
                   title="Lipid / leaflet", title_fontsize=10, frameon=False)
        if omitted_zero_series:
            fig.text(.5, .025,
                     f"{omitted_zero_series} zero-only replicate series omitted from display; "
                     "saved numerical products retained", ha="center", fontsize=9)
        fig.subplots_adjust(left=.09, right=.76, top=.88, bottom=.12,
                            hspace=.42, wspace=.22)
        _save(fig, cfg, out, system, stem)


def _map_figures(cfg, system, reps, out):
    active = [stem for stem in MAP_STEMS
              if figure_enabled(cfg, f"sdc4:{system}:{stem}")]
    if not active:
        return
    maps = {}
    for rep in reps:
        path = Path(rep) / "extensions/sdc4_mechanisms/lipid_density_maps.npz"
        if not path.is_file():
            continue
        with np.load(path, allow_pickle=False) as data:
            required = {"species", "leaflet", "enrichment_fold",
                        "density_per_nm2", "x_edges_A", "y_edges_A"}
            if not required.issubset(data.files):
                raise ValueError(f"SDC4 map is missing required arrays: {path}")
            species, leaflets = data["species"].astype(str), data["leaflet"].astype(str)
            entries = list(zip(species, leaflets))
            if len(set(entries)) != len(entries):
                raise ValueError(f"Duplicate SDC4 species/leaflet in {path}")
            maps[Path(rep).name] = {
                "keys": entries, "x": data["x_edges_A"].copy(),
                "y": data["y_edges_A"].copy(),
                "density": data["density_per_nm2"].copy(),
                "enrichment": data["enrichment_fold"].copy(),
            }
    if not maps:
        return
    names = [Path(rep).name for rep in reps]
    keys = sorted({key for item in maps.values() for key in item["keys"]})
    for stem in active:
        field = "density" if stem == "lipid_density_maps" else "enrichment"
        values = [
            array for item in maps.values() for array in item[field]
            if np.isfinite(array).any()
        ]
        if not values:
            continue
        limit = max(1e-12, float(max(np.nanmax(array) for array in values)))
        fig, axes = plt.subplots(len(names), len(keys),
                                 figsize=(4.2 * len(keys), 3.5 * len(names)),
                                 squeeze=False, constrained_layout=True)
        image = None
        for i, repname in enumerate(names):
            item = maps.get(repname)
            for j, key in enumerate(keys):
                ax = axes[i, j]
                ax.set_title(f"{repname} / {key[0]} / {key[1]}", fontsize=9)
                if item is None or key not in item["keys"]:
                    ax.text(.5, .5, "Not measured", transform=ax.transAxes,
                            ha="center", va="center")
                    continue
                index = item["keys"].index(key)
                values2d = item[field][index]
                if values2d.shape != (len(item["x"]) - 1, len(item["y"]) - 1):
                    raise ValueError(f"SDC4 map array/edges mismatch for {system}/{repname}/{key}")
                image = ax.pcolormesh(item["x"], item["y"], values2d.T, shading="auto",
                                      vmin=0, vmax=limit, cmap="viridis")
                ax.set_aspect("equal")
                ax.set_xlabel("TM-centered x (Å), A → B")
                ax.set_ylabel("TM-centered y (Å)")
        if image is not None:
            fig.colorbar(image, ax=list(axes.flat),
                         label=("Lipids / nm²" if field == "density"
                                else "Local / bulk lipid density"))
            found = sum(name in maps for name in names)
            fig.suptitle(f"{system}: {field} maps — {found}/{len(names)} independent replicas")
            _save(fig, cfg, out, system, stem)
        else:
            plt.close(fig)


def _series_figures(cfg, system, reps, out):
    active = [name for name in SERIES_STEMS
              if figure_enabled(cfg, f"sdc4:{system}:{name}")]
    if not active:
        return
    tables = _tables(reps, "timeseries.parquet")
    if not tables:
        return
    names = [Path(rep).name for rep in reps]
    if "nxip_height_occlusion_joint" in active:
        fig, axs = plt.subplots(1, 2, figsize=(10, 4), squeeze=False,
                                constrained_layout=True)
        present = False
        for ax, chain in zip(axs.flat, ("A", "B")):
            xcol = f"sdc4_ectodomain_{chain}_height_A"
            for rep in names:
                tab = tables.get(rep)
                if tab is None or xcol not in tab:
                    continue
                matches = [col for col in tab
                           if col.startswith(f"sdc4_motif_NxIP_PRO{chain}_")
                           and col.endswith("_steric_occlusion_fraction")]
                if len(matches) != 1:
                    continue
                x, y = tab[xcol].to_numpy(float), tab[matches[0]].to_numpy(float)
                valid = np.isfinite(x) & np.isfinite(y)
                if valid.any():
                    ax.scatter(x[valid], y[valid], s=7, alpha=.23, rasterized=True,
                               label=rep)
                    present = True
            ax.set_xlabel(f"Chain {chain} ectodomain height (Å)")
            ax.set_ylabel("NxIP steric occlusion fraction")
            ax.set_ylim(-.03, 1.03)
            if ax.has_data():
                ax.legend(fontsize=8, frameon=False)
            else:
                ax.text(.5, .5, "NxIP mapping unavailable",
                        transform=ax.transAxes, ha="center", va="center")
        if present:
            fig.suptitle(f"{system}: NxIP–height, descriptive within-replica points only")
            _save(fig, cfg, out, system, "nxip_height_occlusion_joint", dpi=190)
        else:
            plt.close(fig)

    if "vvhplvpl_proline_ramachandran" in active:
        prefixes = sorted({
            col.removesuffix("_phi_deg")
            for tab in tables.values() for col in tab
            if col.startswith("sdc4_proline_") and col.endswith("_phi_deg")
            and col.removesuffix("_phi_deg") + "_psi_deg" in tab
        })
        if not prefixes:
            return
        fig, axs = plt.subplots(1, len(prefixes), squeeze=False,
                                figsize=(3.5 * len(prefixes), 3.8),
                                constrained_layout=True)
        for ax, prefix in zip(axs.flat, prefixes):
            for rep in names:
                tab = tables.get(rep)
                if tab is None or (prefix + "_phi_deg") not in tab:
                    continue
                phi = tab[prefix + "_phi_deg"].to_numpy(float)
                psi = tab[prefix + "_psi_deg"].to_numpy(float)
                finite = np.isfinite(phi) & np.isfinite(psi)
                ax.scatter(phi[finite], psi[finite], s=6, alpha=.20,
                           rasterized=True, label=rep)
            ax.set(xlim=(-180, 180), ylim=(-180, 180),
                   xlabel="φ (degrees)", ylabel="ψ (degrees)",
                   title=prefix.removeprefix("sdc4_proline_"))
            ax.set_aspect("equal", adjustable="box")
            if ax.has_data():
                ax.legend(fontsize=8, frameon=False)
        fig.suptitle(f"{system}: proline Ramachandran distributions by independent replica")
        _save(fig, cfg, out, system, "vvhplvpl_proline_ramachandran")


def plot_sdc4_condition_figures(cfg, root, system_repouts):
    """One condition figure per metric; no PNGs are stored under replica dirs."""
    for system, reps in system_repouts.items():
        if not reps:
            continue
        if not any(figure_enabled(cfg, f"sdc4:{system}:{stem}")
                   for stem in CONDITION_STEMS):
            continue
        output = Path(root) / "project_summary/figures/sdc4" / system
        _region_figures(cfg, system, reps, output)
        _map_figures(cfg, system, reps, output)
        _series_figures(cfg, system, reps, output)
