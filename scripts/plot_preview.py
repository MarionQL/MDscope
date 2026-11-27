#!/usr/bin/env python3
"""Render deterministic fake MD results through the real MDscope plotters.

Usage: python scripts/plot_preview.py --output /tmp/mdscope-plot-preview
No XTC/TPR, statistics calculation, or scientific outputs are required.
Run this *same script* with alternate PYTHONPATHs to compare renderer versions.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

from mdscope.reporting.plotting import plot_all

STAMP = "SYNTHETIC EXAMPLE — NOT MD RESULTS"


def create_fixture(root: Path):
    """Write small, deterministic examples of the *existing* saved-data schema."""
    times = np.linspace(0, 120_000, 61)
    system_repouts = {}
    for system_index, system in enumerate(("condition_A", "condition_B")):
        paths = []
        for rep_index in range(2):
            rep = f"rep{rep_index + 1}"
            repout = root / system / rep
            summary = repout / "summary"
            summary.mkdir(parents=True, exist_ok=False)
            phase = 0.22 * rep_index + 0.4 * system_index
            t = times / 1000.0
            pd.DataFrame({
                "time_ps": times,
                "rmsd_A": (0.5 + 0.15 * system_index
                           + 0.009 * t + 0.16 * np.sin(t / 11 + phase)),
                "radius_of_gyration_A": (14.0 + 0.4 * system_index
                                       + 0.15 * np.sin(t / 17 + phase)),
                "membrane_thickness_A": (31.0 + 1.2 * system_index
                                         + 0.6 * np.cos(t / 15 + phase)),
            }).to_parquet(summary / "timeseries.parquet", index=False)
            by_group = []
            for group_index, group in enumerate(("chain_A", "chain_B")):
                by_group.append(pd.DataFrame({
                    "time_ps": times,
                    "group": group,
                    "rmsd_A": (0.3 + 0.25 * group_index
                               + 0.008 * t + 0.2 * np.sin(t / 10 + phase + group_index)),
                }))
            pd.concat(by_group, ignore_index=True).to_parquet(
                summary / "group_timeseries.parquet", index=False
            )
            paths.append(repout)
        system_repouts[system] = paths
    cfg = {
        "project": {"name": "synthetic_figure_preview"},
        "systems": {
            name: {"label": f"Condition {i + 1}",
                   "replicates": [{"label": f"rep{j+1}"} for j in range(2)]}
            for i, name in enumerate(system_repouts)
        },
        "plots": {
            "enabled": True,
            "families": {"time_series": True, "profiles": False, "interactions": False},
            "distributions": False,
            "autocorrelation": False,
            "bootstrap_tables": False,
            "dpi": 120,
        },
        "statistics": {"enabled": False},
    }
    # Each file below matches an existing *saved-product* consumer. Intervals
    # are invented for visualization only and must never enter real inference.
    summary_root = root / "project_summary"
    statistics = summary_root / "statistics"
    (statistics / "profiles").mkdir(parents=True)
    (statistics / "time_series").mkdir(parents=True)
    (statistics / "solvation").mkdir(parents=True)
    (summary_root / "diffusion").mkdir(parents=True)
    (summary_root / "shared_states").mkdir(parents=True)

    (statistics / "manifest.json").write_text(json.dumps({
        "status": "complete", "methods": {"synthetic_preview": True},
        "families": {family: "enabled" for family in
                     ("time_series", "profiles", "interactions", "solvation_rdfs",
                      "diffusion", "coupling")}
    }))
    profile_rows, condition_rows, contrast_rows, enrichment_rows = [], [], [], []
    population_rows, diagnostics_rows, diffusion_rows, rdf_rows = [], [], [], []
    positions = np.arange(1, 13, dtype=float)
    for system_index, (system, reps) in enumerate(system_repouts.items()):
        for rep_index, repout in enumerate(reps):
            rep = repout.name
            values = 0.24 + 0.025 * system_index + 0.012 * rep_index + (
                0.065 * np.sin(positions / 2.5 + 0.3 * rep_index))
            profile_rows.extend({
                "system": system, "replicate": rep,
                "profile": "lipid_order__POPC__sn1", "position": float(p),
                "value": float(val),
            } for p, val in zip(positions, values))
            enrichment_rows.extend({
                "system": system,
                "metric": f"lipid_footprint_{species}_{leaflet}_0-6A_{suffix}",
                "mean": float(0.65 + 0.35 * system_index + 0.16 * rep_index
                              + 0.12 * species_index + 0.06 * leaflet_index),
            } for suffix in ("fold", "area_density_ratio")
              for species_index, species in enumerate(("POPC", "DOPC"))
              for leaflet_index, leaflet in enumerate(("upper", "lower")))
            diffusion_rows.append({
                "system": system, "metric": "lateral_2d", "replicate": rep,
                "diffusion_nm2_ns": 0.022 + 0.004 * system_index + 0.002 * rep_index,
            })
            frac = np.array([0.24 + 0.05*system_index,
                             0.42 + 0.025*rep_index])
            frac = np.append(frac, 1 - frac.sum())
            population_rows.extend({
                "system": system, "replicate": rep, "state": state,
                "fraction": float(share),
            } for state, share in enumerate(frac))
            diagnostics_rows.append({
                "system": system, "replicate": rep,
                "first_second_half_js_divergence": 0.018 + 0.012*rep_index,
            })

            interaction_dir = repout / "interactions/salt_bridges"
            interaction_dir.mkdir(parents=True)
            pd.DataFrame({
                "residue_1": ["A:ASP42", "A:GLU55", "B:ASP17"],
                "residue_2": ["B:LYS20", "B:ARG21", "A:ARG69"],
                "interaction_class": ["salt_bridge"] * 3,
                "occupancy": [0.40 + .05*system_index, .31, .14+.06*rep_index],
            }).to_parquet(interaction_dir / "occupancy.parquet", index=False)
            pd.DataFrame({
                "lifetime_ps": [280, 600, 1300, 2000, 3500],
                "completed": [True] * 5,
            }).to_parquet(interaction_dir / "lifetimes.parquet", index=False)

            solvation_dir = repout / "extensions/solvation"
            solvation_dir.mkdir(parents=True)
            pd.DataFrame({
                "time_ps": times, "solvent": "water",
                "target_class": "protein", "site_class": "acidic",
                "mean_n_solvent": (4.0 + system_index + 0.3 * rep_index
                                   + 0.4 * np.sin(times / 13500)),
            }).to_parquet(solvation_dir / "species_timeseries.parquet", index=False)
            pd.DataFrame({
                "solvent": "water", "target_id": "protein/acidic/1",
                "duration_ps": [0., 200., 400., 600., 800., 1000.],
                "survival_probability": [1.0, .81, .56, .39, .20, 0.],
            }).to_parquet(solvation_dir / "survival.parquet", index=False)

            contacts = repout / "extensions/lipid_contacts"
            contacts.mkdir(parents=True)
            pd.DataFrame([
                {
                    "time_ps": float(sample_time),
                    "leaflet": "upper", "species_a": a, "species_b": b,
                    "neighbors_per_molecule": 2.2 + .2*system_index + .12*i + .09*j,
                    "enrichment": .88 + .10*system_index + .07*i - .03*j,
                    "pair_count_undirected": 40 + i + j,
                    "molecules_a": 25, "molecules_b": 25,
                }
                for sample_time in (0., 60_000., 120_000.)
                for i, a in enumerate(("POPC", "DOPC"))
                for j, b in enumerate(("POPC", "DOPC"))
            ]).to_parquet(contacts / "matrix.parquet", index=False)

            # Real SDC4 diagnostic renderers consume columns from this
            # *fictional* representative saved table.
            mechanism_dir = repout / "extensions/sdc4_mechanisms"
            mechanism_dir.mkdir(parents=True)
            t = np.linspace(0, 6*np.pi, len(times))
            pd.DataFrame({
                "time_ps": times,
                "sdc4_ectodomain_A_height_A": 18 + 2*np.sin(t + .35*rep_index),
                "sdc4_ectodomain_B_height_A": 20 + 2*np.cos(t + .30*rep_index),
                "sdc4_motif_NxIP_PROA_1_steric_occlusion_fraction": np.clip(
                    .45 + .05*rep_index + .3*np.sin(t+.5+.2*rep_index), 0, 1),
                "sdc4_motif_NxIP_PROB_1_steric_occlusion_fraction": np.clip(
                    .55 - .05*rep_index + .2*np.cos(t+.3+.2*rep_index), 0, 1),
                "sdc4_proline_A_1_phi_deg": -90 + 10*rep_index + 40*np.sin(t),
                "sdc4_proline_A_1_psi_deg": 110 - 9*rep_index + 35*np.cos(t),
            }).to_parquet(mechanism_dir / "timeseries.parquet", index=False)

            grid = np.linspace(-24., 24., 13)
            mesh_x, mesh_y = np.meshgrid(
                (grid[:-1] + grid[1:]) / 2,
                (grid[:-1] + grid[1:]) / 2, indexing="ij"
            )
            layer = np.exp(-(mesh_x**2 + mesh_y**2) / 240)
            density = np.stack([0.13 + .016*rep_index + layer / 8,
                                0.10 + .012*rep_index + layer / 12]).astype(np.float32)
            fold = np.stack([0.85 + .07*rep_index + 0.4*layer,
                             0.95 + .05*rep_index + 0.3*layer]).astype(np.float32)
            np.savez_compressed(
                mechanism_dir / "lipid_density_maps.npz",
                species=np.array(["POPC", "DOPC"]),
                leaflet=np.array(["upper", "lower"]),
                density_per_nm2=density,
                enrichment_fold=fold,
                x_edges_A=grid, y_edges_A=grid,
            )
            pd.DataFrame([
                {
                    "region": region, "species": lipid, "leaflet": "upper",
                    "moiety": moiety, "contact_fraction": 0.12 + 0.05*ridx
                    + 0.07*lidx + 0.025*rep_index + 0.01*system_index,
                    "local_composition_enrichment": 0.8 + 0.1*ridx + .15*lidx + .05*rep_index,
                    "matched_local_contact_fraction": 0.09 + .06*ridx + .025*rep_index,
                    "matched_local_composition_enrichment": 0.7 + .14*ridx + .065*rep_index,
                }
                for ridx, region in enumerate(("TM helix", "Acidic motif", "PRM"))
                for lidx, lipid in enumerate(("POPC", "DOPC"))
                for moiety in ("any", "glycan", "ceramide_head", "ceramide_tail")
            ]).to_parquet(
                mechanism_dir / "lipid_region_specificity.parquet", index=False
            )

        vals = 0.24 + 0.025 * system_index + .006 + 0.065 * np.sin(positions / 2.5)
        condition_rows.extend({
            "system": system, "profile": "lipid_order__POPC__sn1",
            "position": float(p), "position_label": f"C{int(p)}", "mean": float(val),
            "mean_ci_low": float(val - .025),
            "mean_ci_high": float(val + .025),
        } for p, val in zip(positions, vals))
        radii = np.arange(0, 8., .4)
        rdf_rows.extend({
            "system": system, "solvent": "water",
            "target_class": "protein", "site_class": "acidic",
            "r_inner_A": float(r), "r_outer_A": float(r+.4),
            "mean_g_bulk": float(1 + 1.2*np.exp(-((r-2.5)/.72)**2)
                                  + .12*system_index),
            "sd_g_bulk": .12, "n_replicates": 2,
        } for r in radii)

    delta = .025 + .006*np.cos(positions/2.5)
    contrast_rows.extend({
        "profile": "lipid_order__POPC__sn1",
        "system_1": "condition_A", "system_2": "condition_B",
        "position": float(p), "position_label": f"C{int(p)}",
        "difference_system2_minus_system1": float(d),
        "difference_ci_low": float(d-.018),
        "difference_ci_high": float(d+.018),
        "probability_difference_gt_0": .88, "p_value": .13,
    } for p, d in zip(positions, delta))

    pd.DataFrame(profile_rows).to_parquet(
        statistics / "profiles/replicate_profiles.parquet", index=False)
    pd.DataFrame(condition_rows).to_parquet(
        statistics / "profiles/system_profiles.parquet", index=False)
    pd.DataFrame(contrast_rows).to_parquet(
        statistics / "profiles/system_comparisons.parquet", index=False)
    pd.DataFrame(enrichment_rows).to_parquet(
        statistics / "time_series/replicate_statistics.parquet", index=False)
    pd.DataFrame(rdf_rows).to_parquet(
        statistics / "solvation/system_rdf.parquet", index=False)
    pd.DataFrame(diffusion_rows).to_parquet(
        summary_root / "diffusion/replicate_coefficients.parquet", index=False)
    pd.DataFrame(population_rows).to_parquet(
        summary_root / "shared_states/replicate_populations.parquet", index=False)
    pd.DataFrame(diagnostics_rows).to_parquet(
        summary_root / "shared_states/sampling_diagnostics.parquet", index=False)
    # Explicitly invented saved intervals exercise the production forest renderer;
    # these are fixture constants, never estimates from the synthetic traces.
    pd.DataFrame([
        dict(metric='rmsd_A', group=group, system_1='condition_A', system_2='condition_B',
             mean_difference_system2_minus_system1=value, mean_difference_ci_low=low,
             mean_difference_ci_high=high, mean_probability_difference_gt_0=.6,
             p_value=.5, q_value=.9, ci_level=.95)
        for group, value, low, high in ((None, .15, .08, .22), ('chain_A', -.05, -.12, .02), ('chain_B', 0., 0., 0.))
    ]).to_parquet(statistics / 'time_series/system_comparisons.parquet', index=False)
    cfg["statistics"]["enabled"] = True
    cfg["analyses"] = {
        "lipid_contacts": {"enabled": True},
        "sdc4_mechanisms": {"enabled": True},
    }
    cfg["plots"]["families"].update({
        family: True for family in (
            "profiles", "interactions", "solvation", "sdc4", "shared_states",
            "diffusion", "lipid_enrichment", "lipid_contacts"
        )
    })
    return cfg, system_repouts


def mark_synthetic(path: Path):
    """Add a visible caption *outside* the renderer's axes and title."""
    from PIL import ImageFont

    with Image.open(path) as original:
        figure = original.convert("RGB")
    caption_height = 40
    captioned = Image.new("RGB", (figure.width, figure.height + caption_height), "white")
    captioned.paste(figure, (0, caption_height))
    draw = ImageDraw.Draw(captioned)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 17)
    except OSError:
        font = ImageFont.load_default()
    draw.text((12, 9), STAMP, fill="black", font=font)
    draw.line((0, caption_height - 1, figure.width, caption_height - 1),
              fill=(140, 140, 140), width=1)
    captioned.save(path)


def render_preview(root: Path, *, include: list[str] | None = None, presets=None, generic=False):
    root = Path(root)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Refusing to replace preview or results in nonempty {root}")
    root.mkdir(parents=True, exist_ok=True)
    cfg, system_repouts = create_fixture(root)
    cfg["plots"]["table_reports_pdf"] = True
    if generic:
        cfg['analyses']['sdc4_mechanisms']['enabled'] = False
        cfg['plots']['families']['sdc4'] = False
    if presets is not None:
        from mdscope.config import DEFAULTS, merge
        cfg = merge(DEFAULTS, cfg)
        for name in ('rmsd', 'radius_of_gyration', 'membrane', 'membrane_order', 'salt_bridges', 'solvation', 'shared_states'):
            cfg['analyses'][name]['enabled'] = True
        cfg['selections']['groups'] = {name: {'selection': 'protein', 'label': name}
                                       for name in ('chain_A', 'chain_B')}
        for system in cfg['systems'].values():
            system['topology'] = 'SYNTHETIC_NO_TOPOLOGY.pdb'
            for replica in system['replicates']:
                replica['trajectory'] = 'SYNTHETIC_NO_TRAJECTORY.xtc'
        if generic:
            cfg['analyses']['sdc4_mechanisms']['enabled'] = False
            cfg['plots']['families']['sdc4'] = False
        cfg['plots']['presets'] = list(presets)
        from mdscope.planning import resolve_figure_presets
        cfg, omissions = resolve_figure_presets(cfg, root, root / 'synthetic.yaml')
        (root / 'preset_omissions.json').write_text(json.dumps(omissions, indent=2))
    if include is not None:
        cfg["plots"]["include"] = list(include)
    # These are the same plotting functions used by runner.finalize_project.
    from mdscope.reporting.membrane_figures import (
        plot_lipid_contacts,
        plot_lipid_enrichment_matrix,
        plot_project_diffusion,
    )
    from mdscope.reporting.sdc4_figures import plot_sdc4_mechanisms
    from mdscope.reporting.solvation_figures import plot_solvation

    plot_all(cfg, root, system_repouts)
    plot_lipid_enrichment_matrix(cfg, root)
    plot_sdc4_mechanisms(cfg, root, system_repouts)
    plot_lipid_contacts(cfg, root, system_repouts)
    plot_project_diffusion(cfg, root)
    plot_solvation(cfg, root, system_repouts)
    figures = sorted(
        path for path in root.rglob("*.png")
        if ("project_summary" in path.parts or
            ("extensions" in path.parts and "sdc4_mechanisms" in path.parts))
    )
    minimum = 1 if include is not None or presets is not None else 15
    if len(figures) < minimum:
        raise RuntimeError(f"Expected at least {minimum} figure preview(s); got {len(figures)}")
    for path in figures:
        mark_synthetic(path)
    # Export only schema and a bounded first-row sample, never full MD tables.
    from runpy import run_path
    export = run_path(str(Path(__file__).with_name("preview_tables.py")))[
        "export_preview_tables"
    ]
    tables = export(root, sample_rows=5)
    import mdscope
    record = {
        "notice": STAMP,
        "renderer_module": str(Path(mdscope.__file__).resolve()),
        "input_data_kind": "deterministic synthetic Parquet",
        "requested_plot_ids": include,
        "plots": [str(path.relative_to(root)) for path in figures],
        "reports": [path.relative_to(root).as_posix() for path in sorted(root.rglob("statistics_tables.pdf"))],
        "table_previews": tables,
        "table_inventory_csv": "table_previews/table_inventory.csv",
        "plot_ids": {},
    }
    try:
        from mdscope.reporting.plot_ids import id_from_figure_path
    except ImportError:
        id_from_figure_path = None
    if id_from_figure_path is not None:
        record["plot_ids"] = {
            str(path.relative_to(root)): identifier
            for path in figures
            if (identifier := id_from_figure_path(path)) is not None
        }
    (root / "preview_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--include", action="append", help="Select exact plot ID (repeatable)")
    parser.add_argument("--preset", action="append", choices=("primary", "secondary", "tertiary"))
    parser.add_argument("--generic", action="store_true", help="Exclude SDC4 biological presentation")
    args = parser.parse_args()
    record = render_preview(args.output, include=args.include, presets=args.preset, generic=args.generic)
    print(f"Rendered {len(record['plots'])} marked synthetic previews in {args.output}")


if __name__ == "__main__":
    main()
