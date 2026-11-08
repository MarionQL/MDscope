"""PDF presentation of persisted evidence; ReportLab is imported only on export."""
from __future__ import annotations
from datetime import datetime, timezone
from html import escape
from io import BytesIO
import json
from pathlib import Path
import numpy as np
import pandas as pd
from ..config import presentation_enabled, system_label
from .report_inventory import inference_status, _read, _summary_inventory, _metric_groups, read_report_table
from .figure_inventory import figure_inventory, report_coverage, metric_coverage, current_statistics, saved_statistics_methods
from .solvation_figures import (solvation_inventory, solvation_rdf_figure,
    solvation_timeseries_figure, solvation_survival_figure,
    solvation_stat_rows, solvation_comparison_rows)

KEY_METRICS = (
    "sdc4_ectodomain_A_height_A", "sdc4_ectodomain_B_height_A",
    "sdc4_ectodomain_A_tilt_deg", "sdc4_ectodomain_B_tilt_deg",
    "sdc4_ectodomain_A_membrane_contact_state",
    "sdc4_ectodomain_B_membrane_contact_state",
    "sdc4_tm_crossing_angle_deg", "sdc4_tm_interface_residue_pairs",
    "sdc4_gxxxg_A_to_other_tm_residue_pairs",
    "sdc4_gxxxg_B_to_other_tm_residue_pairs",
    "sdc4_gxxxg_interface_residue_pairs",
    "sdc4_local_thickness_A", "sdc4_local_upper_curvature_Ainv",
    "sdc4_local_lower_curvature_Ainv", "sdc4_local_minus_SCH",
)


def _fmt(value, digits=3):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "Unavailable"
    if not np.isfinite(x):
        return "Unavailable"
    if abs(x) >= 10000 or (0 < abs(x) < .001):
        return f"{x:.2g}"
    return f"{x:.{digits}f}"


def _safe(value):
    return escape(str(value), quote=False).replace("\n", " ")


def _maps_figure(system_repouts, species, field, dpi=135):
    """Replicate-mean density/enrichment; never average absent species as zero."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    systems = list(system_repouts)
    fig, axes = plt.subplots(
        len(systems), 2, figsize=(9.1, 3.15*len(systems)),
        squeeze=False, constrained_layout=True,
    )
    found = False
    for i, (system, paths) in enumerate(system_repouts.items()):
        for j, leaflet in enumerate(("upper", "lower")):
            maps = []
            extent = None
            for repout in paths:
                path = repout / "extensions/sdc4_mechanisms/lipid_density_maps.npz"
                if not path.exists():
                    continue
                with np.load(path, allow_pickle=False) as z:
                    indices = np.flatnonzero(
                        (z["species"].astype(str) == species)
                        & (z["leaflet"].astype(str) == leaflet)
                    )
                    if len(indices) != 1:
                        continue
                    arr = np.asarray(z[field][indices[0]], float)
                    if not np.any(np.isfinite(arr)):
                        continue
                    maps.append(arr)
                    extent = (float(z["x_edges_A"][0]), float(z["x_edges_A"][-1]),
                              float(z["y_edges_A"][0]), float(z["y_edges_A"][-1]))
            ax = axes[i, j]
            if maps:
                # A species may be absent in a condition. Never replace
                # absence with a zero map, and don't combine incompatible grids.
                if any(a.shape != maps[0].shape for a in maps):
                    raise ValueError(f"Incompatible lipid map grids for {species}")
                img = ax.imshow(
                    np.nanmean(np.stack(maps), axis=0).T, origin="lower",
                    extent=extent, interpolation="nearest", aspect="equal",
                    cmap="viridis", vmin=0,
                )
                fig.colorbar(img, ax=ax, fraction=.045, pad=.03)
                found = True
                title = f"{system} / {leaflet} (n={len(maps)} reps)"
            else:
                ax.text(.5, .5, "Not observed / unavailable",
                        transform=ax.transAxes, ha="center", va="center")
                title = f"{system} / {leaflet}"
            ax.set_title(title, fontsize=9)
            ax.set_xlabel("TM-centered x (Å)", fontsize=8)
            ax.set_ylabel("TM-centered y (Å)", fontsize=8)
    fig.suptitle(f"{species}: {field.replace('_', ' ')}", fontsize=13)
    if not found:
        plt.close(fig)
        return None
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi)
    plt.close(fig)
    buffer.seek(0)
    return buffer


def _build_pdf(path, cfg, root, system_repouts, options):
    try:
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_CENTER
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
            PageBreak, Image, KeepTogether, HRFlowable,
        )
    except ImportError as exc:
        raise RuntimeError(
            "PDF reporting requires reportlab. Install with: pip install reportlab"
        ) from exc

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="ReportTitle", parent=styles["Title"], fontName="Helvetica-Bold",
        fontSize=21, leading=26, textColor=colors.HexColor("#17324D"),
        spaceAfter=12,
    ))
    styles.add(ParagraphStyle(
        name="ReportHeading", parent=styles["Heading1"],
        fontSize=13, leading=17, textColor=colors.HexColor("#17324D"),
        spaceBefore=12, spaceAfter=6,
    ))
    styles.add(ParagraphStyle(
        name="ReportSub", parent=styles["Heading2"],
        fontSize=10.2, leading=13, spaceBefore=9, spaceAfter=4,
    ))
    styles.add(ParagraphStyle(
        name="ReportBody", parent=styles["BodyText"],
        fontSize=8.6, leading=12, spaceAfter=6,
    ))
    styles.add(ParagraphStyle(
        name="ReportSmall", parent=styles["BodyText"],
        fontSize=7.1, leading=9.4, spaceAfter=3,
    ))
    styles.add(ParagraphStyle(
        name="ReportTable", parent=styles["BodyText"],
        fontSize=8, leading=10.4, splitLongWords=False,
    ))
    from ..config import statistics_family_config
    methods = saved_statistics_methods(root)
    ci_label = "CI (level unrecorded)"
    if methods is not None:
        saved_cfg = {**cfg, 'statistics': methods}
        scalar_cfg = statistics_family_config(saved_cfg, 'time_series') or saved_cfg
        ci_label = f"{100*float(scalar_cfg['statistics'].get('ci', .95)):g}% CI"
    doc = SimpleDocTemplate(
        str(path), pagesize=letter, rightMargin=43, leftMargin=43,
        topMargin=48, bottomMargin=45,
        title=f"MDscope evidence report - {cfg['project']['name']}",
        author="MDscope",
    )
    story = []
    displayed_figures = set()
    table_sources = set()
    verified_figures = set(figure_inventory(root,system_repouts,cfg).figure)
    active_section = True
    current_family = 'overview'
    section_families = {
        'Systems and analyzed sampling':'overview',
        'Condition-level scalar results':'time_series',
        'Replicate-level consistency: primary mechanistic metrics':'time_series',
        'Shared conformational-state model':'shared_states',
        'Matched local lipid specificity':'sdc4',
        'Solvation, hydration and ion coordination':'solvation',
        'Selection audits and methodological limitations':'methods',
        'Protein-centered lipid enrichment':'lipid_enrichment',
        'Scientific figure atlas and coverage audit':'atlas',
        'Files to share for deeper analysis':'methods',
    }
    def read_table(path):
        data = read_report_table(path, cfg, root)
        if not data.empty:
            table_sources.add(str(Path(path).relative_to(root)))
        return data
    def display_cell(value):
        if isinstance(value,str):
            names = {'target_id':'Target', 'mean_n_solvent':'Mean solvent count',
                     'hydrated_fraction':'Hydrated fraction', 'right_censored':'Right censored',
                     'left_censored':'Left censored', 'duration_ps':'Reported duration (ps)',
                     'survival_probability':'Survival probability', 'g_bulk':'g(r), cell-density reference'}
            if value in names:
                return names[value]
            from .labels import _metric_meta, METRIC_LABELS
            if value.startswith('sdc4_') or value in METRIC_LABELS:
                label, units = _metric_meta(value)
                return label + (f' ({units})' if units else '')
            if value in cfg.get('systems',{}):
                return system_label(cfg,value)
        return value
    def para(text, style="ReportBody"):
        return Paragraph(_safe(text), styles[style])
    def header(text):
        nonlocal active_section, current_family
        current_family = section_families.get(text, 'lipid_enrichment' if text.endswith(' maps') else current_family)
        active_section = presentation_enabled(cfg,'report',current_family)
        if not active_section:
            return
        story.append(Paragraph(_safe(text), styles["ReportHeading"]))
    def sub(text):
        text = display_cell(text)
        if not active_section:
            return
        story.append(Paragraph(_safe(text), styles["ReportSub"]))
    def body(text):
        if not active_section:
            return
        story.append(para(text))
    def table(rows, widths=None):
        if not rows or not active_section:
            return
        rendered = [
            [Paragraph("<br/>".join(_safe(part) for part in str(display_cell(cell)).split("\n")), styles["ReportTable"]) for cell in row]
            for row in rows
        ]
        obj = Table(rendered, colWidths=widths, repeatRows=1, hAlign="LEFT")
        obj.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E9EFF4")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#17324D")),
            ("LINEBELOW", (0, 0), (-1, 0), .6, colors.HexColor("#A6BAC9")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1),
             [colors.white, colors.HexColor("#F7F9FB")]),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ]))
        story.append(obj)
        story.append(para(f"Source family: {current_family}. "
                          "Sampling and source-file provenance are recorded in report_provenance.json; raw identifiers are retained in the bundle.", "ReportSmall"))
        story.append(Spacer(1, 8))
    def figure(path_or_buffer, caption, max_width=7.05*inch, max_height=7.6*inch):
        if path_or_buffer is None or not active_section:
            return
        if isinstance(path_or_buffer, (str, Path)) and Path(path_or_buffer).exists():
            relative = str(Path(path_or_buffer).relative_to(root))
            if relative not in verified_figures or relative in displayed_figures:
                return
            displayed_figures.add(relative)
        caption += f" Source family: {current_family}. Sampling and source paths are in the companion manifest."
        try:
            from PIL import Image as PILImage
            if isinstance(path_or_buffer, (str, Path)):
                with PILImage.open(path_or_buffer) as im:
                    w, h = im.size
            else:
                path_or_buffer.seek(0)
                with PILImage.open(path_or_buffer) as im:
                    w, h = im.size
                path_or_buffer.seek(0)
            ratio = min(max_width/w, max_height/h)
            img = Image(path_or_buffer, width=w*ratio, height=h*ratio)
            img.hAlign = "CENTER"
            story.append(img)
            story.append(para(caption, "ReportSmall"))
            story.append(Spacer(1, 7))
        except Exception as exc:
            body(f"Figure unavailable ({caption}): {exc}")
    def footer(canvas, document):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#CBD6E0"))
        canvas.line(43, 38, letter[0]-43, 38)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor("#52687B"))
        canvas.drawString(43, 27, "MDscope | Automated descriptive evidence report")
        canvas.drawRightString(letter[0]-43, 27, f"Page {document.page}")
        canvas.restoreState()

    title = cfg["project"].get("name", "MDscope project")
    story.append(Paragraph("Molecular dynamics evidence report",
                           styles["ReportTitle"]))
    body(f"Project: {title}")
    body(f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    body("Purpose: a compact, self-contained PDF for reviewing condition-dependent "
         "MD observables and discussing mechanistic hypotheses. "
         "All numbers are read from completed analysis outputs.")
    meta_file = root / "run_metadata.json"
    if meta_file.exists():
        meta = json.loads(meta_file.read_text())
        body(f"MDscope version recorded in run metadata: {meta.get('mdscope_version', '?')}")
    header("Systems and analyzed sampling")
    inventory = _summary_inventory(root, cfg, system_repouts)
    if inventory.empty:
        body("No replicate output tables found.")
    else:
        rows = [["System / replicate", "Frames", "First / last time (ns)",
                 "NxIP", "HS", "Maps"]]
        for item in inventory.itertuples():
            rows.append([
                f"{item.system} / {item.replicate}", item.frames,
                f"{_fmt(item.first_ps/1000)} / {_fmt(item.last_ps/1000)}",
                item.nxip_motifs, item.hs_chains,
                "yes" if item.lipid_maps else "no",
            ])
        table(rows, [175, 50, 122, 40, 36, 42])
    if not inventory.empty and any(inventory[col].nunique(dropna=False)>1 for col in ('frames','first_ps','last_ps')):
        body("UNEQUAL SAMPLING: replicas have different analyzed frame counts or time windows. System estimates weight replicas equally; these are not matched-window comparisons. Inspect the sampling inventory before interpreting contrasts.")
    tr = cfg.get("trajectory", {})
    body(f"Trajectory frame slice: start={tr.get('start')}, "
         f"stop={tr.get('stop')}, stride={tr.get('stride')}; "
         "reported time windows are those actually sampled.")
    if not inventory.empty and (inventory["hs_chains"] == 0).any():
        body("IMPORTANT: At least one replicate has no mapped HS chains. "
             "HS-specific observables for that replicate are absent, not zero. "
             "Confirm topology covalent attachments before interpreting HS effects.")
    body("Interpretation guardrails: frames are temporally correlated; "
         "independent trajectories, not frames, are the unit of condition-level replication. "
         "Permutation p-values with only a few replicates have limited resolution. "
         "Report sections are descriptive and do not establish causality or equilibrium.")
    story.append(PageBreak())

    stats_dir = root / "project_summary/statistics/time_series"
    systems = read_table(stats_dir / "system_statistics.parquet")
    comparisons = read_table(stats_dir / "system_comparisons.parquet")
    reps = read_table(stats_dir / "replicate_statistics.parquet")
    max_rows = int(options.get("max_metrics_per_section", 45))
    header("Condition-level scalar results")
    from .figure_inventory import statistics_status
    status = statistics_status(cfg,root,'time_series')
    if status == 'legacy_unverified':
        body("Provenance status: legacy/unverified statistical inputs. Original input/config identity cannot be verified; regenerate statistics to establish current provenance.")
    elif status not in {'valid', 'disabled'}:
        body(f"Statistical inputs excluded: {status}. Complete statistics with the active configuration before interpreting condition-level results.")
    body("System mean and bootstrap confidence intervals are taken directly from "
         "MDscope's hierarchical statistics. Contrasts are system 2 minus system 1; "
         "q-values use the configured multiple-testing correction. "
         "Unavailable denotes missing or non-finite data; a measured zero remains zero.")
    if systems.empty:
        body("System statistics not available. Enable statistics and complete finalization.")
    else:
        system_ids = list(system_repouts)
        if len(system_ids) != 2:
            body("The compact contrast table shows the first two configured systems; "
                 "all pairwise comparisons are included in the ZIP.")
        first, second = (system_ids + ["not_configured", "not_configured"])[:2]
        body(f"Contrast direction: {system_label(cfg,second)} minus {system_label(cfg,first)}.")
        for category, metrics in _metric_groups(systems):
            sub(category)
            total = len(metrics)
            metrics = metrics[:max_rows]
            if total > len(metrics):
                body(f"Showing {len(metrics)} of {total} metrics; the companion "
                     "review ZIP contains all full-precision statistics.")
            rows = [[
                "Metric", "Units", system_label(cfg,first) + "\nMean ["+ci_label+"]", system_label(cfg,second) + "\nMean ["+ci_label+"]",
                "Second minus first\n["+ci_label+"]; p / q",
            ]]
            for metric in metrics:
                slice_ = systems[systems["metric"] == metric]
                groups = (slice_["group"].fillna("").astype(str).unique()
                          if "group" in slice_ else [""])
                for group in groups:
                    part = (slice_[slice_["group"].fillna("").astype(str) == group]
                            if "group" in slice_ else slice_)
                    from .labels import _metric_meta
                    label, units = _metric_meta(metric)
                    if units == metric:
                        units = "See schema"
                    if group:
                        label += " [" + group.replace("_", " ").replace("/", " / ") + "]"
                    values = {}
                    for item in part.itertuples():
                        values[item.system] = (
                            f"{_fmt(item.mean)}\n[{_fmt(item.mean_ci_low)}, "
                            f"{_fmt(item.mean_ci_high)}]"
                        )
                    delta = "Unavailable"
                    if not comparisons.empty:
                        contrast = comparisons[
                            (comparisons["metric"] == metric)
                            & (comparisons["system_1"] == first)
                            & (comparisons["system_2"] == second)
                        ]
                        if "group" in contrast:
                            contrast = contrast[
                                contrast["group"].fillna("").astype(str) == group
                            ]
                        if len(contrast):
                            c = contrast.iloc[0]
                            delta = (
                                f"{_fmt(c['mean_difference_system2_minus_system1'])} "
                                f"\n[{_fmt(c['mean_difference_ci_low'])}, "
                                f"{_fmt(c['mean_difference_ci_high'])}]\n"
                                f"p={_fmt(c.get('mean_p_value'))} / "
                                f"q={_fmt(c.get('mean_q_value'))}\n"
                                f"{inference_status(c.get('mean_q_value'),c.get('mean_p_value'))}"
                            )
                    rows.append([
                        label, units.replace("degrees","deg") or "-", values.get(first, "Unavailable"), values.get(second, "Unavailable"),
                        delta,
                    ])
            table(rows, [145, 45, 100, 100, 120])
    story.append(PageBreak())

    header("Replicate-level consistency: primary mechanistic metrics")
    if reps.empty:
        body("Replicate-level statistics unavailable.")
    else:
        candidates = list(KEY_METRICS)
        hs_metrics = sorted(
            m for m in reps["metric"].dropna().astype(str).unique()
            if m.startswith("sdc4_hs_") and (
                m.endswith("_orientation_cos")
                or m.endswith("_lipid_contact_lipids")
                or m.endswith("_hydration_water_per_heavy_atom")
            )
        )
        candidates.extend(hs_metrics[:int(options.get("max_hs_key_metrics", 20))])
        for metric in candidates:
            subset = reps[reps["metric"] == metric]
            if subset.empty:
                continue
            sub(metric)
            rows = [["System", "Replicate", "Mean", "Block-bootstrap CI"]]
            for item in subset.itertuples():
                rows.append([
                    item.system, item.replicate, _fmt(item.mean),
                    f"{_fmt(item.mean_ci_low)} to {_fmt(item.mean_ci_high)}",
                ])
            table(rows, [145, 82, 72, 166])

    state_dir = root / "project_summary/shared_states"
    model_path = state_dir / "model.json"
    if model_path.exists() and presentation_enabled(cfg,'report','shared_states'):
        story.append(PageBreak())
        header("Shared conformational-state model")
        model = json.loads(model_path.read_text())
        body(f"States: {model.get('n_states')}; training replicates: "
             f"{model.get('n_replicates')}; "
             f"training frames per replicate: "
             f"{model.get('training_frames_per_replicate')}.")
        body("Features: " + ", ".join(model.get("features", [])))
        if model.get("missing_frame_policy"):
            body(model["missing_frame_policy"])
        body("One common scaler and clustering model are used across conditions. "
             "State labels are comparable between systems, but their populations "
             "are not direct free-energy estimates.")
        centroids = read_table(state_dir / "state_centroids.parquet")
        if not centroids.empty:
            sub("State centroids (original feature units)")
            for row in centroids.itertuples(index=False):
                state = getattr(row, "state")
                rows = [["Feature", "Centroid"]]
                for feature in model.get("features", []):
                    rows.append([feature, _fmt(getattr(row, feature))])
                story.append(para(f"State {state}", "ReportSmall"))
                table(rows, [330, 135])
        population = read_table(state_dir / "system_populations.parquet")
        if not population.empty:
            sub("Replicate-averaged state occupancy")
            rows = [["System", "State", "Mean fraction", "Replicate CI", "n reps"]]
            for item in population.itertuples():
                rows.append([
                    item.system, item.state, _fmt(item.mean_replicate_fraction),
                    f"{_fmt(item.replicate_bootstrap_ci_low)} to "
                    f"{_fmt(item.replicate_bootstrap_ci_high)}",
                    item.n_replicates,
                ])
            table(rows, [120, 40, 80, 155, 65])
        diagnostics = read_table(state_dir / "sampling_diagnostics.parquet")
        if not diagnostics.empty:
            sub("Sampling diagnostics")
            rows = [["System", "Replicate", "Half-vs-half JS", "Transitions", "Frames"]]
            for item in diagnostics.itertuples():
                rows.append([
                    item.system, item.replicate,
                    _fmt(item.first_second_half_js_divergence),
                    item.observed_transitions, item.frames,
                ])
            table(rows, [130, 82, 95, 80, 65])
            if "missing_feature_frames" in diagnostics:
                body("Missing state-feature frames by replica: " + "; ".join(
                    f"{item.system}/{item.replicate}: {item.missing_feature_frames} "
                    f"of {item.input_frames}; {item.temporal_resampling_status}"
                    for item in diagnostics.itertuples()))
            body("First/second-half Jensen-Shannon divergence is a descriptive "
                 "sampling diagnostic; a low value is not proof of convergence.")
        figure(state_dir / "state_populations.png",
               "Shared state occupancies; black points denote independent replicates.")
        figure(state_dir / "state_sampling_diagnostics.png",
               "First/second-half state-population divergence.")
        conditional = read_table(state_dir / "state_conditioned_hs_solvation.parquet")
        if not conditional.empty:
            sub("State- and contact-conditioned HS solvation")
            max_cond = int(options.get("max_conditioned_rows", 80))
            rows = [["System / rep", "HS", "State", "Contact", "Frames", "Waters/atom"]]
            for item in conditional.head(max_cond).itertuples():
                rows.append([
                    f"{item.system}/{item.replicate}", item.hs_chain,
                    item.state, item.membrane_contact_state, item.frames,
                    _fmt(item.hydration_water_per_heavy_atom_mean),
                ])
            table(rows, [130, 60, 42, 49, 49, 115])
            if len(conditional) > max_cond:
                body(f"{len(conditional)-max_cond} additional conditioned rows "
                     "are included in the companion ZIP.")

    lipid_stats = read_table(state_dir / "lipid_specificity_system_statistics.parquet")
    if not lipid_stats.empty and presentation_enabled(cfg,'report','sdc4'):
        matched = lipid_stats.loc[lipid_stats.measure.isin((
            "mean_matched_local_contact_fraction", "mean_matched_local_enrichment"))]
        if not matched.empty:
            story.append(PageBreak())
            header("Matched local lipid specificity")
            body("Each system mean gives equal weight to independent replicas. "
                 "Intervals resample replicas and are descriptive with small n. "
                 "The fraction and enrichment condition on local anchor membership; "
                 "they are not affinities. Full tables and per-replica valid-frame "
                 "counts are included in the evidence bundle.")
            limit = int(options.get("max_conditioned_rows", 80))
            for metric, label in (
                ("mean_matched_local_contact_fraction", "Local contact fraction"),
                ("mean_matched_local_enrichment", "Local composition enrichment"),
            ):
                panel = matched.loc[matched.measure.eq(metric)]
                if panel.empty:
                    continue
                sub(label)
                rows = [["System", "Region", "Species / leaflet", "Mean", "Replicate CI", "n reps"]]
                for item in panel.head(limit).itertuples():
                    rows.append([item.system, item.region,
                                 f"{item.species} / {item.leaflet}",
                                 _fmt(item.replicate_mean),
                                 f"{_fmt(item.replicate_bootstrap_ci_low)} to "
                                 f"{_fmt(item.replicate_bootstrap_ci_high)}",
                                 item.n_replicates])
                table(rows, [90, 105, 105, 50, 100, 60])
                if len(panel) > limit:
                    body(f"{len(panel)-limit} additional rows are in the evidence bundle.")

    # Solvation results are reported explicitly, rather than hidden among
    # generic scalar metrics. Large per-target/event tables are indexed here
    # and included at full precision in the companion ZIP.
    sol_inventory = solvation_inventory(system_repouts) if presentation_enabled(cfg,'report','solvation') else pd.DataFrame()
    if not sol_inventory.empty:
        story.append(PageBreak())
        header("Solvation, hydration and ion coordination")
        body("Solvation counts are unique solvent entities within configured shells. "
             "The hydration fraction is the fraction of targets with at least one "
             "solvent in the shell. RDF g(r) uses selected solvent count / full periodic cell volume "
             "without excluded-volume correction. Replicates, not frames or "
             "molecular sites, are independent sampling units.")
        summary = solvation_stat_rows(root) if current_statistics(cfg,root,'time_series') else pd.DataFrame()
        if not summary.empty:
            sub("System-level hydration and coordination statistics")
            body("Moving-block and hierarchical bootstrap confidence intervals "
                 "reuse the project-wide statistics settings.")
            rows = [["System", "Solvent / target class / site", "Metric",
                     "Mean ["+ci_label+"]", "n reps"]]
            for item in summary.head(55).itertuples():
                rows.append([item.system, str(item.group).replace("solvation/",""),
                    item.metric, f"{_fmt(item.mean)}\n[{_fmt(item.mean_ci_low)}, "
                    f"{_fmt(item.mean_ci_high)}]", item.n_replicates])
            table(rows,[77,174,91,125,43])
            if len(summary)>55:
                body(f"Showing 55 of {len(summary)} system/metric rows. "
                     "All rows are preserved in the review ZIP.")
        comparisons_sol=solvation_comparison_rows(root) if current_statistics(cfg,root,'time_series') else pd.DataFrame()
        if not comparisons_sol.empty:
            sub("Between-condition solvation differences")
            rows=[["Comparison","Solvent / target class / site","Metric",
                   "Difference ["+ci_label+"]"]]
            for item in comparisons_sol.head(35).itertuples():
                rows.append([f"{item.system_2} - {item.system_1}",
                    str(item.group).replace("solvation/",""),item.metric,
                    f"{_fmt(item.mean_difference_system2_minus_system1)} "
                    f"[{_fmt(item.mean_difference_ci_low)}, "
                    f"{_fmt(item.mean_difference_ci_high)}]"])
            table(rows,[100,166,92,152])
            if len(comparisons_sol)>35:
                body(f"Showing 35 of {len(comparisons_sol)} comparisons; "
                     "full results including p/q values are in the ZIP.")
        rdf_figure=solvation_rdf_figure(root,system_repouts) if current_statistics(cfg,root,'solvation_rdfs') and presentation_enabled(cfg,'plots','solvation') else None
        if rdf_figure:
            sub("Radial distribution functions")
            figure(rdf_figure[0],"Replicate-mean RDFs with ±1 replicate SD bands. "
                   f"Showing {rdf_figure[1]} of {rdf_figure[2]} chemical classes; "
                   "all numerical curves are in the review ZIP.")
        series_figure=solvation_timeseries_figure(system_repouts) if presentation_enabled(cfg,'plots','solvation') else None
        if series_figure:
            sub("Coordination over simulation time")
            figure(series_figure[0],"Per-replicate coordination time series. "
                   f"Showing {series_figure[1]} of {series_figure[2]} chemical classes.")
        survival_figure=solvation_survival_figure(system_repouts) if presentation_enabled(cfg,'plots','solvation') else None
        if survival_figure:
            sub("Observed solvent-contact survival")
            figure(survival_figure[0],"Descriptive survival curves from sampled "
                   "contact episodes, with recorded censoring. "
                   f"Showing up to {survival_figure[1]} curves; no uncertainty "
                   "bands are inferred.")
        # Compact previews of the remaining scientific observables. The
        # full data remain lossless in the bundle; no invented CI is shown.
        for filename,title_,columns in (
            ("target_summary","Per-target hydration and ion coordination",
             ["solvent","target_id","mean_n_solvent","hydrated_fraction"]),
            ("shell_candidates","Candidate hydration-shell boundaries",
             ["solvent","target_class","site_class","suggested_cutoff_A","status"]),
            ("shell_stability","Hydration-shell convergence diagnostics",
             ["solvent","target_class","site_class","median_cutoff_A","status"]),
            ("survival","Water/ion contact survival",
             ["solvent","target_id","time_ps","survival_probability"]),
            ("hydration_sasa","Hydration and solvent accessibility",
             ["time_ps","target_id","n_solvent","sasa_nm2","sasa_coverage"]),
            ("water_density_3d","Local 3D water-density voxel values",
             ["target_id","x_A","y_A","z_A","g_bulk"]),
            ("bridges","Solvent-mediated target bridges",
             ["time_ps","solvent","target_a","target_b"]),
            ("water_hbond_bridges","Geometry-qualified water bridges",
             ["time_ps","target_a","target_b"]),
            ("episodes","Sampled solvent residence episodes",
             ["solvent","target_id","duration_ps","right_censored"]),
        ):
            examples=[]
            total=0
            for system_,paths_ in system_repouts.items():
                for rep_ in paths_:
                    path_=rep_/"extensions"/"solvation"/(filename+".parquet")
                    if not path_.exists():
                        continue
                    import pyarrow.parquet as pq
                    pf=pq.ParquetFile(path_)
                    total+=pf.metadata.num_rows
                    if len(examples)>=8:
                        continue
                    available=[c for c in columns if c in pf.schema_arrow.names]
                    if not available:
                        continue
                    batch=next(pf.iter_batches(batch_size=8,columns=available),None)
                    if batch is None:
                        continue
                    for record in batch.to_pylist()[:max(0,8-len(examples))]:
                        examples.append((f"{system_}/{rep_.name}",record))
            if not total:
                continue
            sub(title_)
            body(f"{total} rows available across replicas. "
                 "The following are representative table entries, not "
                 "statistical summaries; complete results are in the ZIP.")
            keys=list(examples[0][1]) if examples else []
            if keys:
                rows=[["System / replicate"]+keys]
                for origin,record in examples:
                    rows.append([origin]+[
                        _fmt(record[k]) if isinstance(record[k],(int,float,np.number))
                        else str(record[k])[:38] for k in keys])
                width=510/(len(keys)+1)
                table(rows,[width]*(len(keys)+1))
        sub("Solvation output inventory")
        body("All generated solvation datasets are indexed below. Large "
             "per-frame, per-target, residence-event and 3D voxel tables are "
             "not expanded into thousands of PDF pages. The companion ZIP "
             "contains all full-precision solvation Parquet files, plus bounded "
             "CSV previews for convenient inspection.")
        rows=[["System / replicate","Dataset","Rows"]]
        for item in sol_inventory.itertuples():
            rows.append([f"{item.system}/{item.replicate}",item.table,item.rows])
        # Split long inventories to keep ReportLab tables manageable.
        for start in range(1,len(rows),45):
            table([rows[0]]+rows[start:start+45],[175,230,105])
        body("Residence survival curves are descriptive and account for "
             "recorded censoring; no confidence interval is implied. "
             "3D density grids are sparse, target-local and only sampled "
             "within the documented spherical domain. Missing or disabled "
             "outputs are not treated as zero.")

    story.append(PageBreak())
    header("Selection audits and methodological limitations")
    for system, paths in system_repouts.items():
        for repout in paths:
            path_ = repout / "extensions/sdc4_mechanisms/selection_audit.json"
            if not path_.exists():
                continue
            audit = json.loads(path_.read_text())
            sub(f"{system} / {repout.name}")
            body("NxIP motif topology residue indices: "
                 + json.dumps(audit.get("motif_resindices", {}), sort_keys=True))
            body("GXXXG motif topology residue indices: "
                 + json.dumps(audit.get("gxxxg_resindices", {}), sort_keys=True))
            body("Mapped HS chain atom counts: "
                 + json.dumps(audit.get("hs_atom_counts", {}), sort_keys=True))
            body("Lipid anchor counts: "
                 + json.dumps(audit.get("lipid_counts", {}), sort_keys=True))
    body("NxIP steric occlusion is a contact proxy, not a binding free energy. "
         "Motif SASA depends on the selection used in the underlying SASA analysis. "
         "HS hydration and ion coordination use distance cutoffs, not electrostatic "
         "free energies. Membrane curvature is a small-slope quadratic estimate; "
         "lipid density maps are not corrected for protein-excluded area.")
    body("Local lipid specificity uses matched anchor-in-disk contact and availability "
         "populations. Legacy unmatched local enrichment is retained only for "
         "comparison. Species absent from a local population have undefined "
         "contact fractions; frames with no local contacts have undefined "
         "composition enrichment. Means exclude these frames; inspect the "
         "reported valid-frame counts. Enrichment is not binding affinity.")
    body("Residence durations are sampled-contact estimators. Lipid duration_ps "
         "includes the next interval (terminal median at the right edge); "
         "observed_span_ps excludes imputation. Solvation duration_ps retains "
         "the observed span and sampled_duration_ps records the imputed value. "
         "Censoring and trajectory stride limit kinetic interpretation.")
    body("Do not infer a membrane-recruitment mechanism from a difference in "
         "one observable alone. Compare replicate-level effects, state occupancy, "
         "orientation and contacts, and temporal sampling together.")

    if presentation_enabled(cfg,'plots','lipid_enrichment') and presentation_enabled(cfg,'report','lipid_enrichment'):
        fig = root / "project_summary/figures/membrane/lipid_enrichment_heatmap.png"
        if fig.exists():
            story.append(PageBreak())
            header("Protein-centered lipid enrichment")
            figure(fig, "Annular lipid enrichment/depletion, calculated using "
                   "leaflet-specific composition baselines.")
        if bool(options.get("include_lipid_maps", True)):
            available = set()
            for paths in system_repouts.values():
                for repout in paths:
                    p = repout / "extensions/sdc4_mechanisms/lipid_density_maps.npz"
                    if p.exists():
                        with np.load(p, allow_pickle=False) as z:
                            available.update(z["species"].astype(str).tolist())
            requested = options.get("map_species") or [
                "GM1a", "BGAL", "CHL1", "CER160"
            ]
            for species in requested:
                if species not in available:
                    continue
                for field in ("density_per_nm2", "enrichment_fold"):
                    try:
                        buf = _maps_figure(system_repouts, species, field)
                    except Exception as exc:
                        body(f"Unable to render {species} {field}: {exc}")
                        continue
                    if buf is not None:
                        story.append(PageBreak())
                        header(f"{species} - {'density' if field.startswith('density') else 'enrichment'} maps")
                        figure(buf,
                               "Replicate-mean, TM-centered leaflet maps. "
                               "A missing species is not plotted as zero enrichment.")

    # Standardized, analysis-agnostic figure gallery. Every primary figure
    # generated by the plotting pipeline is eligible for the report. Keep a
    # configurable bound on PDF pages; the ZIP always carries all figures.
    story.append(PageBreak())
    header("Scientific figure atlas and coverage audit")
    coverage=report_coverage(root,system_repouts,cfg)
    if not coverage.empty:
        rows=[["Analysis family","Data present","Visual present","Coverage"]]
        for item in coverage.itertuples():
            rows.append([item.analysis,"yes" if item.data_available else "no",
                         "yes" if item.figure_available else "no",item.status])
        table(rows,[162,93,93,162])
        missing=coverage[coverage.status.isin(["missing_figure","invalid_figure","invalid_data","missing_data"])]
        if not missing.empty:
            body("WARNING: Analysis outputs without a generated figure: "+
                 ", ".join(missing.analysis.astype(str))+
                 ". Review the figure manifest and plot configuration.")
    metric_audit=metric_coverage(root,cfg,system_repouts)
    if not metric_audit.empty:
        missing_metrics=metric_audit[metric_audit.status!="covered"]
        body(f"Per-metric visual coverage: {len(metric_audit)-len(missing_metrics)} "
             f"of {len(metric_audit)} scalar statistical observables have "
             "their expected standardized time-series figure.")
        if not missing_metrics.empty:
            body("Missing per-metric figures (first 15): "+
                 ", ".join(str(display_cell(metric)) for metric in missing_metrics.metric.head(15))+
                 ". The complete list is in figure_metric_coverage.csv.")
    figures=figure_inventory(root,system_repouts,cfg)
    max_gallery=int(options.get("max_gallery_figures",100))
    if max_gallery<0:
        raise ValueError("report.max_gallery_figures must be >= 0")
    if figures.empty:
        body("No generated scientific figures were found. "
             "Enable plots and run the plotting pipeline before exporting.")
    else:
        primary=figures[figures.kind=="primary"]
        # Deterministic ordering and captions preserve the source location.
        body(f"{len(figures)} generated figures found; "
             f"{len(primary)} are primary analysis figures. "
             "The complete manifest and all original PNGs are in the review ZIP.")
        for item in primary.head(max_gallery).itertuples():
            path_=root/item.figure
            if item.figure in displayed_figures:
                continue
            sub(item.display_label)
            figure(path_,f"Figure ID: {item.figure_id}; source paths in figure_manifest.csv. "
                   "See the accompanying statistical tables for uncertainty "
                   "and the relevant analysis metadata for selections.")
        if len(primary)>max_gallery:
            body(f"{len(primary)-max_gallery} additional primary figures are "
                 "indexed in the complete manifest and preserved in the ZIP. "
                 "Increase report.max_gallery_figures to include more in the PDF.")
        body("Supplementary distribution, autocorrelation and bootstrap-table "
             "figures are preserved in the ZIP rather than repeated in the PDF.")
    story.append(PageBreak())
    header("Files to share for deeper analysis")
    body("Send this PDF for an initial scientific review. For detailed "
         "re-analysis, also send MDscope_review_bundle.zip: it contains "
         "full-precision statistical tables, shared-state tables, topology "
         "selection audits, lipid-map arrays and bounded time-series samples. "
         "It intentionally excludes raw XTC/TPR trajectories and large caches.")
    body("The report is generated deterministically from stored analysis "
         "results, except for the timestamp and any plots whose rendering "
         "libraries differ. It does not generate LLM-written interpretations.")

    # Disabled sections can leave separators; never emit blank separator pages.
    cleaned = []
    for item in story:
        if isinstance(item,PageBreak) and (not cleaned or isinstance(cleaned[-1],PageBreak)):
            continue
        cleaned.append(item)
    while cleaned and isinstance(cleaned[-1],PageBreak):
        cleaned.pop()
    doc.build(cleaned, onFirstPage=footer, onLaterPages=footer)
    return {'figures': sorted(displayed_figures), 'tables': sorted(table_sources)}
