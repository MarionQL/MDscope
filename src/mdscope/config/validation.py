"""Resolved YAML and scientific selection validation (no persistence I/O)."""
from __future__ import annotations

import math
import re
from copy import deepcopy
from pathlib import Path

import numpy as np

from ..errors import InputPathError
from .defaults import DEFAULTS
from .runs import analysis_enabled, analysis_runs

def _validate_time_window(window, label):
    if window is not None and (
            not isinstance(window, (list, tuple)) or len(window) != 2
            or any(type(value) not in (int, float) or not math.isfinite(value) for value in window)
            or window[0] >= window[1]):
        raise ValueError(f"{label} requires two increasing finite times")


def _validate_tail_spec(tail_spec, label):
    if not isinstance(tail_spec, dict) or not (
        tail_spec.get("carbon_pattern") or tail_spec.get("carbon_names")
    ):
        raise ValueError(f"{label} requires carbon_pattern or carbon_names")
    if tail_spec.get("carbon_pattern"):
        try:
            re.compile(str(tail_spec["carbon_pattern"]))
        except re.error as exc:
            raise ValueError(f"{label}.carbon_pattern is not a valid regular expression: {exc}") from exc
    names, positions = tail_spec.get("carbon_names"), tail_spec.get("positions")
    if names and positions is not None and len(names) != len(positions):
        raise ValueError(f"{label}.positions must match carbon_names")


def _validate_tail_naming(order):
    """Tail naming shared by membrane_order and SDC4 local order."""
    tails = order.get("tails", {}) or {}
    if not tails:
        raise ValueError("membrane_order.tails must define at least one tail")
    for tail, tail_spec in tails.items():
        _validate_tail_spec(tail_spec, f"membrane_order.tails.{tail}")
    for resname, overrides in (order.get("lipid_overrides") or {}).items():
        if not isinstance(overrides, dict):
            raise ValueError(f"membrane_order.lipid_overrides.{resname} must map tail names")
        for tail, tail_spec in overrides.items():
            label = f"membrane_order.lipid_overrides.{resname}.{tail}"
            if tail not in tails:
                raise ValueError(f"{label}: tail {tail!r} is not defined in membrane_order.tails")
            if not isinstance(tail_spec, dict):
                raise ValueError(f"{label} must be a mapping")
            if tail_spec.get("carbon_pattern") or tail_spec.get("carbon_names"):
                _validate_tail_spec(tail_spec, label)


def validate_config(cfg: dict, base: Path = Path("."), *, check_paths: bool = True) -> None:
    if analysis_enabled(cfg, "potassium_bridge_audit"):
        analyses = cfg["analyses"]
        solvation = analyses["solvation"]
        motif = solvation["targets"].get("motif", {}).get("acidic_patch", {})
        hs = analyses["sdc4_mechanisms"]["hs_auto"]
        if not (solvation.get("water_selection") and isinstance(motif, dict)
                and motif.get("sequence") and motif.get("selection")
                and hs.get("glycan_selection") and hs.get("attachment_selection")):
            raise ValueError("potassium_bridge_audit requires solvation water and acidic_patch "
                             "selectors and sdc4_mechanisms.hs_auto selections")
    contacts = cfg.get('analyses', {}).get('lipid_contacts', {})
    if contacts.get('enabled', False):
        if not contacts.get('selection') or not contacts.get('species_selections'):
            raise ValueError('lipid_contacts requires selection and species_selections')
        if float(contacts.get('cutoff_A', 5.)) <= 0 or int(contacts.get('batch_atoms', 4096)) < 1:
            raise ValueError('lipid_contacts cutoff_A and batch_atoms must be positive')
    systems = cfg.get("systems")
    if not systems:
        raise ValueError("Config must define at least one system under systems:")
    for field in ("stride", "chunk_frames", "start", "stop"):
        value = cfg["trajectory"].get(field)
        if value is not None and (isinstance(value, bool) or not isinstance(value, int)):
            raise ValueError(f"trajectory.{field} must be an integer")
    if cfg["trajectory"].get("stride", 1) < 1:
        raise ValueError("trajectory.stride must be >= 1")
    if cfg["trajectory"].get("chunk_frames", 1000) < 1:
        raise ValueError("trajectory.chunk_frames must be >= 1")
    start, stop = cfg['trajectory'].get('start'), cfg['trajectory'].get('stop')
    # Preserve supported Python negative indexing, which needs trajectory
    # length to resolve. Reject an obviously empty forward slice early.
    if start is not None and stop is not None and int(start) >= 0 and int(stop) >= 0 and int(stop) <= int(start):
        raise ValueError('trajectory.stop must exceed start for a forward frame window')
    time_spec = cfg.get("time", {}) or {}
    if time_spec.get("dt_ps") is not None:
        dt_ps = float(time_spec["dt_ps"])
        if not np.isfinite(dt_ps) or dt_ps <= 0:
            raise ValueError("time.dt_ps must be a finite value > 0")
    origin_ps = float(time_spec.get("origin_ps", 0.0))
    if not np.isfinite(origin_ps):
        raise ValueError("time.origin_ps must be finite")
    for name, spec in systems.items():
        if not isinstance(name, str):
            raise ValueError("System names must be strings")
        _validate_path_component(name, "system name")
        if not isinstance(spec, dict):
            raise ValueError(f"System {name!r} must be a mapping")
        if not spec.get("topology") or not spec.get("replicates"):
            raise ValueError(
                f"System {name!r} must define topology and at least one replicate"
            )
        top = Path(spec["topology"])
        top = top if top.is_absolute() else (base / top).resolve()
        if check_paths and not top.is_file():
            raise InputPathError(f"System {name!r} topology does not exist: {top}")
        labels = []
        for index, rep in enumerate(spec["replicates"], 1):
            if isinstance(rep, dict):
                if not rep.get("trajectory"):
                    raise ValueError(
                        f"System {name!r} replicate {index} is missing trajectory"
                    )
                trajectory_value = rep["trajectory"]
                label = str(rep.get("label", f"rep{index}"))
            else:
                trajectory_value = rep
                label = f"rep{index}"
            _validate_path_component(
                label,
                f"System {name!r} replicate {index} label",
            )
            labels.append(label)
            traj = Path(trajectory_value)
            traj = traj if traj.is_absolute() else (base / traj).resolve()
            if check_paths and not traj.is_file():
                raise InputPathError(
                    f"System {name!r} replicate {label!r} trajectory does not exist: {traj}"
                )
            if isinstance(rep, dict) and rep.get("msd_trajectory"):
                msd_path = Path(rep["msd_trajectory"])
                if not msd_path.is_absolute():
                    msd_path = (base / msd_path).resolve()
                if check_paths and not msd_path.is_file():
                    raise InputPathError(
                        f"Unmodified membrane MSD trajectory does not exist: {msd_path}"
                    )
        duplicates = sorted({label for label in labels if labels.count(label) > 1})
        if duplicates:
            raise ValueError(
                f"System {name!r} has duplicate replicate labels: "
                + ", ".join(duplicates)
            )

    align = cfg.get("analyses", {}).get("alignment", {})
    if align.get(
        "enabled",
        True
    ) and not (align.get("selection") or cfg.get("selections", {}).get("alignment")):
        raise ValueError("alignment is enabled but no alignment selection was provided")
    for group_id, spec in analysis_groups(cfg).items():
        if not spec.get("selection"):
            raise ValueError(f"selections.groups.{group_id} is missing selection")
    for name in cfg.get("analyses", {}):
        analysis_runs(cfg, name)
    for metric in ("rmsd", "radius_of_gyration"):
        groups = cfg.get("analyses", {}).get(metric, {}).get("groups")
        if groups is not None and (not isinstance(groups, list) or any(not isinstance(x, str) for x in groups)
                                   or set(groups) - set(analysis_groups(cfg))):
            raise ValueError(f"analyses.{metric}.groups must list configured group IDs")
    for run in analysis_runs(cfg, "sasa"):
        spec = run.spec
        if str(spec.get("algorithm", "shrake_rupley")).lower() not in {"shrake_rupley", "shrake-rupley", "sr", "lee_richards", "lee-richards", "lr"}:
            raise ValueError(f"sasa[{run.run_id}].algorithm is unsupported")
        for field in ("n_points", "n_slices", "n_threads"):
            value = spec.get(field)
            if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
                raise ValueError(f"sasa[{run.run_id}].{field} must be a positive integer")
        radius_options = spec.get("element_radii_A") or {}
        if not isinstance(radius_options, dict):
            raise ValueError(f"sasa[{run.run_id}].element_radii_A must be a mapping")
        radii = list(radius_options.values())
        if spec.get("probe_radius_A") is not None:
            radii.append(spec["probe_radius_A"])
        if any(not np.isfinite(float(x)) or float(x) <= 0 for x in radii):
            raise ValueError(f"sasa[{run.run_id}] radii must be finite and positive")

    for run in analysis_runs(cfg, "density_1d"):
        if str(run.spec.get("axis", "z")).lower() not in {"x", "y", "z"}:
            raise ValueError(
                f"density_1d[{run.run_id}].axis must be x, y, or z"
            )
        reference = run.spec.get("membrane_center_selection")
        if reference is not None:
            if (not isinstance(reference, str) or not reference.strip()
                    or run.spec.get("center_selection") or run.run_id == DEFAULT_INSTANCE_ID):
                raise ValueError("Membrane density needs a named run and a single membrane center selection")
            width = float(run.spec.get("bin_width_A", 1.0))
            extent = float(run.spec.get("half_width_A", 0.0))
            if (not np.isfinite(width) or width <= 0 or not np.isfinite(extent)
                    or extent < width or not np.isclose(extent / width, round(extent / width))):
                raise ValueError("Membrane density half_width_A must be a positive multiple of bin_width_A")
        selections = run.spec.get("selections") or {}
        if float(run.spec.get("bin_width_A", 1.0)) <= 0 or not selections:
            raise ValueError(
                f"density_1d[{run.run_id}] requires positive bin_width_A and "
                "non-empty selections"
            )
        if not isinstance(selections, dict):
            raise ValueError(
                f"density_1d[{run.run_id}].selections must be a mapping"
            )
        for selection_id, selection_spec in selections.items():
            _validate_path_component(
                str(selection_id),
                f"density_1d[{run.run_id}] selection ID",
            )
            if isinstance(selection_spec, str):
                continue
            if not isinstance(selection_spec, dict):
                raise ValueError(
                    f"density_1d[{run.run_id}].selections.{selection_id} "
                    "must be a selection string or mapping"
                )
            if not selection_spec.get("selection"):
                raise ValueError(
                    f"density_1d[{run.run_id}].selections.{selection_id} "
                    "is missing selection"
                )
            mode = str(selection_spec.get("mode", "number")).lower()
            if mode not in {"number", "mass"}:
                raise ValueError(
                    f"density_1d[{run.run_id}].selections.{selection_id}.mode "
                    "must be number or mass"
                )

    mem = cfg.get("analyses", {}).get("membrane", {})
    if mem.get("enabled"):
        if str(mem.get("axis", "z")).lower() not in {
            "x",
            "y",
            "z"
        } or not mem.get("headgroup_selection"):
            raise ValueError(
                "enabled membrane analysis requires axis x/y/z and headgroup_selection"
            )
        frac = float((mem.get("msd", {}) or {}).get("slope_start_fraction", 0.2))
        end_frac = float((mem.get("msd", {}) or {}).get("slope_end_fraction", 0.8))
        if not 0 <= frac < end_frac <= 1:
            raise ValueError("membrane.msd slope fractions require 0 <= start < end <= 1")
        multi = mem.get("msd", {}) or {}
        _validate_time_window(multi.get("multi_origin_window_ps"), "membrane.msd.multi_origin_window_ps")
        if multi.get("multi_origin", False):
            if not multi.get("enabled", False):
                raise ValueError("membrane.msd.multi_origin requires msd.enabled")
            lag = multi.get("multi_origin_max_lag_frames")
            if lag is not None and int(lag) < 4:
                raise ValueError("multi_origin_max_lag_frames must be at least 4")
            if int(multi.get("multi_origin_particle_batch", 64)) < 1:
                raise ValueError("multi_origin_particle_batch must be positive")

    order = cfg.get("analyses", {}).get("membrane_order", {})
    if order.get("enabled"):
        if str(order.get("axis", "z")).lower() not in {"x", "y", "z"}:
            raise ValueError("membrane_order.axis must be x, y, or z")
        if float(order.get("hydrogen_cutoff_A", 1.25)) <= 0:
            raise ValueError("membrane_order.hydrogen_cutoff_A must be > 0")
        if int(order.get("min_carbons_per_tail", 4)) < 1:
            raise ValueError("membrane_order.min_carbons_per_tail must be >= 1")
    mechanisms = cfg["analyses"].get("sdc4_mechanisms", {})
    if order.get("enabled") or (mechanisms.get("enabled") and
                                (mechanisms.get("local_order") or {}).get("enabled")):
        _validate_tail_naming(order)

    if cfg["analyses"]["ligand_contacts"].get("enabled") and not cfg["selections"].get("ligand"):
        raise ValueError("ligand_contacts is enabled but selections.ligand is null")
    for run in analysis_runs(cfg, "residue_membrane_distance"):
        if float(run.spec.get("initial_search_cutoff_A", 6.0)) <= 0:
            raise ValueError(
                f"residue_membrane_distance[{run.run_id}]."
                "initial_search_cutoff_A must be > 0"
            )
        if float(run.spec.get("search_growth_factor", 2.0)) <= 1:
            raise ValueError(
                f"residue_membrane_distance[{run.run_id}]."
                "search_growth_factor must be > 1"
            )
        if int(run.spec.get("fallback_block_atoms", 4096)) < 1:
            raise ValueError(
                f"residue_membrane_distance[{run.run_id}]."
                "fallback_block_atoms must be >= 1"
            )
        if not (run.spec.get("membrane_selection") or run.spec.get("lipid_resnames")):
            raise ValueError(
                f"residue_membrane_distance[{run.run_id}] requires "
                "membrane_selection or lipid_resnames"
            )

    for analysis_name in ("lipid_enrichment","bilayer_normal"):
        ext = cfg["analyses"].get(analysis_name, {})
        if ext.get("enabled"):
            if not ext.get("species_selections") or not isinstance(ext.get("species_selections"), dict):
                raise ValueError(f"{analysis_name}.species_selections must define lipid species anchor selections")
            if str(ext.get("axis","z")).lower() not in {"x","y","z"}:
                raise ValueError(f"{analysis_name}.axis must be x, y, or z")
    enrich = cfg["analyses"].get("lipid_enrichment",{})
    if enrich.get("enabled"):
        legacy = {"protein_selection", "inner_radius_A", "outer_radius_A", "area_reference"} & enrich.keys()
        if legacy:
            raise ValueError(f"lipid_enrichment uses footprint_selection and shell_edges_A; remove {sorted(legacy)}")
        if not enrich.get("footprint_selection"):
            raise ValueError("lipid_enrichment.footprint_selection must select the TM/juxtamembrane atoms")
        edges = np.asarray(enrich["shell_edges_A"], float)
        if (edges.ndim != 1 or len(edges) < 2 or edges[0] != 0
                or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0)):
            raise ValueError("lipid_enrichment.shell_edges_A must increase from zero")
        distances = [float(enrich[key]) for key in ("exclusion_radius_A", "slab_half_width_A", "contact_cutoff_A")]
        if not np.isfinite(distances).all() or min(distances) <= 0:
            raise ValueError("lipid_enrichment distances must be finite and positive")
        grid_size = enrich["area_grid_size"]
        if type(grid_size) is not int or not 16 <= grid_size <= 256:
            raise ValueError("lipid_enrichment.area_grid_size must be an integer from 16 to 256")
    normal = cfg["analyses"].get("bilayer_normal",{})
    if normal.get("enabled") and (float(normal.get("neighbor_radius_A",15))<=0 or int(normal.get("min_neighbors",5))<3):
        raise ValueError("bilayer_normal requires positive radius and >=3 neighbors")
    lipid_hb = cfg["analyses"].get("lipid_hydrogen_bonds",{})
    if lipid_hb.get("enabled") and not lipid_hb.get("lipid_selection"):
        raise ValueError("lipid_hydrogen_bonds.lipid_selection is required")

    from ..core.motifs import validate_motifs
    validate_motifs(cfg["selections"].get("motifs", {}))
    local_motifs = cfg["analyses"].get("local_environment", {}).get("motifs")
    if local_motifs is not None:
        if not isinstance(local_motifs, dict):
            raise ValueError("local_environment.motifs must map labels to a regex or {motif: name}")
        for label, value in local_motifs.items():
            if isinstance(value, dict):
                if set(value) != {"motif"} or value["motif"] not in cfg["selections"].get("motifs", {}):
                    raise ValueError(f"local_environment motif {label}: reference must be "
                                     "{motif: <selections.motifs name>}")
            else:
                validate_motifs({label: {"pattern": value}})
    solvation = cfg["analyses"].get("solvation", {})
    if solvation.get("enabled"):
        if solvation.get("preset", "standard") not in ("standard", "detailed", "kinetics"):
            raise ValueError("solvation.preset must be standard, detailed or kinetics")
        if not isinstance(solvation.get("targets"), dict) or not solvation["targets"]:
            raise ValueError("solvation.targets must define explicit chemical-site selections")
        if solvation.get("ion_selection") and solvation.get("ion_selections"):
            raise ValueError("solvation: specify ion_selection OR ion_selections")
        if not (solvation.get("water_selection") or solvation.get("ion_selection")
                or solvation.get("ion_selections")):
            raise ValueError("solvation requires water_selection or ion selections")
        for kind, groups in solvation["targets"].items():
            if not isinstance(groups, dict) or not groups:
                raise ValueError(f"solvation.targets.{kind} must map labels to sites")
            for label, selection in groups.items():
                if isinstance(selection, dict) and "motif" in selection:
                    if (not isinstance(selection["motif"], str)
                            or selection["motif"] not in cfg["selections"].get("motifs", {})):
                        raise ValueError(f"solvation target {kind}/{label}: unknown motif reference")
                    if any(key in selection for key in ("sequence", "selection", "group_by")):
                        raise ValueError("solvation motif references cannot override sequence, selection or group_by")
                    continue
                if isinstance(selection, dict) and "sequence" in selection:
                    import re
                    sequence = selection["sequence"]
                    if not isinstance(sequence, str) or not re.fullmatch(
                        "[ACDEFGHIKLMNPQRSTVWY]+", sequence.upper()
                    ):
                        raise ValueError(
                            f"solvation.targets.{kind}.{label}.sequence must be an exact protein sequence")
                    if "selection" in selection and (
                        not isinstance(selection["selection"], str)
                        or not selection["selection"].strip()
                    ):
                        raise ValueError(
                            f"solvation.targets.{kind}.{label}.selection must be nonempty")
                    if selection.get("group_by", "selection") != "selection":
                        raise ValueError(
                            f"solvation.targets.{kind}.{label}: sequence targets are grouped by motif")
                    continue
                actual = selection if isinstance(selection, str) else (
                    selection.get("selection") if isinstance(selection, dict) else None)
                if not isinstance(actual, str) or not actual.strip():
                    raise ValueError(f"solvation.targets.{kind}.{label} needs a selection")
        if solvation.get("ion_selections") and not isinstance(solvation["ion_selections"], dict):
            raise ValueError("solvation.ion_selections must be a mapping")
        radius = float(solvation.get("search_radius_A", 8.0))
        width = float(solvation.get("rdf_bin_width_A", .2))
        if not (np.isfinite(radius) and np.isfinite(width) and 0 < width <= radius):
            raise ValueError("solvation: invalid search radius or RDF bin width")
        if any(not 0 < float(value) <= radius
               for value in (solvation.get("shell_cutoffs_A") or {}).values()):
            raise ValueError("solvation shell cutoffs must be positive and <= search radius")

    from ..core.glycans import validate_glycan_spec, validate_glycans
    glycans = cfg["selections"].get("glycans", {})
    validate_glycans(glycans)
    hydration = cfg["analyses"].get("hydration", {})
    references = hydration.get("glycans", [])
    if (not isinstance(references, list) or len(set(map(str, references))) != len(references)
            or any(not isinstance(name, str) or name not in glycans for name in references)):
        raise ValueError("hydration.glycans must list distinct selections.glycans names ([] selects none)")
    if hydration.get("enabled"):
        if not (hydration.get("water_oxygen_selection")
                or cfg["analyses"].get("sdc4_mechanisms", {}).get("water_oxygen_selection")
                or solvation.get("water_selection")):
            raise ValueError("hydration requires a water oxygen selection")
        if (not np.isfinite(float(hydration.get("cutoff_A", 3.5)))
                or float(hydration.get("cutoff_A", 3.5)) <= 0):
            raise ValueError("hydration.cutoff_A must be finite and positive")
    mechanisms = cfg["analyses"].get("sdc4_mechanisms", {})
    sugar_cutoff = mechanisms.get("hs_sugar_water_cutoff_A")
    if sugar_cutoff is not None and (
            not mechanisms.get("enabled") or not (mechanisms.get("hs_auto") or {}).get("enabled")
            or not mechanisms.get("water_oxygen_selection")
            or not np.isfinite(float(sugar_cutoff)) or float(sugar_cutoff) <= 0):
        raise ValueError("sdc4_mechanisms.hs_sugar_water_cutoff_A requires enabled SDC4 and hs_auto, water_oxygen_selection and a positive cutoff")
    water_cutoff = mechanisms.get("residue_water_cutoff_A")
    if water_cutoff is not None and (not mechanisms.get("residue_occlusion")
            or not mechanisms.get("water_oxygen_selection")
            or not np.isfinite(float(water_cutoff)) or float(water_cutoff) <= 0):
        raise ValueError("sdc4_mechanisms.residue_water_cutoff_A requires residue_occlusion, water_oxygen_selection and a positive cutoff")
    if not isinstance(mechanisms.get("residue_occlusion", False), bool):
        raise ValueError("sdc4_mechanisms.residue_occlusion must be a boolean")
    surface_sites = mechanisms.get("residue_surface_sites", {})
    if not isinstance(surface_sites, dict) or any(
            not isinstance(name, str) or not name.strip()
            or not isinstance(selection, str) or not selection.strip()
            for name, selection in surface_sites.items()):
        raise ValueError("sdc4_mechanisms.residue_surface_sites must map site names to selections")
    environment_run = mechanisms.get("environment_sasa_run")
    if mechanisms.get("enabled") and environment_run is not None:
        sasa_runs = {run.run_id: run.spec for run in analysis_runs(cfg, "sasa")}
        reference = mechanisms.get("sasa_run", "default")
        if (not isinstance(environment_run, str) or environment_run not in sasa_runs
                or sasa_runs[environment_run].get("reference_run") != reference):
            raise ValueError("sdc4_mechanisms.environment_sasa_run must reference its sasa_run")
    if mechanisms.get("enabled"):
        tm = mechanisms.get("tm_chains") or {}
        ecto = mechanisms.get("ectodomains") or {}
        hs = mechanisms.get("hs_chains") or {}
        sites = mechanisms.get("hs_attachment_sites") or {}
        if not isinstance(tm, dict) or len(tm) != 2 or set(tm) != set(ecto):
            raise ValueError("sdc4_mechanisms needs two matching tm_chains and ectodomains")
        hs_auto = mechanisms.get("hs_auto") or {}
        if hs_auto.get("enabled"):
            if hs or sites:
                raise ValueError("Use hs_auto OR explicit hs_chains/hs_attachment_sites, not both")
            validate_glycan_spec(hs_auto, "sdc4_mechanisms.hs_auto", min_residues=2)
        elif not isinstance(hs, dict) or not isinstance(sites, dict) or set(hs) != set(sites):
            raise ValueError("sdc4_mechanisms.hs_chains and hs_attachment_sites must match")
        for group in (tm, ecto, hs, sites):
            for key, selection in group.items():
                _validate_path_component(str(key), "sdc4 mechanism group")
                if not isinstance(selection, str) or not selection.strip():
                    raise ValueError("SDC4 group selections must be nonempty strings")
        if not mechanisms.get("lipid_anchors") or not mechanisms.get("lipid_contact_selection"):
            raise ValueError("sdc4_mechanisms requires lipid_anchors and lipid_contact_selection")
        if not mechanisms.get("surface_anchor_species") or not isinstance(mechanisms["surface_anchor_species"], list):
            raise ValueError("sdc4_mechanisms.surface_anchor_species must be a nonempty list")
        if (hs or hs_auto.get("enabled")) and not mechanisms.get("water_oxygen_selection"):
            raise ValueError("HS hydration requires water_oxygen_selection")
        if str(mechanisms.get("axis", "z")).lower() != "z":
            raise ValueError("sdc4_mechanisms supports z-normal membranes only")
        if float(mechanisms.get("map_extent_A", 40)) <= 0 or float(mechanisms.get("map_bin_A", 2)) <= 0:
            raise ValueError("sdc4_mechanisms map extent and bin size must be positive")
        if int(mechanisms.get("local_min_lipids", 8)) < 8:
            raise ValueError("sdc4_mechanisms.local_min_lipids must be >= 8")
        if mechanisms.get("local_order", {}).get("enabled") and not mechanisms["local_order"].get("lipid_selection"):
            raise ValueError("sdc4_mechanisms.local_order.lipid_selection is required")
    from ..state_feature_spec import validate_state_features
    shared = cfg["analyses"].get("shared_states", {})
    validate_state_features(shared)
    if shared.get("enabled") and shared.get("preset", "sdc4") == "sdc4" and not mechanisms.get("enabled"):
        raise ValueError("shared_states preset sdc4 requires sdc4_mechanisms.enabled")
    allowed = {"time_series", "profiles", "interactions", "solvation", "sdc4",
               "shared_states", "diffusion", "lipid_enrichment", "lipid_contacts", "other"}
    plots = cfg.get("plots", {})
    if not isinstance(plots.get("table_reports_pdf", False), bool):
        raise ValueError("plots.table_reports_pdf must be boolean")
    window_ns = plots.get("time_series_rolling_window_ns")
    if window_ns is not None and (
        type(window_ns) not in (int, float)
        or not math.isfinite(window_ns) or window_ns <= 0
    ):
        raise ValueError("plots.time_series_rolling_window_ns must be null or a positive finite number")
    families = plots.get("families", {})
    if (not isinstance(families, dict) or set(families) - allowed
            or any(not isinstance(value, bool) for value in families.values())):
        raise ValueError("plots.families must map known families to booleans")
    from ..reporting.plot_ids import validate_figure_options, validate_plot_ids
    validate_plot_ids(plots.get("include"))
    validate_figure_options(plots.get("options", {}))
    excluded = plots.get("interaction_exclude_kinds", [])
    if not isinstance(excluded, list) or any(not isinstance(value, str) for value in excluded):
        raise ValueError("plots.interaction_exclude_kinds must be a list")
    stats = cfg.get("statistics", {})
    residue_response_settings(cfg)
    from ..statistics_selection import validate_metric_selection
    validate_metric_selection(stats)
    for field in ("block", "hypothesis_tests"):
        if stats.get(field) is not None and not isinstance(stats[field], dict):
            raise ValueError(f"statistics.{field} must be a mapping")
    if int(stats.get("bootstrap_draws", 2000)) < 100:
        raise ValueError("statistics.bootstrap_draws must be >= 100")
    if not 0 < float(stats.get("ci", 0.95)) < 1:
        raise ValueError("statistics.ci must be between 0 and 1")
    for key in ("interaction_include_kinds", "interaction_exclude_kinds"):
        value = stats.get(key)
        if value is not None and (not isinstance(value, list) or any(not isinstance(x, str) for x in value)):
            raise ValueError(f"statistics.{key} must be a list of kinds")
    if int(stats.get("interaction_progress_every", 1000)) < 1:
        raise ValueError("statistics.interaction_progress_every must be positive")
    _validate_time_window(stats.get("scalar_window_ps"), "statistics.scalar_window_ps")
    if stats.get("matched_window_ps") != "auto":
        _validate_time_window(stats.get("matched_window_ps"), "statistics.matched_window_ps")
    families = stats.get("families", {})
    if not isinstance(families, dict) or set(families) - {"time_series", "profiles", "interactions", "solvation_rdfs", "diffusion", "coupling"}:
        raise ValueError("statistics.families has unknown families or is not a mapping")
    for name, option in families.items():
        if isinstance(option, bool):
            continue
        if not isinstance(option, dict) or set(option) - (
                {"enabled", "bootstrap", "bootstrap_draws", "ci", "hypothesis_tests", "hierarchy"}
                | ({"adaptive_hierarchy"} if name == "time_series" else set())):
            raise ValueError(f"statistics.families.{name}: unsupported options")
        if "adaptive_hierarchy" in option:
            adaptive = option["adaptive_hierarchy"]
            counts = {"min_draws", "max_draws", "batch_size", "stable_checks"}
            tolerances = {"absolute_tolerance", "relative_tolerance"}
            if (not isinstance(adaptive, dict) or not counts <= adaptive.keys()
                    or adaptive.keys() - counts - tolerances
                    or any(type(adaptive[k]) is not int or adaptive[k] < 1 for k in counts)
                    or not 100 <= adaptive["min_draws"] <= adaptive["max_draws"]
                    or any(type(adaptive.get(k, 0)) not in (int, float)
                           or not math.isfinite(adaptive.get(k, 0)) or adaptive.get(k, 0) < 0
                           for k in tolerances)
                    or not any(adaptive.get(k, 0) > 0 for k in tolerances)):
                raise ValueError("statistics.families.time_series.adaptive_hierarchy requires valid "
                                 "draw bounds, batch_size, stable_checks and a positive tolerance")
        for key in ("enabled", "bootstrap", "hypothesis_tests"):
            if key in option and not isinstance(option[key], bool):
                raise ValueError(f"statistics.families.{name}.{key} must be boolean")
        if "bootstrap_draws" in option and int(option["bootstrap_draws"]) < 100:
            raise ValueError(f"statistics.families.{name}.bootstrap_draws must be >= 100")
        if "ci" in option and not 0 < float(option["ci"]) < 1:
            raise ValueError(f"statistics.families.{name}.ci must be between 0 and 1")
        hierarchy = option.get("hierarchy", {})
        if (not isinstance(hierarchy, dict) or set(hierarchy) - {"enabled", "block_within_replicate", "resample_replicates"}
                or any(not isinstance(v, bool) for v in hierarchy.values())):
            raise ValueError(f"statistics.families.{name}.hierarchy must contain boolean controls")
    if not isinstance(stats.get("bootstrap", True), bool):
        raise ValueError("statistics.bootstrap must be boolean")
    hierarchy = stats.get("hierarchy", {})
    if (not isinstance(hierarchy, dict) or set(hierarchy) - {"enabled", "block_within_replicate", "resample_replicates"}
            or any(not isinstance(v, bool) for v in hierarchy.values())):
        raise ValueError("statistics.hierarchy must contain boolean controls")
    hypotheses = stats.get("hypothesis_families", {})
    if not isinstance(hypotheses, dict):
        raise ValueError("statistics.hypothesis_families must map names to metric lists")
    seen = set()
    for name, metrics in hypotheses.items():
        if not isinstance(name, str) or not name.strip() or name == "exploratory" or not isinstance(metrics, list) or not metrics or any(not isinstance(x, str) for x in metrics):
            raise ValueError("hypothesis families require named nonempty metric lists; exploratory is reserved")
        if seen.intersection(metrics) or len(set(metrics)) != len(metrics):
            raise ValueError("A metric cannot belong to multiple hypothesis families")
        seen.update(metrics)
    block = stats.get("block", {}) or {}
    if str(block.get("method", "moving")) != "moving":
        raise ValueError("statistics.block.method currently supports only 'moving'")
    length = block.get("length", "auto")
    if length != "auto" and int(length) < 1:
        raise ValueError("statistics.block.length must be 'auto' or an integer >= 1")
    if float(block.get("tau_multiplier", 2.0)) <= 0:
        raise ValueError("statistics.block.tau_multiplier must be > 0")
    max_acf_lag = block.get("max_acf_lag")
    if max_acf_lag is not None and int(max_acf_lag) < 1:
        raise ValueError("statistics.block.max_acf_lag must be null or >= 1")
    tests = stats.get("hypothesis_tests", {}) or {}
    tests_active = tests.get("enabled", True) or any(
        isinstance(option, dict) and option.get("enabled", True) and option.get("hypothesis_tests", False)
        for option in families.values())
    if tests_active:
        design = str(tests.get("design", "independent"))
        if design not in {"independent", "paired"}:
            raise ValueError(
                "statistics.hypothesis_tests.design must be independent or "
                "paired"
            )
        if design == "paired" and len(cfg.get("systems", {})) > 1:
            label_sets = []
            for system in cfg["systems"]:
                labels = {
                    spec["id"]
                    for spec in replicate_specs(cfg, system)
                }
                label_sets.append((system, labels))
            reference_system, reference_labels = label_sets[0]
            for system, labels in label_sets[1:]:
                if labels != reference_labels:
                    raise ValueError(
                        "paired hypothesis tests require identical replicate "
                        f"labels across systems; {reference_system!r} and "
                        f"{system!r} differ"
                    )
        if str(tests.get("method", "permutation")) != "permutation":
            raise ValueError(
                "statistics.hypothesis_tests.method currently supports only permutation"
            )
        if int(tests.get("permutation_draws", 10000)) < 100:
            raise ValueError(
                "statistics.hypothesis_tests.permutation_draws must be >= 100"
            )
        if int(tests.get("exact_max_assignments", 100000)) < 1:
            raise ValueError(
                "statistics.hypothesis_tests.exact_max_assignments must be >= 1"
            )
        if str(tests.get("multiple_testing", "fdr_bh")) not in {"fdr_bh", "none"}:
            raise ValueError("statistics.hypothesis_tests.multiple_testing must be fdr_bh or none")





def residue_response_settings(cfg):
    """Resolve and validate the saved-data exploratory screen, independently of UI."""
    supplied = cfg.get("statistics", {}).get("residue_response", {})
    if not isinstance(supplied, dict):
        raise ValueError("statistics.residue_response must be a mapping")
    defaults = DEFAULTS["statistics"]["residue_response"]
    unknown = set(supplied) - set(defaults)
    if unknown:
        raise ValueError(f"Unknown statistics.residue_response options: {sorted(unknown)}")
    settings = {**deepcopy(defaults), **deepcopy(supplied)}
    if type(settings["enabled"]) is not bool:
        raise ValueError("statistics.residue_response.enabled must be boolean")
    for key in ("profiles", "chains"):
        values = settings[key]
        if values is not None and (not isinstance(values, list) or
                any(not isinstance(v, str) or not v.strip() for v in values) or
                len(set(values)) != len(values)):
            raise ValueError(f"statistics.residue_response.{key} must be null or distinct nonempty strings")
    from ..statistics_contracts import RESIDUE_RESPONSE_PROFILE
    unsupported = [name for name in settings["profiles"] or []
                   if not RESIDUE_RESPONSE_PROFILE.match(name)]
    if unsupported:
        raise ValueError(f"Unsupported statistics.residue_response profiles: {unsupported}")
    length = settings["min_length"]
    if type(length) is not int or length < 2:
        raise ValueError("statistics.residue_response.min_length must be an integer >= 2")
    cosine = settings["minimum_cosine"]
    if (isinstance(cosine, bool) or not isinstance(cosine, (int, float)) or
            not math.isfinite(cosine) or not 0 <= cosine <= 1):
        raise ValueError("statistics.residue_response.minimum_cosine must be finite between 0 and 1")
    return settings


