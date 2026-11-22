"""Validation rules for membrane configuration."""
import numpy as np
from ..runs import analysis_runs
from .common import _validate_time_window, _validate_tail_naming


def validate_membrane(cfg):
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
