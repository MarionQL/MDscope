"""Validation rules for structure configuration."""
import numpy as np
from ..defaults import DEFAULT_INSTANCE_ID
from ..runs import analysis_groups, analysis_runs, _validate_path_component


def validate_structure(cfg):
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
