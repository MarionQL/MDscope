"""Configuration authoring without Streamlit, trajectory reads or result writes."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import yaml

from ..config import DEFAULTS, resolve_config
from ..config_schema import option_schema
from ..planning import build_plan
from ..statistics_selection import studio_metric_selection_defaults


def analysis_options():
    """Supported option names come from the CLI schema, never a GUI copy."""
    return option_schema(DEFAULTS)["analyses"]


def new_document():
    return {
        "project": {"name": "mdscope_project", "output": "mdscope_results"},
        "systems": {"condition1": {"topology": "topology.tpr", "replicates": [
            {"label": "rep1", "trajectory": "rep1.xtc"}]}},
        "selections": {"groups": {}},
        "analyses": {name: {"enabled": False} for name in analysis_options()},
        "statistics": {"enabled": False,
                       "metric_selection": studio_metric_selection_defaults()},
        "plots": {"enabled": False, "include": []},
    }


def import_document(text):
    raw = yaml.safe_load(text)
    resolve_config(raw)
    return raw


def change_group(document, old_id, new_id=None, *, duplicate=False):
    """Rename/duplicate/delete while migrating only actual shared-group callers."""
    groups = document.setdefault("selections", {}).setdefault("groups", {})
    if new_id == old_id:
        return
    if new_id is not None:
        if new_id in groups:
            raise ValueError(f"Group {new_id!r} already exists")
        candidate = deepcopy(document)
        candidate["selections"]["groups"][new_id] = deepcopy(groups[old_id])
        # Use the backend's filesystem-safe group ID contract.
        from ..config import analysis_groups
        analysis_groups(candidate)
    replaced = {}
    for key, spec in groups.items():
        if key != old_id or duplicate:
            replaced[key] = spec
        if key == old_id and new_id is not None:
            replaced[new_id] = deepcopy(spec)
    document["selections"]["groups"] = replaced
    if not duplicate:
        for name in ("rmsd", "radius_of_gyration"):
            spec = document.get("analyses", {}).get(name, {})
            if isinstance(spec, dict) and spec.get("groups") is not None:
                spec["groups"] = [new_id if key == old_id else key
                                  for key in spec["groups"]
                                  if key != old_id or new_id is not None]

    if not duplicate:
        selection = document.get("statistics", {}).get("metric_selection") or {}
        for family, field, old_value, new_value in (
            ("time_series", "group", old_id, new_id),
            ("profiles", "profile", f"rmsf__group_{old_id}",
             f"rmsf__group_{new_id}" if new_id else None),
        ):
            if family in selection:
                updated = []
                for entry in selection[family]:
                    if entry.get(field) == old_value:
                        if new_value is None:
                            continue
                        entry = {**entry, field: new_value}
                    updated.append(entry)
                selection[family] = updated


def review(document, base, choices=None):
    """Return validated YAML + production plan, filtering inactive figure choices.

    choices maps analysis names or 'statistics' to exact plot IDs. The editable
    choices live outside this function and survive disabling their owner.
    """
    raw = deepcopy(document)
    cfg = resolve_config(raw, base=Path(base))
    if choices is not None:
        from ..registry import STRUCTURAL_SCALAR_METRICS, scalar_contracts
        from ..registry import configured_profile_variables
        catalog_cfg = deepcopy(cfg)
        for name in ("rmsd", "radius_of_gyration", "rmsf", "membrane", "dssp", "sasa", "density_1d", "residue_membrane_distance"):
            catalog_cfg["analyses"][name]["enabled"] = True
        for spec in catalog_cfg["analyses"].values():
            for run_id, run in (spec.get("runs") or {}).items():
                spec["runs"][run_id] = {**run, "enabled": True} if isinstance(run, dict) else {"enabled": True}
        catalog_variables = scalar_contracts(catalog_cfg)
        catalog_profiles = configured_profile_variables(catalog_cfg)
        variables = scalar_contracts(cfg)
        active_scalars = {(var["metric"], var["group"]) for var in variables}
        known_scalars = {(var["metric"], var["group"]) for var in catalog_variables}
        active_profiles = {var["profile"] for var in configured_profile_variables(cfg)}
        known_profiles = {var["profile"] for var in catalog_profiles}
        selectors = raw.get("statistics", {}).get("metric_selection")
        if selectors is not None:
            for family, known, active, identity in (
                ("time_series", known_scalars, active_scalars, lambda item: (item["metric"], item.get("group"))),
                ("profiles", known_profiles, active_profiles, lambda item: item["profile"]),
            ):
                if family in selectors:
                    selectors[family] = [item for item in selectors[family]
                                         if identity(item) not in known or identity(item) in active]
            cfg = resolve_config(raw, base=Path(base))
        eligible = {
            name: {fig["id"] for var in variables if var["analysis"] == name
                   for fig in var["figures"] if not fig["requires_statistics"]}
            for name in {var["analysis"] for var in variables}
        }
        from ..planning import analysis_figure_candidates
        for name in cfg["analyses"]:
            eligible.setdefault(name, set()).update(fig["id"] for fig in analysis_figure_candidates(cfg, name))
        selected_stats = cfg["statistics"].get("metric_selection", {}).get("time_series")
        eligible["statistics"] = {
            fig["id"] for var in variables
            if selected_stats is None or any(
                item["metric"] == var["metric"] and item.get("group") == var["group"]
                for item in selected_stats)
            for fig in var["figures"] if fig["requires_statistics"]
        }
        from ..planning import statistical_figure_candidates
        eligible["statistics"].update(fig["id"] for fig in statistical_figure_candidates(cfg))
        ids = []
        for owner, selected in choices.items():
            active = (cfg["statistics"].get("enabled") if owner == "statistics"
                      else cfg["analyses"].get(owner, {}).get("enabled"))
            if active:
                ids.extend(identifier for identifier in selected
                           if owner not in eligible or identifier in eligible[owner])
        # Preserve imported IDs outside the guided scalar catalog. They remain
        # explicit configuration, with prerequisites visible in the review plan.
        original = ((raw.get("plots", {}).get("include") or [])
                    if raw.get("plots", {}).get("enabled", False) else [])
        managed_ids = {fig["id"] for var in catalog_variables for fig in var["figures"]}
        from ..reporting.plot_ids import safe_plot_stem
        managed_ids.update(f"profiles:{safe_plot_stem(var['profile'])}:{variant}"
                           for var in catalog_profiles for variant in
                           ("profile", "distribution", "system_difference", "comparison_table", "bootstrap_table"))
        def managed(identifier):
            return identifier in managed_ids
        ids.extend(identifier for identifier in original if not managed(identifier))
        plots = raw.setdefault("plots", {})
        plots.update(enabled=bool(ids), include=list(dict.fromkeys(ids)))
        for variant, flag in (("distribution", "distributions"),
                              ("autocorrelation", "autocorrelation"),
                              ("bootstrap_table", "bootstrap_tables")):
            if any(identifier.endswith(":" + variant) for identifier in ids):
                plots[flag] = True
        cfg = resolve_config(raw, base=Path(base))
    plan = build_plan(cfg, Path(base) / "analysis.yaml")
    return yaml.safe_dump(raw, sort_keys=False, allow_unicode=True), plan


def topology_preflight(document, base):
    """Read topology only through production preflight; never open trajectories."""
    from ..core.validation import validate_selections
    cfg = resolve_config(document, base=Path(base))
    return validate_selections(cfg, Path(base))


def statistics_selector_fields(family):
    """Guided table columns from backend identity and estimator contracts."""
    from ..statistics_selection import SIMPLE_FAMILY_KEYS, ESTIMATOR_FLAGS
    keys = (("kind", "interaction", "metric") if family == "interactions"
            else SIMPLE_FAMILY_KEYS[family])
    return [*keys, "inference", *sorted(ESTIMATOR_FLAGS[family])]


def delete_named_run(document, name, run_id):
    """Delete a run and its catalogued selectors; reject dangling references."""
    from ..config import MULTI_RUN_ANALYSES, merge
    from ..registry import configured_profile_variables
    from ..reporting.plot_ids import safe_plot_stem
    if name not in MULTI_RUN_ANALYSES:
        raise ValueError(f"{name} does not support named runs")
    spec = document["analyses"][name]
    runs = spec["runs"]
    if run_id not in runs:
        raise ValueError(f"Unknown run {run_id!r}")
    if len(runs) == 1:
        raise ValueError("Last run: turn off Configure named runs to return to the default analysis")
    references = [key for key, value in runs.items() if key != run_id
                  and isinstance(value, dict) and value.get("reference_run") == run_id]
    mechanism = document["analyses"].get("sdc4_mechanisms", {})
    if name == "sasa":
        references += [f"sdc4_mechanisms.{key}" for key in ("sasa_run", "environment_sasa_run")
                       if mechanism.get(key) == run_id]
    if references:
        raise ValueError(f"Run {run_id!r} is referenced by {', '.join(references)}; change those references before deletion")
    catalog = merge(DEFAULTS, document)
    catalog["analyses"][name]["enabled"] = True
    run = catalog["analyses"][name]["runs"][run_id]
    catalog["analyses"][name]["runs"][run_id] = {**run, "enabled": True} if isinstance(run, dict) else {"enabled": True}
    profiles = {item["profile"] for item in configured_profile_variables(catalog)
                if item["analysis"] == name and item["run_id"] == run_id}
    del runs[run_id]
    selectors = document.get("statistics", {}).get("metric_selection") or {}
    if "profiles" in selectors:
        selectors["profiles"] = [item for item in selectors["profiles"] if item["profile"] not in profiles]
    plots = document.get("plots", {})
    if plots.get("include") is not None:
        stems = {safe_plot_stem(profile) for profile in profiles}
        if name == "residue_membrane_distance":
            stems.add(f"residue_membrane_distance__{safe_plot_stem(run_id)}")
        plots["include"] = [identifier for identifier in plots["include"]
                            if not (identifier.startswith("profiles:") and identifier.split(":")[1] in stems)]
