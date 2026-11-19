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


def review(document, base, choices=None):
    """Return validated YAML + production plan, filtering inactive figure choices.

    choices maps analysis names or 'statistics' to exact plot IDs. The editable
    choices live outside this function and survive disabling their owner.
    """
    raw = deepcopy(document)
    cfg = resolve_config(raw, base=Path(base))
    if choices is not None:
        from ..registry import STRUCTURAL_SCALAR_METRICS, structural_scalar_contracts
        variables = structural_scalar_contracts(cfg)
        eligible = {
            name: {fig["id"] for var in variables if var["analysis"] == name
                   for fig in var["figures"] if not fig["requires_statistics"]}
            for name in ("rmsd", "radius_of_gyration")
        }
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
        def managed(identifier):
            parts = identifier.split(":")
            return len(parts) == 3 and parts[0] == "time_series" and any(
                parts[1] == metric or parts[1].startswith(metric + "__group_")
                for metric in STRUCTURAL_SCALAR_METRICS.values())
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
