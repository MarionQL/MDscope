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


def change_group(document, old_id, new_id=None, *, duplicate=False, choices=None):
    """Rename/duplicate/delete while migrating only actual shared-group callers."""
    from ..config import merge
    from ..planning import authoring_catalog
    before = authoring_catalog(merge(DEFAULTS, document), include_inactive=True)
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
        if "time_series" in selection:
            selection["time_series"] = [
                {**item, "group": new_id} if item.get("group") == old_id else item
                for item in selection["time_series"]
                if item.get("group") != old_id or new_id is not None]

    if not duplicate:
        after = authoring_catalog(merge(DEFAULTS, document), include_inactive=True)
        old_profiles = [item["profile"] for item in before["profiles"]]
        new_profiles = [item["profile"] for item in after["profiles"]]
        removed = [item for item in old_profiles if item not in new_profiles]
        added = [item for item in new_profiles if item not in old_profiles]
        mapping = dict(zip(removed, added)) if new_id else dict.fromkeys(removed)
        if "profiles" in selection:
            selection["profiles"] = [dict(item, profile=mapping.get(item["profile"], item["profile"]))
                                     for item in selection["profiles"]
                                     if item["profile"] not in mapping or mapping[item["profile"]] is not None]
        from ..reporting.plot_ids import safe_plot_stem
        stems = {safe_plot_stem(old): safe_plot_stem(new) if new else None for old, new in mapping.items()}
        for variable in before["scalars"]:
            if variable["group"] == old_id:
                for figure in variable["figures"]:
                    stem = figure["id"].split(":")[1]
                    replacement = next((item for item in after["scalars"]
                                        if item["analysis"] == variable["analysis"] and item["group"] == new_id), None) if new_id else None
                    stems[stem] = replacement["figures"][0]["id"].split(":")[1] if replacement else None
        def migrate(ids):
            result = []
            for identifier in ids:
                parts = identifier.split(":")
                if len(parts) > 2 and parts[0] in {"profiles", "time_series"} and parts[1] in stems:
                    replacement = stems[parts[1]]
                    if replacement is None:
                        continue
                    parts[1] = replacement
                    identifier = ":".join(parts)
                result.append(identifier)
            return result
        plots = document.get("plots", {})
        if plots.get("include") is not None:
            plots["include"] = migrate(plots["include"])
        if choices is not None:
            for owner in choices:
                choices[owner] = migrate(choices[owner])


def review(document, base, choices=None):
    """Return validated YAML + production plan, filtering inactive figure choices.

    choices maps analysis names or 'statistics' to exact plot IDs. The editable
    choices live outside this function and survive disabling their owner.
    """
    raw = deepcopy(document)
    cfg = resolve_config(raw, base=Path(base))
    if choices is not None:
        from ..planning import authoring_catalog
        catalog = authoring_catalog(cfg, include_inactive=True)
        active = authoring_catalog(cfg)
        variables = active["scalars"]
        active_scalars = {(var["metric"], var["group"]) for var in variables}
        known_scalars = {(var["metric"], var["group"]) for var in catalog["scalars"]}
        active_profiles = {var["profile"] for var in active["profiles"]}
        known_profiles = {var["profile"] for var in catalog["profiles"]}
        selectors = raw.get("statistics", {}).get("metric_selection")
        if selectors is not None:
            for family, known, active_identities, identity in (
                ("time_series", known_scalars, active_scalars, lambda item: (item["metric"], item.get("group"))),
                ("profiles", known_profiles, active_profiles, lambda item: item["profile"]),
            ):
                if family in selectors:
                    selectors[family] = [item for item in selectors[family]
                                         if identity(item) not in known or identity(item) in active_identities]
            cfg = resolve_config(raw, base=Path(base))
        eligible = {owner: {figure["id"] for figure in figures}
                    for owner, figures in authoring_catalog(cfg)["figures"].items()}
        ids = []
        for owner, selected in choices.items():
            active = (cfg["statistics"].get("enabled") if owner == "statistics"
                      else cfg["analyses"].get(owner, {}).get("enabled"))
            if active:
                ids.extend(identifier for identifier in selected
                           if owner not in eligible or identifier in eligible[owner])
        # Exact mode must not silently carry unowned imported figures forward.
        original = ((raw.get("plots", {}).get("include") or [])
                    if raw.get("plots", {}).get("enabled", False) else [])
        managed_ids = {fig["id"] for figures in catalog["figures"].values() for fig in figures}
        unknown = [identifier for identifier in original if identifier not in managed_ids]
        if unknown:
            raise ValueError("Imported figures lack authoring ownership contracts: " + ", ".join(unknown) +
                             ". Keep legacy figure mode or explicitly remove these IDs before exact export.")
        plots = raw.setdefault("plots", {})
        plots.update(enabled=bool(ids), include=list(dict.fromkeys(ids)))
        # Exact figure choices supersede imported legacy switches that would
        # otherwise suppress their production renderers. Unrelated switches
        # and the original editable document remain unchanged.
        for family in {identifier.split(":", 1)[0] for identifier in ids}:
            if family in plots:
                plots[family] = True
            if family in plots.get("families", {}):
                plots["families"][family] = True
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


def delete_named_run(document, name, run_id, *, choices=None):
    """Delete a run and its catalogued selectors; reject dangling references."""
    from ..config import MULTI_RUN_ANALYSES
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
    if not isinstance(mechanism, dict):
        mechanism = {}
    if name == "sasa":
        references += [f"sdc4_mechanisms.{key}" for key in ("sasa_run", "environment_sasa_run")
                       if mechanism.get(key) == run_id]
    if references:
        raise ValueError(f"Run {run_id!r} is referenced by {', '.join(references)}; change those references before deletion")
    _discard_run_outputs(document, name, {run_id}, choices)
    del runs[run_id]


def clear_named_runs(document, name, *, choices=None):
    """Return to the default run without retaining removed-run selections."""
    from ..config import MULTI_RUN_ANALYSES
    if name not in MULTI_RUN_ANALYSES:
        raise ValueError(f"{name} does not support named runs")
    spec = document["analyses"][name]
    run_ids = set(spec.get("runs") or {})
    mechanism = document["analyses"].get("sdc4_mechanisms", {})
    if name == "sasa" and isinstance(mechanism, dict):
        references = [key for key in ("sasa_run", "environment_sasa_run") if mechanism.get(key) in run_ids]
        if references:
            raise ValueError("Named SASA runs are referenced by sdc4_mechanisms: " + ", ".join(references))
    _discard_run_outputs(document, name, run_ids, choices)
    spec.pop("runs", None)


def _discard_run_outputs(document, name, run_ids, choices):
    from ..config import merge
    from ..planning import authoring_catalog
    from ..reporting.plot_ids import safe_plot_stem
    catalog = authoring_catalog(merge(DEFAULTS, document), include_inactive=True)
    profiles = {item["profile"] for item in catalog["profiles"]
                if item["analysis"] == name and item["run_id"] in run_ids}
    selectors = document.get("statistics", {}).get("metric_selection") or {}
    if "profiles" in selectors:
        selectors["profiles"] = [item for item in selectors["profiles"] if item["profile"] not in profiles]
    plots = document.get("plots", {})
    stems = {safe_plot_stem(profile) for profile in profiles}
    if name == "residue_membrane_distance":
        stems.update(f"residue_membrane_distance__{safe_plot_stem(run_id)}" for run_id in run_ids)
    def retain(identifiers):
        return [identifier for identifier in identifiers
                if not (identifier.startswith("profiles:") and identifier.split(":")[1] in stems)]
    if plots.get("include") is not None:
        plots["include"] = retain(plots["include"])
    if choices is not None:
        for owner in choices:
            choices[owner] = retain(choices[owner])


def commit_draft(document, draft):
    """Commit editable input atomically; incomplete drafts remain editable."""
    replacement = deepcopy(draft)
    if not isinstance(replacement, dict):
        raise ValueError("Authoring document must be a mapping")
    document.clear()
    document.update(replacement)


def add_named_run(document, name, run_id):
    """Validate run identity before mutating the draft."""
    from ..config import MULTI_RUN_ANALYSES, analysis_runs, merge
    if name not in MULTI_RUN_ANALYSES:
        raise ValueError(f"{name} does not support named runs")
    runs = document["analyses"][name].get("runs", {})
    if not run_id or run_id in runs:
        raise ValueError("Enter a unique run ID")
    candidate = deepcopy(document)
    candidate["analyses"][name].setdefault("runs", {})[run_id] = {}
    candidate["analyses"][name]["enabled"] = True
    analysis_runs(merge(DEFAULTS, candidate), name)
    document["analyses"][name].setdefault("runs", {})[run_id] = {}


def add_group(document, group_id):
    from ..config import analysis_groups
    groups = document.setdefault("selections", {}).setdefault("groups", {})
    if not group_id or group_id in groups:
        raise ValueError("Enter a unique filesystem-safe group ID")
    candidate = deepcopy(document)
    candidate["selections"]["groups"][group_id] = {"selection": "all"}
    analysis_groups(candidate)
    groups[group_id] = {"label": group_id, "selection": "all"}


def add_system(document, system_id):
    """Add a condition using the backend's identifier contract."""
    from ..config import _validate_path_component
    _validate_path_component(system_id, "system name")
    systems = document.setdefault("systems", {})
    if not system_id or system_id in systems:
        raise ValueError("Enter a unique filesystem-safe condition ID")
    systems[system_id] = {"topology": "topology.tpr", "replicates": [
        {"label": "rep1", "trajectory": "rep1.xtc"}]}


def delete_system(document, system_id, *, choices=None):
    """Remove condition-owned figures without changing shared figure choices."""
    from ..config import merge
    from ..planning import authoring_catalog
    systems = document["systems"]
    if system_id not in systems:
        raise ValueError(f"Unknown condition {system_id!r}")
    if len(systems) == 1:
        raise ValueError("Keep at least one condition")
    def identifiers(doc):
        catalog = authoring_catalog(merge(DEFAULTS, doc), include_inactive=True)
        return {fig["id"] for figures in catalog["figures"].values() for fig in figures}
    candidate = deepcopy(document)
    del candidate["systems"][system_id]
    removed = identifiers(document) - identifiers(candidate)
    del systems[system_id]
    plots = document.get("plots", {})
    if plots.get("include") is not None:
        plots["include"] = [key for key in plots["include"] if key not in removed]
    if choices is not None:
        for owner in choices:
            choices[owner] = [key for key in choices[owner] if key not in removed]
