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



def _atomic_edit(document, choices, operation, *args, **kwargs):
    """Apply a lifecycle edit only after all validation/reconciliation succeeds."""
    draft = deepcopy(document)
    draft_choices = deepcopy(choices) if choices is not None else None
    operation(draft, *args, choices=draft_choices, **kwargs)
    commit_draft(document, draft)
    if choices is not None:
        choices.clear()
        choices.update(draft_choices)


def replace_global_selections(document, edited):
    """Replace global keys while retaining named groups; fail without mutation."""
    if not isinstance(edited, dict):
        raise ValueError("Selections must be a mapping")
    if "groups" in edited:
        raise ValueError("Use the named group editor for groups")
    candidate = deepcopy(document)
    previous = candidate.setdefault("selections", {}).get("groups", {})
    candidate["selections"] = {**deepcopy(edited), "groups": deepcopy(previous)}
    commit_draft(document, candidate)

def set_figure_selection(document, catalog, identifiers):
    """Commit exact hierarchy choices using the existing catalog's ownership."""
    selected = set(identifiers)
    known = {
        figure["id"] for figures in catalog.values()
        for figure in figures
    }
    unknown = sorted(selected - known)
    if unknown:
        raise ValueError("Unknown figure IDs in hierarchy: " + ", ".join(unknown))
    choices = {
        owner: [figure["id"] for figure in figures
                if figure["id"] in selected]
        for owner, figures in catalog.items()
    }
    ordered = list(dict.fromkeys(
        identifier for owned in choices.values() for identifier in owned
    ))
    candidate = deepcopy(document)
    candidate.setdefault("plots", {})["include"] = ordered
    candidate["plots"]["enabled"] = bool(ordered)
    commit_draft(document, candidate)
    return choices


def change_group(document, old_id, new_id=None, *, duplicate=False, choices=None):
    return _atomic_edit(document, choices, _change_group, old_id, new_id,
                        duplicate=duplicate)


def delete_named_run(document, name, run_id, *, choices=None):
    return _atomic_edit(document, choices, _delete_named_run, name, run_id)


def clear_named_runs(document, name, *, choices=None):
    return _atomic_edit(document, choices, _clear_named_runs, name)


def delete_system(document, system_id, *, choices=None):
    return _atomic_edit(document, choices, _delete_system, system_id)


def _change_group(document, old_id, new_id=None, *, duplicate=False, choices=None):
    """Rename/duplicate/delete while migrating only actual shared-group callers."""
    from ..config import merge
    from ..planning import authoring_catalog
    before = authoring_catalog(merge(DEFAULTS, document), include_inactive=True)
    groups = document.setdefault("selections", {}).setdefault("groups", {})
    if old_id not in groups:
        raise ValueError(f"Unknown group {old_id!r}")
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
        # Use the catalog's ownership to map changed figures, rather than
        # rewriting plot-ID fragments with a second identity registry.
        # Compare the registered profile identities before and after the
        # group edit. Registry ordering identifies corresponding variants
        # without encoding analysis-specific name patterns in Studio.
        before_profiles = {(item["analysis"], item["run_id"]): []
                           for item in before["profiles"]}
        after_profiles = {(item["analysis"], item["run_id"]): []
                          for item in after["profiles"]}
        for item in before["profiles"]:
            before_profiles[(item["analysis"], item["run_id"])].append(item["profile"])
        for item in after["profiles"]:
            after_profiles[(item["analysis"], item["run_id"])].append(item["profile"])
        profile_mapping = {}
        for identity, original_profiles in before_profiles.items():
            updated_profiles = after_profiles.get(identity, [])
            removed_profiles = [profile for profile in original_profiles
                                if profile not in updated_profiles]
            added_profiles = [profile for profile in updated_profiles
                              if profile not in original_profiles]
            if new_id is not None and len(removed_profiles) != len(added_profiles):
                raise ValueError("Cannot reconcile group-owned profile identities for "
                                 + repr(identity))
            profile_mapping.update(zip(
                removed_profiles,
                added_profiles if new_id is not None else [None] * len(removed_profiles),
            ))
        if "profiles" in selection:
            selection["profiles"] = [
                dict(item, profile=profile_mapping[item["profile"]])
                if item["profile"] in profile_mapping else item
                for item in selection["profiles"]
                if item["profile"] not in profile_mapping
                or profile_mapping[item["profile"]] is not None
            ]
        from ..reporting.plot_ids import safe_plot_stem
        old_scalars = {
            figure["id"]: (item["analysis"], figure["id"].rsplit(":", 1)[-1])
            for item in before["scalars"] if item["group"] == old_id
            for figure in item["figures"]
        }
        new_scalars = {
            (item["analysis"], figure["id"].rsplit(":", 1)[-1]): figure["id"]
            for item in after["scalars"] if item["group"] == new_id
            for figure in item["figures"]
        }
        known = {
            figure["id"] for figures in before["figures"].values()
            for figure in figures
        }
        after_ids = {
            figure["id"] for figures in after["figures"].values()
            for figure in figures
        }
        mapping = {
            identifier: (new_scalars.get(identity) if new_id is not None else None)
            for identifier, identity in old_scalars.items()
        }
        for old_profile, new_profile in profile_mapping.items():
            old_stem = safe_plot_stem(old_profile)
            new_stem = safe_plot_stem(new_profile) if new_profile else None
            prefix = "profiles:" + old_stem + ":"
            for identifier in known:
                if not identifier.startswith(prefix):
                    continue
                replacement = (
                    "profiles:" + new_stem + identifier[len(prefix) - 1:]
                    if new_stem else None
                )
                mapping[identifier] = (
                    replacement if replacement in after_ids else None
                )

        def migrate(ids):
            return [mapping.get(identifier, identifier) for identifier in ids
                    if identifier not in mapping or mapping[identifier] is not None]
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
        # Exact-ID export must not silently carry unowned imported figures forward.
        original = ((raw.get("plots", {}).get("include") or [])
                    if raw.get("plots", {}).get("enabled", False) else [])
        managed_ids = {fig["id"] for figures in catalog["figures"].values() for fig in figures}
        unknown = [identifier for identifier in original if identifier not in managed_ids]
        if unknown:
            raise ValueError("Imported figures lack authoring ownership contracts: " + ", ".join(unknown) +
                             ". Leave the imported selection unchanged or explicitly remove these IDs before exact export.")
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


def planning_snapshot(document, base):
    """Resolve one consistent draft for both plan and active figure catalog."""
    from ..planning import authoring_catalog
    cfg = resolve_config(document, base=Path(base))
    return (build_plan(cfg, Path(base) / "analysis.yaml"),
            authoring_catalog(cfg)["figures"])


def topology_preflight(document, base):
    """Read topology only through production preflight; never open trajectories."""
    from ..core.validation import validate_selections
    cfg = resolve_config(document, base=Path(base))
    return validate_selections(cfg, Path(base))


def inspect_results_directory(document, base, results_directory):
    """Read-only saved-output inventory delegated to the production pipeline."""
    from ..pipeline.runner import inspect_results_directory as inspect
    return inspect(document, base, results_directory)


def results_coverage(inventory):
    """Group exact analysis/run/condition/replica identities without collapsing runs."""
    grouped = {}
    for entry in inventory:
        key = (entry["analysis"], entry["run_id"], entry["condition"])
        grouped.setdefault(key, {})[entry["replica"]] = entry["status"]
    return [
        {"analysis": analysis, "run_id": run_id, "condition": condition,
         "replicas": dict(sorted(replicas.items())),
         "verified": sum(v == "verified_existing" for v in replicas.values()),
         "total": len(replicas),
         "status": ("verified_existing" if all(
             v == "verified_existing" for v in replicas.values()) else
             "partial" if any(v == "verified_existing" for v in replicas.values())
             else "unavailable")}
        for (analysis, run_id, condition), replicas in sorted(grouped.items())
    ]


def planned_statistics_sources(plan):
    """Known analysis-to-family dependencies; never imply saved-data success."""
    from ..planning import ANALYSIS_METADATA
    sources = {family["family"]: [] for family in plan["statistics_families"]}
    for analysis in plan["analyses"]:
        for family in ANALYSIS_METADATA.get(analysis["analysis"], ("", (), ()))[1]:
            sources.setdefault(family, []).append(
                f"{analysis['analysis']}:{analysis['run_id']}")
    return sources


def statistics_source_status(plan, inventory=None):
    """Expose planned producers and verified saved coverage without conflation."""
    from ..planning import ANALYSIS_METADATA
    scheduled = planned_statistics_sources(plan)
    coverage = results_coverage(inventory or [])
    return {
        family["family"]: {
            "status": "scheduled" if scheduled.get(family["family"]) else
                      "verified_existing" if (
                          any(family["family"] in ANALYSIS_METADATA.get(
                              row["analysis"], ("", (), ()))[1] for row in coverage)
                          and all(row["status"] == "verified_existing"
                                  for row in coverage if family["family"] in
                                  ANALYSIS_METADATA.get(row["analysis"], ("", (), ()))[1])
                      ) else "partial_existing" if any(
                          family["family"] in ANALYSIS_METADATA.get(
                              row["analysis"], ("", (), ()))[1] for row in coverage
                      ) else "unavailable",
            "scheduled_sources": scheduled.get(family["family"], []),
            "existing_coverage": [
                row for row in coverage
                if family["family"] in ANALYSIS_METADATA.get(
                    row["analysis"], ("", (), ()))[1]
            ],
        }
        for family in plan["statistics_families"]
    }


def unavailable_exact_statistics(document, plan):
    """Reject explicit requested identities with no scheduled input producer.

    Unrestricted families are not errors: the backend may skip unsupported
    metrics, and existing results require separate verified reuse contracts.
    """
    if not document.get("statistics", {}).get("enabled", False):
        return []
    status = statistics_source_status(plan)
    return [
        family for family, info in status.items()
        if info["status"] == "unavailable"
        and statistics_family_mode(document, family) == "exact"
    ]


def statistics_family_mode(document, family):
    """Return 'all', 'none', or 'exact' without enumerating dynamic identities.

    None and an omitted family both mean all eligible backend identities.
    Explicit empty lists mean none; nonempty lists are exact allowlists.
    """
    from ..statistics_selection import studio_metric_selection_defaults
    if family not in studio_metric_selection_defaults():
        raise ValueError(f"Unknown statistics family: {family}")
    selectors = document.get("statistics", {}).get("metric_selection")
    if selectors is None or family not in selectors:
        return "all"
    return "none" if not selectors[family] else "exact"


def set_statistics_family_mode(document, family, mode):
    """Update a whole family without expanding unverified dynamic variables.

    A family-wide 'all' is represented by an omitted family key, preserving
    backend semantics for present and future dynamic identities.
    """
    if mode not in {"all", "none"}:
        raise ValueError("Statistics family mode must be 'all' or 'none'")
    statistics_family_mode(document, family)
    stats = document.setdefault("statistics", {})
    selectors = stats.get("metric_selection")
    if selectors is None:
        if mode == "all":
            return
        selectors = {}
        stats["metric_selection"] = selectors
    if mode == "all":
        selectors.pop(family, None)
    else:
        selectors[family] = []


def set_statistics_exact_family(document, family, entries, *, catalog_complete=False):
    """Require proven complete identity coverage when narrowing legacy-all.

    Explicit exact selections may be authored directly from an already exact
    family, but converting an unrestricted family would otherwise silently
    exclude dynamic identities absent from the Studio planning catalog.
    """
    from ..statistics_selection import validate_metric_selection
    before = statistics_family_mode(document, family)
    if before == "all" and not catalog_complete:
        raise ValueError(
            f"Cannot convert unrestricted {family} statistics to exact IDs: "
            "the planning catalog does not establish complete dynamic identity "
            "coverage. Keep Select All or supply a verified complete catalog."
        )
    candidate = deepcopy(document.get("statistics", {}))
    selectors = candidate.get("metric_selection")
    if selectors is None:
        selectors = {}
    else:
        selectors = deepcopy(selectors)
    selectors[family] = deepcopy(entries)
    candidate["metric_selection"] = selectors
    validate_metric_selection(candidate)
    document.setdefault("statistics", {})["metric_selection"] = selectors


def statistics_selector_fields(family):
    """Guided table columns from backend identity and estimator contracts."""
    from ..statistics_selection import SIMPLE_FAMILY_KEYS, ESTIMATOR_FLAGS
    keys = (("kind", "interaction", "metric") if family == "interactions"
            else SIMPLE_FAMILY_KEYS[family])
    return [*keys, "inference", *sorted(ESTIMATOR_FLAGS[family])]


def _delete_named_run(document, name, run_id, *, choices=None):
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
    previous = deepcopy(document)
    del runs[run_id]
    _prune_removed_catalog_outputs(previous, document, choices)


def _clear_named_runs(document, name, *, choices=None):
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
    previous = deepcopy(document)
    spec.pop("runs", None)
    _prune_removed_catalog_outputs(previous, document, choices)


def _prune_removed_catalog_outputs(previous, document, choices):
    """Discard only identities actually removed from the shared catalog."""
    from ..config import merge
    from ..planning import authoring_catalog

    old = authoring_catalog(merge(DEFAULTS, previous), include_inactive=True)
    new = authoring_catalog(merge(DEFAULTS, document), include_inactive=True)
    valid_profiles = {item["profile"] for item in new["profiles"]}
    removed_profiles = {
        item["profile"] for item in old["profiles"]
        if item["profile"] not in valid_profiles
    }
    selectors = document.get("statistics", {}).get("metric_selection")
    if isinstance(selectors, dict) and "profiles" in selectors:
        selectors["profiles"] = [
            item for item in selectors["profiles"]
            if item["profile"] not in removed_profiles
        ]
        # Exact selector entries generate statistical figure candidates.
        # Rebuild after pruning stale selectors, or the retired figure
        # incorrectly appears active and survives the catalog difference.
        new = authoring_catalog(merge(DEFAULTS, document),
                                include_inactive=True)
    old_ids = {fig["id"] for figures in old["figures"].values()
               for fig in figures}
    new_ids = {fig["id"] for figures in new["figures"].values()
               for fig in figures}
    removed_figures = old_ids - new_ids
    plots = document.get("plots", {})
    if plots.get("include") is not None:
        plots["include"] = [
            item for item in plots["include"] if item not in removed_figures
        ]
    if choices is not None:
        for owner in choices:
            choices[owner] = [
                item for item in choices[owner]
                if item not in removed_figures
            ]


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
    commit_draft(document, candidate)


def add_group(document, group_id):
    """Validate a detached group addition before modifying the draft."""
    from ..config import analysis_groups
    candidate = deepcopy(document)
    groups = candidate.setdefault("selections", {}).setdefault("groups", {})
    if not group_id or group_id in groups:
        raise ValueError("Enter a unique filesystem-safe group ID")
    groups[group_id] = {"label": group_id, "selection": "all"}
    analysis_groups(candidate)
    commit_draft(document, candidate)


def add_system(document, system_id):
    """Validate a detached condition addition before modifying the draft."""
    from ..config import _validate_path_component
    if not system_id or system_id in document.get("systems", {}):
        raise ValueError("Enter a unique filesystem-safe condition ID")
    _validate_path_component(system_id, "system name")
    candidate = deepcopy(document)
    candidate.setdefault("systems", {})[system_id] = {
        "topology": "topology.tpr",
        "replicates": [{"label": "rep1", "trajectory": "rep1.xtc"}],
    }
    commit_draft(document, candidate)


def add_replica(document, system_id):
    """Add a unique, validated replica to a detached draft."""
    from ..config import _validate_path_component
    candidate = deepcopy(document)
    if system_id not in candidate.get("systems", {}):
        raise ValueError(f"Unknown condition {system_id!r}")
    replicas = candidate["systems"][system_id]["replicates"]
    used = {(rep.get("label") or f"rep{i}") if isinstance(rep, dict) else f"rep{i}"
            for i, rep in enumerate(replicas, 1)}
    number = 1
    while f"rep{number}" in used:
        number += 1
    label = f"rep{number}"
    _validate_path_component(label, "replica label")
    replicas.append({"label": label, "trajectory": f"{label}.xtc"})
    commit_draft(document, candidate)


def remove_replica(document, system_id, index):
    """Reject invalid removals without changing the draft or figure choices."""
    candidate = deepcopy(document)
    if system_id not in candidate.get("systems", {}):
        raise ValueError(f"Unknown condition {system_id!r}")
    replicas = candidate["systems"][system_id]["replicates"]
    if len(replicas) <= 1:
        raise ValueError("Keep at least one replica per condition")
    if not 0 <= index < len(replicas):
        raise ValueError(f"Unknown replica index {index}")
    replicas.pop(index)
    commit_draft(document, candidate)


def _delete_system(document, system_id, *, choices=None):
    """Remove condition-owned figures without changing shared figure choices."""
    systems = document["systems"]
    if system_id not in systems:
        raise ValueError(f"Unknown condition {system_id!r}")
    if len(systems) == 1:
        raise ValueError("Keep at least one condition")
    previous = deepcopy(document)
    del systems[system_id]
    _prune_removed_catalog_outputs(previous, document, choices)
