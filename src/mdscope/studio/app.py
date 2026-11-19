"""Optional localhost configuration wizard. Run via mdscope-studio."""
from copy import deepcopy
from pathlib import Path

import streamlit as st
import yaml

from mdscope.config import DEFAULTS, analysis_groups, resolve_config
from mdscope.studio.model import (
    analysis_options, change_group, import_document, new_document, review,
    topology_preflight,
)
from mdscope.statistics_selection import ESTIMATOR_FLAGS, studio_metric_selection_defaults


def yaml_editor(label, value, key, *, height=140, help=None):
    text = st.text_area(label, yaml.safe_dump(value, sort_keys=False), key=key, height=height, help=help)
    return yaml.safe_load(text)


def analysis_settings(name, spec, document):
    """Edit backend defaults with typed widgets; optional fields stay YAML."""
    defaults = DEFAULTS["analyses"].get(name, {})
    from mdscope.config_schema import option_help
    result = {}
    for key, default in defaults.items():
        if key == "enabled" or (key == "groups" and name in {"rmsd", "radius_of_gyration"}):
            continue
        value = spec.get(key, default)
        label = key.replace("_", " ")
        widget_key = f"edit:setting:{name}:{key}"
        help_text = option_help(name, key)
        if type(default) is bool and type(value) is bool:
            edited = st.checkbox(label, value=value, key=widget_key, help=help_text)
        elif type(default) in {int, float} and type(value) in {int, float}:
            edited = st.number_input(label, value=value, key=widget_key, help=help_text)
        elif isinstance(default, str) and isinstance(value, str):
            edited = st.text_input(label, value=value, key=widget_key, help=help_text)
        else:
            edited = yaml_editor(label + " (YAML; null = backend default)", value, widget_key, height=80, help=help_text)
        if key in spec or edited != default:
            result[key] = edited
    from mdscope.config import MULTI_RUN_ANALYSES
    excluded = {"enabled"} | ({"groups"} if name in {"rmsd", "radius_of_gyration"} else set())
    if name in MULTI_RUN_ANALYSES:
        # A run action reruns immediately, before this editor returns. Preserve
        # preceding widget edits in the detached draft before that commit.
        spec.update(result)
        excluded.add("runs")
        use_runs = st.checkbox("Configure named runs", value="runs" in spec, key=f"edit:runs:{name}:enabled")
        if use_runs:
            runs = spec.setdefault("runs", {})
            for run_id in list(runs):
                with st.expander(f"Run · {run_id}", expanded=True):
                    current = runs[run_id]
                    current = {"enabled": current} if isinstance(current, bool) else current
                    run_options = yaml_editor("Run options (inherits analysis settings)", current,
                                              f"edit:run:{name}:{run_id}")
                    if not isinstance(run_options, dict):
                        raise ValueError("Run options must be a mapping")
                    runs[run_id] = run_options
                    if st.button("Delete run", key=f"edit:run:{name}:{run_id}:delete"):
                        from mdscope.studio.model import delete_named_run
                        delete_named_run(document, name, run_id, choices=st.session_state.choices)
                        restart_editor(document)
            new_id = st.text_input("New run ID", key=f"edit:runs:{name}:new")
            if st.button("Add run", key=f"edit:runs:{name}:add"):
                from mdscope.studio.model import add_named_run
                add_named_run(document, name, new_id)
                restart_editor(document)
            result["runs"] = runs
        else:
            if "runs" in spec:
                from mdscope.studio.model import clear_named_runs
                clear_named_runs(document, name, choices=st.session_state.choices)
    optional = {key: value for key, value in spec.items()
                if key not in defaults and key not in excluded}
    edited = yaml_editor("Additional supported options (YAML)", optional,
                         f"edit:analysis:{name}") or {}
    if not isinstance(edited, dict) or set(edited) & (set(defaults) | excluded):
        raise ValueError("Use the controls above for default options, enabled and groups")
    return {**result, **edited}


def restart_editor(document):
    from mdscope.studio.model import commit_draft
    commit_draft(st.session_state.document, document)
    # Rebuild widgets from the committed draft and explicit figure choices.
    # Deleted/recreated IDs must not inherit stale widget values.
    for key in list(st.session_state):
        if key.startswith(("edit:", "figure:", "stats:", "enable:")):
            del st.session_state[key]
    st.rerun()


@st.cache_data(show_spinner=False, max_entries=32)
def preview_document(text, base):
    """Cache pure schema-only planning by draft contents; edits invalidate it."""
    return review(yaml.safe_load(text), base)


def show_outputs(plan, name, document, base):
    st.caption("ILLUSTRATIVE — NO MD RESULTS. Placeholders describe future saved data.")
    for analysis in plan["analyses"]:
        if analysis["analysis"] == name:
            st.caption(analysis["purpose"] + " · coordinate source: " + analysis["coordinate_source"])
            if analysis.get("required_topology"):
                st.caption("Prerequisites: " + "; ".join(analysis["required_topology"]))
    for variable in plan["observable_contracts"]:
        if variable["analysis"] != name:
            continue
        st.write(f"{variable['metric']} · {variable['label']}")
        st.caption(variable["row_grain"] + "; columns are a subset of a shared table.")
        rows = [{column: (f"t{i} (ps)" if column == "time_ps" else
                          variable["group"] if column == "group" else f"example ({variable['units'].get(column, 'unknown')})")
                 for column in variable["columns"]} for i in range(2)]
        st.table(rows)
        st.json({key: variable[key] for key in (
            "measurement_selection", "group_selection", "selection_semantics",
            "alignment_selection", "coordinate_source", "units")})
        for system, spec in document["systems"].items():
            for i, rep in enumerate(spec["replicates"], 1):
                label = rep.get("label", f"rep{i}") if isinstance(rep, dict) else f"rep{i}"
                st.code(str(Path(base) / plan["project"]["output"] / system / label /
                            variable["path_within_replica"]))
    outputs = [item for item in (plan["registered_replica_outputs"] + plan["structural_profile_outputs"]) if item["analysis"] == name]
    if outputs:
        for output in outputs:
            st.write(f"{output['run_id']} · {output['role']} · {output['format']}")
            if output.get("statistics_requirement"):
                st.info(output["statistics_requirement"])
            st.caption(output.get("row_grain") or "Row grain / array dimensions depend on the producer contract; unresolved fields remain unknown.")
            if output.get("columns"):
                st.table([{column: "illustrative placeholder" for column in output["columns"]}
                          for _ in range(2)])
                st.json({"units": output.get("units"), "columns_status": output["columns_status"]})
            elif output.get("array_keys") or output.get("shape_axes"):
                st.json({key: output.get(key) for key in ("array_keys", "shape_axes", "dtype", "units")})
            else:
                st.caption(output.get("unresolved_reason") or "Metadata schema depends on selections/topology.")
            for replica in plan["replicates"]:
                st.code(str(Path(base) / plan["project"]["output"] / replica["system"] /
                            replica["id"] / output["path_within_replica"]))
    elif not any(item["analysis"] == name for item in plan["observable_contracts"]):
        st.info("Output schema depends on selections/topology or is not yet registered. See review plan.")


def main():
    st.set_page_config(page_title="MDscope Studio", layout="wide")
    st.title("MDscope Studio")
    st.caption("Configure analyses → choose optional figures/statistics → export YAML. No trajectory jobs run here.")
    state = st.session_state
    state.setdefault("document", new_document())
    state.setdefault("choices", {})
    base = st.text_input("Configuration directory (relative paths resolve here)", str(Path.cwd()))
    results_dir = st.text_input(
        "Existing MDscope results directory (optional)",
        value=str(Path(base) / state.document.get("project", {}).get(
            "output", DEFAULTS["project"]["output"])),
        key="studio:results_dir",
        help="Inspect saved registered analysis products without reading trajectories. "
             "A file inventory alone cannot verify scientific provenance.")
    if st.button("Inspect existing results", key="studio:inspect_results"):
        try:
            from mdscope.studio.model import inspect_results_directory
            state["studio:results_inventory"] = inspect_results_directory(
                state.document, base, results_dir)
            state["studio:results_inventory_path"] = results_dir
            state["studio:results_inventory_config"] = yaml.safe_dump(state.document, sort_keys=True)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            state.pop("studio:results_inventory", None)
            st.error(str(exc))
    if "studio:results_inventory" in state:
        inventory = state["studio:results_inventory"]
        if (state.get("studio:results_inventory_path") != results_dir or
                state.get("studio:results_inventory_config") !=
                yaml.safe_dump(state.document, sort_keys=True)):
            st.warning("Results inventory is stale after configuration/path changes; inspect again.")
        st.caption("Verification requires matching saved configuration, completed "
                   "analysis stage marker, replica summary and production artifact "
                   "schema. Unregistered outputs remain unavailable.")
        from mdscope.studio.model import results_coverage
        coverage = results_coverage(inventory)
        if coverage:
            st.write("**Verified coverage by analysis, run and condition**")
            st.dataframe(coverage, hide_index=True)
        if inventory:
            st.dataframe(inventory, hide_index=True)
        else:
            st.info("No registered outputs are expected for the currently enabled analyses.")

    uploaded = st.file_uploader("Import existing YAML", type=["yaml", "yml"])
    if uploaded is not None and st.button("Load YAML"):
        try:
            document = import_document(uploaded.getvalue().decode("utf-8"))
            # Recreate widget state as well as the document on import.
            for key in list(state):
                if key.startswith(("edit:", "enable:", "figure:", "stats:")):
                    del state[key]
            state.document = document
            state.choices = {}
            state.pop("preflight", None)
            st.rerun()
        except (ValueError, yaml.YAMLError, UnicodeError) as exc:
            st.error(str(exc))
    doc = deepcopy(state.document)
    choices_before = deepcopy(state.choices)
    plots = doc.get("plots", {})
    # One hierarchy for imported family switches and exact selections.
    # Imported configurations are not converted until a checkbox is changed.
    st.caption("Figures · select all, whole categories, or individual outputs.")
    from mdscope.planning import authoring_catalog
    from mdscope.config import merge, figure_enabled
    catalog = authoring_catalog(merge(DEFAULTS, doc), include_inactive=True)["figures"]
    figures_by_family = {}
    for owner, figures in catalog.items():
        for fig in figures:
            family = fig["id"].split(":", 1)[0]
            figures_by_family.setdefault(family, []).append((owner, fig))
    all_ids = list(dict.fromkeys(
        fig["id"] for entries in figures_by_family.values() for _, fig in entries))
    exact_import = plots.get("include") is not None
    initial = set(plots.get("include") or []) if exact_import else {
        identifier for identifier in all_ids
        if figure_enabled(resolve_config(doc), identifier)
    }
    initial.update(identifier for values in state.choices.values() for identifier in values
                   if exact_import)
    def selected(identifier):
        return state.get(f"figure:tree:{identifier}", identifier in initial)
    def sync_parents():
        for family, entries in figures_by_family.items():
            ids = {fig["id"] for _, fig in entries}
            state[f"figure:tree:family:{family}"] = all(selected(i) for i in ids)
        state["figure:tree:all"] = bool(all_ids) and all(selected(i) for i in all_ids)
    def set_figures(ids, checked):
        for identifier in ids:
            state[f"figure:tree:{identifier}"] = checked
        state["figure:tree:edited"] = True
        sync_parents()
    def changed_figure():
        state["figure:tree:edited"] = True
        sync_parents()
    if all_ids:
        st.checkbox("Select All", value=all(selected(i) for i in all_ids),
                    key="figure:tree:all",
                    on_change=lambda: set_figures(all_ids, state["figure:tree:all"]))
    for family, entries in figures_by_family.items():
        ids = list(dict.fromkeys(fig["id"] for _, fig in entries))
        with st.expander(f"{family.replace('_', ' ').title()} ({len(ids)} figures)"):
            st.checkbox(f"Select all {family.replace('_', ' ')}",
                        value=all(selected(i) for i in ids),
                        key=f"figure:tree:family:{family}",
                        on_change=lambda ids=ids, family=family: set_figures(
                            ids, state[f"figure:tree:family:{family}"]))
            for owner, fig in entries:
                identifier = fig["id"]
                st.checkbox(fig.get("label", identifier),
                            value=identifier in initial, key=f"figure:tree:{identifier}",
                            help=identifier, on_change=changed_figure)
                st.caption(f"{owner} · {fig['destination']}")
                if "studio:results_inventory" in state and owner != "statistics":
                    matches = [item for item in state["studio:results_inventory"]
                               if item["analysis"] == owner]
                    if matches:
                        from mdscope.studio.model import results_coverage
                        for item in results_coverage(matches):
                            st.caption(
                                f"Saved {item['run_id']} / {item['condition']}: "
                                f"{item['verified']}/{item['total']} replicas verified "
                                f"({item['status']}). Unverified inputs require scheduling.")
    exact_figures = exact_import or bool(state.get("figure:tree:edited"))
    if exact_import and not state.get("figure:tree:edited") and not state.choices:
        state.choices = {
            owner: [fig["id"] for fig in figures if fig["id"] in initial]
            for owner, figures in catalog.items()
        }
    if state.get("figure:tree:edited"):
        chosen = {identifier for identifier in all_ids if selected(identifier)}
        state.choices = {
            owner: [fig["id"] for fig in figures if fig["id"] in chosen]
            for owner, figures in catalog.items()
        }
        # Exact IDs take precedence over imported family switches on export.
        doc.setdefault("plots", {})["include"] = list(chosen)
        doc["plots"]["enabled"] = bool(chosen)
        if not exact_import:
            st.warning("Editing individual figures converts category switches to exact IDs. "
                       "Figures whose identities depend on unavailable topology/data "
                       "may not be listed; review the exported figure list.")
    tabs = st.tabs(["1 · Project & systems", "2 · Selections & alignment",
                    "3 · Analyses", "4 · Statistics", "5 · Review & export"])
    errors = []
    with tabs[0]:
        try:
            advanced_project = st.checkbox("Edit project setup as advanced YAML", key="edit:project:advanced")
            if state.get("project_editor_mode") != advanced_project:
                for key in list(state):
                    if key != "edit:project:advanced" and (key in {"edit:project", "edit:systems", "edit:trajectory", "edit:time"}
                            or key.startswith(("edit:project:", "edit:system:", "edit:rep:", "edit:trajectory:", "edit:time:"))):
                        del state[key]
                state.project_editor_mode = advanced_project
            if advanced_project:
                for section in ("project", "systems", "trajectory", "time"):
                    edited = yaml_editor(section, doc.get(section, DEFAULTS.get(section, {})), f"edit:{section}")
                    if not isinstance(edited, dict):
                        raise ValueError(f"{section} must be a mapping")
                    doc[section] = edited
            else:
                from mdscope.studio.project_editor import project_settings
                project_settings(doc, restart_editor)
            st.caption("Use conditions as system IDs; replicate entries accept label, trajectory and optional msd_trajectory. Time units: ps.")
            if st.button("Check input paths (read-only)"):
                resolve_config(doc, base=Path(base), check_paths=True)
                st.success("Configured input paths exist; trajectory contents were not checked.")
        except (ValueError, TypeError, OSError, yaml.YAMLError) as exc:
            errors.append(str(exc)); st.error(str(exc))
    with tabs[1]:
        st.write("Free-form MDAnalysis atom selections")
        st.caption("Examples: protein and segid PROA; resid 10:60 and backbone; resname POPC and name P. Expressions are preserved exactly.")
        try:
            selections = doc.setdefault("selections", {})
            common = {key: value for key, value in selections.items() if key != "groups"}
            common = yaml_editor("Global selections (measurement features and fitting atoms are independent)",
                                 common, "edit:selections")
            if not isinstance(common, dict):
                raise ValueError("Selections must be a mapping")
            if "groups" in common:
                raise ValueError("Use the named group editor for groups")
            # Replace global keys so YAML deletions survive export. Named groups
            # are edited independently and retain their existing identity.
            groups = selections.get("groups", {})
            selections.clear()
            selections.update(common)
            selections["groups"] = groups
            for group_id in list(groups):
                if isinstance(groups[group_id], str):
                    groups[group_id] = {"selection": groups[group_id]}
                with st.expander(group_id):
                    spec = groups[group_id]
                    for field, label in (("label", "Display label"), ("selection", "Measurement group selection"),
                                         ("alignment_selection", "Independent alignment selection (optional)")):
                        value = st.text_input(label, spec.get(field) or "", key=f"edit:group:{group_id}:{field}")
                        if value:
                            spec[field] = value
                        else:
                            spec.pop(field, None)
                    target = st.text_input("New group ID for rename/duplicate", key=f"edit:group:{group_id}:id")
                    for action in ("Rename", "Duplicate", "Delete"):
                        if st.button(action, key=f"edit:group:{group_id}:{action}"):
                            if action != "Delete" and not target:
                                raise ValueError("Enter a new group ID")
                            change_group(doc, group_id, None if action == "Delete" else target,
                                         duplicate=action == "Duplicate", choices=state.choices)
                            restart_editor(doc)
            new_id = st.text_input("Group ID", key="edit:new_group")
            if st.button("Add group"):
                from mdscope.studio.model import add_group
                add_group(doc, new_id)
                restart_editor(doc)
            if st.button("Preview selections from topology (no trajectories)", key="edit:groups:preflight"):
                # An explicit topology-only read, using the same validator as export.
                counts = topology_preflight(doc, base)
                st.success("Topology selections parsed. Trajectories, coordinates and chemistry beyond selection checks remain unverified.")
                for system_id, result in counts.items():
                    st.write(f"Condition · {system_id}")
                    for group_id, info in result.get("groups", {}).items():
                        matched = info["measurement"]
                        st.caption(f"{group_id}: {matched['atoms']} atoms; {matched['residues']} residues")
                    st.json(result)
            else:
                st.info("Selections not checked against topology. Fitted coordinates apply only to supported structural measurements; membrane/diffusion use source coordinates.")
        except (ValueError, TypeError, OSError, yaml.YAMLError) as exc:
            errors.append(str(exc)); st.error(str(exc))
    with tabs[2]:
        for name, schema in analysis_options().items():
            spec = doc.setdefault("analyses", {}).setdefault(name, {"enabled": DEFAULTS["analyses"].get(name, {}).get("enabled", False)})
            if isinstance(spec, bool):
                spec = doc["analyses"][name] = {"enabled": spec}
            enabled = st.checkbox(f"Enable {name}", value=spec.get("enabled", False), key=f"enable:{name}")
            spec["enabled"] = enabled
            if not enabled:
                continue
            with st.expander(name, expanded=True):
                try:
                    if name in {"rmsd", "radius_of_gyration"}:
                        default = spec.get("groups", list(doc["selections"].get("groups", {})))
                        spec["groups"] = st.multiselect("Named measurement groups", list(doc["selections"].get("groups", {})),
                                                         default=default, key=f"edit:groups:{name}")
                        st.caption("Global measurement is also saved; this list filters grouped rows only.")
                    edited = analysis_settings(name, spec, doc)
                    doc["analyses"][name] = {"enabled": True, **edited,
                                              **({"groups": spec["groups"]} if name in {"rmsd", "radius_of_gyration"} and "groups" in spec else {})}
                    st.caption("Supported keys: " + ", ".join(key for key in schema if key != "enabled"))
                    st.caption("Units: _A = Å; _ps = ps; _deg = degrees. Named runs are supported only where runs is listed.")
                    _, plan = preview_document(yaml.safe_dump(doc, sort_keys=False), base)
                    show_outputs(plan, name, doc, base)
                    from mdscope.planning import authoring_catalog
                    figures = authoring_catalog(resolve_config(doc))["figures"].get(name, [])

                    if not figures:
                        st.caption("No analysis-only figure IDs verified for this card yet; no figures selected automatically.")
                except (ValueError, TypeError, KeyError, yaml.YAMLError) as exc:
                    errors.append(str(exc)); st.error(str(exc))
    with tabs[3]:
        stats = doc.setdefault("statistics", {})
        stats["enabled"] = st.checkbox("Enable statistics", value=stats.get("enabled", DEFAULTS["statistics"]["enabled"]), key="stats:master")
        st.caption("Independent units: simulation replicates. Eligibility requires compatible sampling/identities; CI availability is determined by the backend, never by frame count alone.")
        if stats["enabled"]:
            try:
                _, plan = preview_document(yaml.safe_dump(doc, sort_keys=False), base)
                from mdscope.studio.model import (
                    statistics_family_mode, set_statistics_family_mode,
                    set_statistics_exact_family, statistics_selector_fields,
                )
                from mdscope.statistics_selection import validate_metric_selection
                import pandas as pd
                families = tuple(studio_metric_selection_defaults())
                from mdscope.studio.model import statistics_source_status
                inventory_fresh = (
                    state.get("studio:results_inventory_path") == results_dir
                    and state.get("studio:results_inventory_config") ==
                    yaml.safe_dump(state.document, sort_keys=True))
                source_status = statistics_source_status(
                    plan, state.get("studio:results_inventory") if inventory_fresh else None)
                st.caption("Select all statistics, entire families, or individual known variables. "
                           "All includes dynamic identities that cannot be enumerated before execution.")
                def set_family(family, checked):
                    set_statistics_family_mode(state.document, family, "all" if checked else "none")
                    # The widget values are reset when the family changes mode.
                    state.pop(f"stats:table:{family}", None)
                    state.pop(f"stats:tree:convert:{family}", None)
                    state["stats:tree:all"] = all(
                        statistics_family_mode(state.document, item) == "all" for item in families)
                def set_all_families():
                    checked = state["stats:tree:all"]
                    for family in families:
                        set_family(family, checked)
                        state[f"stats:tree:family:{family}"] = checked
                st.checkbox("Select All Statistics",
                            value=all(statistics_family_mode(doc, family) == "all"
                                      for family in families),
                            key="stats:tree:all", on_change=set_all_families)
                for family in families:
                    mode = statistics_family_mode(doc, family)
                    availability = source_status.get(family, {"status": "unavailable", "scheduled_sources": []})
                    with st.expander(f"{family.replace('_', ' ').title()} · {mode.title()} · {availability['status']}"):
                        if availability["scheduled_sources"]:
                            st.caption("Scheduled analysis inputs: " +
                                       ", ".join(availability["scheduled_sources"]))
                        for source in availability.get("existing_coverage", []):
                            st.caption(
                                f"Saved {source['analysis']} / {source['run_id']} / "
                                f"{source['condition']}: {source['verified']}/"
                                f"{source['total']} replicas verified "
                                f"({source['status']})")
                        else:
                            st.warning("No scheduled producer for this family. Existing "
                                       "files must be verified before reuse; the results "
                                       "inventory does not establish that verification.")
                        st.checkbox(f"All {family.replace('_', ' ')} statistics",
                                    value=mode == "all",
                                    key=f"stats:tree:family:{family}",
                                    on_change=lambda family=family: set_family(
                                        family, state[f"stats:tree:family:{family}"]))
                        st.caption("All includes future and data-dependent identities."
                                   if mode == "all" else
                                   "Explicit allowlist; unchecked variables are excluded.")
                        if mode == "all":
                            st.info("This family is unrestricted. To customize individual "
                                    "identities, first deselect the family; the resulting "
                                    "empty allowlist can then be populated explicitly.")
                            continue
                        selectors = stats.get("metric_selection") or {}
                        existing = selectors.get(family, [])
                        if family in {"time_series", "profiles"}:
                            variables = (plan["observable_contracts"] if family == "time_series"
                                         else plan["profile_variables"])
                            known = []
                            updated = []
                            for var in variables:
                                selector = var["statistics_selector"]
                                identity = ((selector["metric"], selector.get("group"))
                                            if family == "time_series" else selector["profile"])
                                known.append(identity)
                                selected = next((item for item in existing if
                                    ((item["metric"], item.get("group"))
                                     if family == "time_series" else item["profile"]) == identity), None)
                                key = (f"stats:{var['metric']}:{var['group']}"
                                       if family == "time_series" else
                                       f"stats:profile:{var['profile']}")
                                if st.checkbox(var.get("label", str(identity)),
                                               value=selected is not None, key=key):
                                    entry = {**(selected or {}), **selector}
                                    for flag in sorted(ESTIMATOR_FLAGS[family]):
                                        entry[flag] = st.checkbox(
                                            flag, value=(selected or {}).get(
                                                flag, (selected or {}).get("inference", True)),
                                            key=f"{key}:{flag}")
                                    updated.append(entry)
                            updated.extend(item for item in existing if
                                ((item["metric"], item.get("group"))
                                 if family == "time_series" else item["profile"]) not in known)
                        else:
                            fields = statistics_selector_fields(family)
                            source = pd.DataFrame(existing, columns=fields)
                            columns = {field: st.column_config.CheckboxColumn(field, default=True)
                                       for field in fields
                                       if field in {"inference", "ci", "hypothesis_tests"}}
                            edited = st.data_editor(source, num_rows="dynamic", hide_index=True,
                                                    column_config=columns,
                                                    key=f"stats:table:{family}")
                            updated = [{key: value for key, value in row.items() if pd.notna(value)}
                                       for row in edited.to_dict("records")]
                        if updated != existing:
                            set_statistics_exact_family(doc, family, updated)
                validate_metric_selection(stats)
                st.markdown("**Interaction statistics · minimum occupancy**")
                cutoff = float(stats.get("interaction_min_occupancy", 0.0))
                if not 0.0 <= cutoff <= 1.0:
                    raise ValueError("statistics.interaction_min_occupancy must be between 0 and 1")
                selected_cutoff = st.slider(
                    "Retain interactions with occupancy in at least one replicate",
                    min_value=0.0, max_value=1.0, value=cutoff, step=0.01,
                    format="%.2f", key="stats:interaction_min_occupancy",
                    help="0 disables filtering. 0.20 retains an interaction if it "
                         "occupies at least 20% of sampled frames in any replicate. "
                         "This filters interaction identities before statistical inference.")
                if "interaction_min_occupancy" in stats or selected_cutoff != 0.0:
                    stats["interaction_min_occupancy"] = selected_cutoff
                controls = yaml_editor("CI / block / test settings", {key: val for key, val in stats.items()
                                                                     if key not in {"enabled", "metric_selection", "interaction_min_occupancy"}}, "stats:settings") or {}
                if set(controls) & {"enabled", "metric_selection", "interaction_min_occupancy"}:
                    raise ValueError("Use the statistics selectors above")
                stats.update(controls)
                st.caption("Dynamic metric identities outside the structural scalar catalog require exact backend IDs. No identity or output schema is invented.")
                with st.expander("Expected statistics outputs and eligibility"):
                    st.caption("ILLUSTRATIVE — NO MD RESULTS. Candidate files are conditional on eligible data; Base schemas below include only verified writer columns; optional provenance remains conditional.")
                    for family in plan["statistics_families"]:
                        schema = family.get("output_schema")
                        if schema:
                            st.write(f"{family['family']} · {schema['row_grain']}")
                            st.code(f"{family['output_directory']}{schema['filename']}")
                            st.table(schema["sample_rows"])
                            st.caption(schema["limitations"])
                        st.json({key: family[key] for key in ("family", "files", "metric_contract")
                                 if key in family})
                from mdscope.planning import authoring_catalog
                figures = authoring_catalog(resolve_config(doc))["figures"].get("statistics", [])

            except (ValueError, TypeError, KeyError, yaml.YAMLError) as exc:
                errors.append(str(exc)); st.error(str(exc))
    with tabs[4]:
        try:
            text, plan = review(doc, base, state.choices if exact_figures else None)
            st.warning("Schema-only validation: topology, chemistry, trajectory alignment and sampling eligibility remain unverified.")
            from mdscope.studio.model import unavailable_exact_statistics
            missing_stats = unavailable_exact_statistics(yaml.safe_load(text), plan)
            if missing_stats:
                st.error("Exact statistics requested without scheduled source analyses: " +
                         ", ".join(missing_stats) +
                         ". Enable a compatible analysis; existing files are not "
                         "yet verified as reusable.")

            require_topology = st.checkbox("Require topology preflight before export")
            signature = text + str(Path(base).resolve())
            if st.button("Validate topology selections (no XTC)"):
                state.pop("preflight", None)
                counts = topology_preflight(yaml.safe_load(text), base)
                state.preflight = (signature, counts)
            preflight = state.get("preflight")
            if preflight and preflight[0] == signature:
                st.success("Configured topology preflight passed for the current document; trajectory contents and unsupported chemistry checks remain unverified.")
                st.json(preflight[1])
                for group in analysis_groups(resolve_config(doc)):
                    identities = {info["groups"][group]["features" if info["groups"][group].get("feature_status") == "required" else "measurement"]["ordered_identity_sha256"]
                                  for info in preflight[1].values()}
                    if len(identities) > 1:
                        st.warning(f"{group}: selected atom identity/order differs across systems. Check comparability before inference.")
            st.json({"analyses": plan["analyses"], "figures": yaml.safe_load(text).get("plots"),
                     "statistics": yaml.safe_load(text).get("statistics")})
            st.code(text, language="yaml")
            valid = not errors and not missing_stats and (not require_topology or bool(preflight and preflight[0] == signature))
            st.download_button("Download analysis.yaml", text, file_name="analysis.yaml", mime="text/yaml", disabled=not valid)
            st.caption("No output files are overwritten. Save YAML in the configuration directory shown above.")
        except (ValueError, TypeError, KeyError, OSError, yaml.YAMLError) as exc:
            st.error(str(exc))

    # A rejected action must not commit unrelated partial edits or figure
    # selections from the same rerun. Keep widget input available for repair.
    if errors:
        state.choices = choices_before
    else:
        from mdscope.studio.model import commit_draft
        commit_draft(state.document, doc)


if __name__ == "__main__":
    main()
