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


def figure_picker(figures, key, existing):
    selected = []
    for fig in figures:
        if st.checkbox(fig.get("label", fig["id"]), value=fig["id"] in existing, key=f"{key}:{fig['id']}", help=fig["id"]):
            selected.append(fig["id"])
        st.caption(f"{fig['destination']} · {'Statistics required' if fig['requires_statistics'] else 'Descriptive'}")
    return selected


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
    exact_figures = st.checkbox("Use exact figure choices", value=(plots.get("include") is not None
                               or not plots.get("enabled", False)), key="edit:exact_figures")
    if not exact_figures:
        st.info("Imported legacy figure switches are preserved. Enable exact figure choices to replace them with an explicit allowlist.")
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
                    if exact_figures:
                        state.choices[name] = figure_picker(figures, f"figure:{name}", state.choices.get(name, (doc.get("plots", {}).get("include") or []) if doc.get("plots", {}).get("enabled") else []))
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
                exact_stats = st.checkbox("Select statistics by exact variable",
                                          value=stats.get("metric_selection") is not None,
                                          key="stats:exact")
                if not exact_stats:
                    # Null is the backend's legacy-all selection. Restore the
                    # prior exact allowlist if the user opts back in this session.
                    if stats.get("metric_selection") is not None:
                        state["stats:exact:saved"] = deepcopy(stats["metric_selection"])
                    stats["metric_selection"] = None
                    st.info("Legacy statistics mode includes all eligible variables. Enable exact variable selection to restore your previous choices or start a new allowlist.")
                else:
                    if stats.get("metric_selection") is None:
                        stats["metric_selection"] = deepcopy(
                            state.get("stats:exact:saved", studio_metric_selection_defaults()))
                    selectors = stats["metric_selection"]
                    entries = []
                    for var in plan["observable_contracts"]:
                        selector = var["statistics_selector"]
                        key = f"stats:{var['metric']}:{var['group']}"
                        selected = next((item for item in selectors.get("time_series", [])
                                         if item.get("group") == selector["group"] and item["metric"] == selector["metric"]), None)
                        if st.checkbox(f"Compute {var['metric']} statistics · {var['label']}", value=selected is not None, key=key):
                            entry = {**(selected or {}), **selector}
                            for flag in sorted(ESTIMATOR_FLAGS["time_series"]):
                                entry[flag] = st.checkbox(flag, value=(selected or {}).get(flag, (selected or {}).get("inference", True)), key=f"{key}:{flag}")
                            entries.append(entry)
                    known = {(var["metric"], var["group"]) for var in plan["observable_contracts"]}
                    selectors["time_series"] = entries + [item for item in selectors.get("time_series", [])
                                                             if (item["metric"], item.get("group")) not in known]
                    guided_profiles = []
                    known_profiles = {var["profile"] for var in plan["profile_variables"]}
                    for var in plan["profile_variables"]:
                        key = f"stats:profile:{var['profile']}"
                        selected = next((item for item in selectors.get("profiles", [])
                                         if item["profile"] == var["profile"]), None)
                        st.caption(f"{var['analysis']} · {var['run_id']} · {var['units']} · {var['prerequisite']}")
                        if st.checkbox(f"Compute {var['profile']} statistics", value=selected is not None, key=key):
                            entry = {**(selected or {}), **var["statistics_selector"]}
                            for flag in sorted(ESTIMATOR_FLAGS["profiles"]):
                                entry[flag] = st.checkbox(flag, value=(selected or {}).get(flag, (selected or {}).get("inference", True)), key=f"{key}:{flag}")
                            guided_profiles.append(entry)
                    from mdscope.studio.model import statistics_selector_fields
                    import pandas as pd
                    advanced = {}
                    for family in studio_metric_selection_defaults():
                        if family == "time_series":
                            continue
                        with st.expander(f"{family} · exact variable selectors"):
                            st.caption("Add one row per exact backend variable identity. Empty table selects none. Dynamic identities remain unverified until supported data/topology inspection.")
                            fields = statistics_selector_fields(family)
                            source_entries = selectors.get(family, [])
                            if family == "profiles":
                                source_entries = [item for item in source_entries if item["profile"] not in known_profiles]
                            source = pd.DataFrame(source_entries, columns=fields)
                            columns = {field: st.column_config.CheckboxColumn(field, default=True)
                                       for field in fields if field in {"inference", "ci", "hypothesis_tests"}}
                            edited = st.data_editor(source, num_rows="dynamic", hide_index=True,
                                                    column_config=columns,
                                                    key=f"stats:table:{family}")
                            advanced[family] = [
                                {key: value for key, value in row.items() if pd.notna(value)}
                                for row in edited.to_dict("records")
                            ]
                    advanced["profiles"] = guided_profiles + [item for item in advanced["profiles"]
                                                               if item["profile"] not in known_profiles]
                    stats["metric_selection"] = {"time_series": selectors["time_series"], **advanced}
                    controls = yaml_editor("CI / block / test settings", {key: val for key, val in stats.items()
                                                                         if key not in {"enabled", "metric_selection"}}, "stats:settings") or {}
                    if set(controls) & {"enabled", "metric_selection"}:
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
                    if exact_figures:
                        state.choices["statistics"] = figure_picker(figures, "figure:statistics", state.choices.get("statistics", (doc.get("plots", {}).get("include") or []) if doc.get("plots", {}).get("enabled") else []))
            except (ValueError, TypeError, KeyError, yaml.YAMLError) as exc:
                errors.append(str(exc)); st.error(str(exc))
    with tabs[4]:
        try:
            text, plan = review(doc, base, state.choices if exact_figures else None)
            st.warning("Schema-only validation: topology, chemistry, trajectory alignment and sampling eligibility remain unverified.")
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
            valid = not errors and (not require_topology or bool(preflight and preflight[0] == signature))
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
