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


def yaml_editor(label, value, key, *, height=140):
    text = st.text_area(label, yaml.safe_dump(value, sort_keys=False), key=key, height=height)
    return yaml.safe_load(text)


def analysis_settings(name, spec):
    """Edit backend defaults with typed widgets; optional fields stay YAML."""
    defaults = DEFAULTS["analyses"].get(name, {})
    result = {}
    for key, default in defaults.items():
        if key in {"enabled", "groups"}:
            continue
        value = spec.get(key, default)
        label = key.replace("_", " ")
        widget_key = f"edit:setting:{name}:{key}"
        if type(default) is bool and type(value) is bool:
            edited = st.checkbox(label, value=value, key=widget_key)
        elif type(default) in {int, float} and type(value) in {int, float}:
            edited = st.number_input(label, value=value, key=widget_key)
        elif isinstance(default, str) and isinstance(value, str):
            edited = st.text_input(label, value=value, key=widget_key)
        else:
            edited = yaml_editor(label + " (YAML; null = backend default)", value, widget_key, height=80)
        if key in spec or edited != default:
            result[key] = edited
    optional = {key: value for key, value in spec.items()
                if key not in defaults and key not in {"enabled", "groups"}}
    edited = yaml_editor("Additional supported options / named runs (YAML)", optional,
                         f"edit:analysis:{name}") or {}
    if not isinstance(edited, dict) or set(edited) & (set(defaults) | {"enabled", "groups"}):
        raise ValueError("Use the controls above for default options, enabled and groups")
    return {**result, **edited}


def figure_picker(figures, key, existing):
    selected = []
    for fig in figures:
        if st.checkbox(fig["id"], value=fig["id"] in existing, key=f"{key}:{fig['id']}"):
            selected.append(fig["id"])
        st.caption(f"{fig['destination']} · {'Statistics required' if fig['requires_statistics'] else 'Descriptive'}")
    return selected


def show_outputs(plan, name, document, base):
    st.caption("ILLUSTRATIVE — NO MD RESULTS. Placeholders describe future saved data.")
    for variable in plan["observable_contracts"]:
        if variable["analysis"] != name:
            continue
        st.write(f"{variable['metric']} · {variable['label']}")
        st.caption(variable["row_grain"] + "; columns are a subset of a shared table.")
        rows = [{column: (f"t{i} (ps)" if column == "time_ps" else
                          variable["group"] if column == "group" else "example (Å)")
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
    outputs = [item for item in plan["registered_replica_outputs"] if item["analysis"] == name]
    if outputs:
        st.dataframe(outputs)
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
    doc = state.document
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
            for section in ("project", "systems", "trajectory", "time"):
                edited = yaml_editor(section, doc.get(section, DEFAULTS.get(section, {})), f"edit:{section}")
                if not isinstance(edited, dict):
                    raise ValueError(f"{section} must be a mapping")
                doc[section] = edited
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
                                 common or {key: val for key, val in DEFAULTS["selections"].items() if key != "groups"},
                                 "edit:selections")
            if not isinstance(common, dict):
                raise ValueError("Selections must be a mapping")
            if "groups" in common:
                raise ValueError("Use the named group editor for groups")
            selections.update(common)
            groups = selections.setdefault("groups", {})
            for group_id in list(groups):
                if isinstance(groups[group_id], str):
                    groups[group_id] = {"selection": groups[group_id]}
                with st.expander(group_id):
                    spec = groups[group_id]
                    for field, label in (("label", "Display label"), ("selection", "Measurement group selection"),
                                         ("alignment_selection", "Independent alignment selection (optional)")):
                        value = st.text_input(label, spec.get(field) or "", key=f"edit:group:{group_id}:{field}")
                        spec[field] = value if value else None
                    target = st.text_input("New group ID for rename/duplicate", key=f"edit:group:{group_id}:id")
                    for action in ("Rename", "Duplicate", "Delete"):
                        if st.button(action, key=f"edit:group:{group_id}:{action}"):
                            if action != "Delete" and not target:
                                raise ValueError("Enter a new group ID")
                            change_group(doc, group_id, None if action == "Delete" else target,
                                         duplicate=action == "Duplicate")
                            for key in list(state):
                                if key.startswith("edit:groups:"):
                                    del state[key]
                            st.rerun()
            new_id = st.text_input("Group ID", key="edit:new_group")
            if st.button("Add group"):
                if not new_id or new_id in groups:
                    raise ValueError("Enter a unique filesystem-safe group ID")
                candidate = deepcopy(doc)
                candidate["selections"]["groups"][new_id] = {"selection": "all"}
                analysis_groups(candidate)
                groups[new_id] = {"label": new_id, "selection": "all"}
                st.rerun()
            st.info("Not validated (topology unavailable or unchecked). Fitted coordinates apply only to supported structural measurements. Membrane/diffusion use source coordinates.")
        except (ValueError, TypeError, yaml.YAMLError) as exc:
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
                    if "groups" in schema:
                        default = spec.get("groups", list(doc["selections"].get("groups", {})))
                        spec["groups"] = st.multiselect("Named measurement groups", list(doc["selections"].get("groups", {})),
                                                         default=default, key=f"edit:groups:{name}")
                        st.caption("Global measurement is also saved; this list filters grouped rows only.")
                    edited = analysis_settings(name, spec)
                    doc["analyses"][name] = {"enabled": True, **edited,
                                              **({"groups": spec["groups"]} if "groups" in spec else {})}
                    st.caption("Supported keys: " + ", ".join(key for key in schema if key != "enabled"))
                    st.caption("Units: _A = Å; _ps = ps; _deg = degrees. Named runs are supported only where runs is listed.")
                    _, plan = review(doc, base)
                    show_outputs(plan, name, doc, base)
                    figures = [fig for var in plan["observable_contracts"] if var["analysis"] == name
                               for fig in var["figures"] if not fig["requires_statistics"]]
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
                _, plan = review(doc, base)
                exact_stats = (stats.get("metric_selection") is not None or
                               st.checkbox("Select statistics by exact variable", key="stats:exact"))
                if not exact_stats:
                    st.info("Imported legacy statistics are preserved (all eligible variables). Enable exact variable selection to replace this behavior.")
                else:
                    if stats.get("metric_selection") is None:
                        stats["metric_selection"] = studio_metric_selection_defaults()
                    selectors = stats["metric_selection"]
                    entries = []
                    figures = []
                    for var in plan["observable_contracts"]:
                        selector = var["statistics_selector"]
                        key = f"stats:{var['metric']}:{var['group']}"
                        selected = next((item for item in selectors.get("time_series", [])
                                         if item.get("group") == selector["group"] and item["metric"] == selector["metric"]), None)
                        if st.checkbox(f"Compute {var['metric']} statistics · {var['label']}", value=selected is not None, key=key):
                            entry = dict(selector)
                            for flag in sorted(ESTIMATOR_FLAGS["time_series"]):
                                entry[flag] = st.checkbox(flag, value=(selected or {}).get(flag, True), key=f"{key}:{flag}")
                            entries.append(entry)
                            figures.extend(fig for fig in var["figures"] if fig["requires_statistics"])
                    known = {(var["metric"], var["group"]) for var in plan["observable_contracts"]}
                    selectors["time_series"] = entries + [item for item in selectors.get("time_series", [])
                                                             if (item["metric"], item.get("group")) not in known]
                    from mdscope.studio.model import statistics_selector_fields
                    import pandas as pd
                    advanced = {}
                    for family in studio_metric_selection_defaults():
                        if family == "time_series":
                            continue
                        with st.expander(f"{family} · exact variable selectors"):
                            st.caption("Add one row per exact backend variable identity. Empty table selects none. Dynamic identities remain unverified until supported data/topology inspection.")
                            fields = statistics_selector_fields(family)
                            source = pd.DataFrame(selectors.get(family, []), columns=fields)
                            columns = {field: st.column_config.CheckboxColumn(field, default=True)
                                       for field in fields if field in {"inference", "ci", "hypothesis_tests"}}
                            edited = st.data_editor(source, num_rows="dynamic", hide_index=True,
                                                    column_config=columns,
                                                    key=f"stats:table:{family}")
                            advanced[family] = [
                                {key: value for key, value in row.items() if pd.notna(value)}
                                for row in edited.to_dict("records")
                            ]
                    stats["metric_selection"] = {"time_series": selectors["time_series"], **advanced}
                    controls = yaml_editor("CI / block / test settings", {key: val for key, val in stats.items()
                                                                         if key not in {"enabled", "metric_selection"}}, "stats:settings") or {}
                    if set(controls) & {"enabled", "metric_selection"}:
                        raise ValueError("Use the statistics selectors above")
                    stats.update(controls)
                    st.caption("Dynamic metric identities outside the structural scalar catalog require exact backend IDs. No identity or output schema is invented.")
                    with st.expander("Expected statistics outputs and eligibility"):
                        st.caption("ILLUSTRATIVE — NO MD RESULTS. Candidate files are conditional on eligible data; complete writer schemas are not yet catalogued.")
                        for family in plan["statistics_families"]:
                            st.json({key: family[key] for key in ("family", "files", "metric_contract")
                                     if key in family})
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
                st.success("Topology preflight passed for the current document; trajectory contents remain unchecked.")
                st.json(preflight[1])
                for group in analysis_groups(resolve_config(doc)):
                    identities = {info["groups"][group]["features"]["ordered_identity_sha256"]
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


if __name__ == "__main__":
    main()
