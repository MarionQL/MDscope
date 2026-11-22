"""Optional localhost configuration wizard. Run via mdscope-studio."""
import base64
from copy import deepcopy
from pathlib import Path

import streamlit as st
import yaml

from mdscope.config import DEFAULTS, analysis_groups, resolve_config
from mdscope.statistics.selection import (
    ESTIMATOR_FLAGS,
    studio_metric_selection_defaults,
)
from mdscope.studio.model import (
    analysis_options,
    change_group,
    import_document,
    new_document,
    review,
    topology_preflight,
)


def _studio_logo_markup(logo_file: Path) -> str:
    """Optional header asset; never make Studio unavailable for missing artwork."""
    try:
        encoded = base64.b64encode(logo_file.read_bytes()).decode("ascii")
    except OSError:
        return ""
    return (
        '<img class="mdscope-mark" alt="Brown Lab logo" '
        f'src="data:image/png;base64,{encoded}">'
    )


def saved_readiness_details(inventory):
    """Use the shared inspection reasons in each downstream selection tree."""
    if not inventory:
        return
    with st.expander("Saved input readiness and reasons"):
        fields = ("analysis", "run_id", "condition", "replica", "scheduled", "status",
                  "reason", "statistics_ready", "statistics_reason")
        st.dataframe([{key: row.get(key) for key in fields if key in row}
                      for row in inventory], hide_index=True)


def display_metric(name):
    """Human-facing name; keep the original backend identifier in YAML."""
    known = {"rmsd": "RMSD", "rmsf": "RMSF", "rdf": "RDF", "sasa": "SASA",
             "msd": "MSD", "rg": "Radius of gyration"}
    units = {"A": "Å", "ps": "ps", "ns": "ns", "deg": "°", "nm": "nm"}
    parts = name.split("_")
    suffix = units.get(parts[-1])
    core = parts[:-1] if suffix else parts
    readable = " ".join(known.get(part.lower(), part.replace("-", " ").capitalize())
                        for part in core)
    return f"{readable} ({suffix})" if suffix else readable



# Presentation-only grouping. Backend schemas/catalogs remain authoritative.
ANALYSIS_CATEGORIES = (
    ("Protein structure & dynamics", ("alignment", "rmsd", "rmsf", "radius_of_gyration", "dssp", "helix", "tilt", "helix_contacts")),
    ("Membrane structure & dynamics", ("membrane", "membrane_order", "leaflet_normal", "lipid_alignment", "lipid_contacts", "lipid_diffusion", "local_geometry", "residue_membrane_distance")),
    ("Solvation & ion coordination", ("solvation", "hydration", "water", "ion", "rdf", "sasa")),
    ("Molecular interactions", ("hydrogen_bonds", "salt_bridges", "lipid_hydrogen_bonds", "lipid_contacts", "lipid_lifetimes", "lifetimes", "interactions")),
    ("Glycans & carbohydrates", ("hs", "glycan", "sugar")),
    ("Advanced & specialized", ("coupling", "coupling_screen", "density_1d", "sdcc_mechanisms", "position_bridge_audit")),
)


def analysis_category(name):
    """One stable presentation category; no backend option is hidden."""
    key = name.lower()
    for category, hints in ANALYSIS_CATEGORIES:
        if key in hints:
            return category
    for category, hints in ANALYSIS_CATEGORIES:
        if any(key.startswith(hint + "_") or key.endswith("_" + hint) for hint in hints):
            return category
    return "Other analyses"


def grouped_analyses(names):
    """Preserve backend ordering within scientific presentation categories."""
    groups = {}
    for name in names:
        groups.setdefault(analysis_category(name), []).append(name)
    return groups


def eligible_statistics_families(source_status, families):
    """Show only families backed by scheduled or verified saved producers."""
    return tuple(family for family in families if
                 source_status.get(family, {}).get("scheduled_sources") or
                 source_status.get(family, {}).get("status") == "verified_existing")


def analysis_title(name):
    return display_metric(name).replace("Radius Of", "Radius of")


def figure_title(fig):
    return fig.get("label", fig["id"]).replace("_A", " (Å)").replace("_ps", " (ps)").replace("_deg", " (°)")


def figure_stat_type(identifier):
    """Classify statistical figures by the exact production plot variant."""
    variant = identifier.rsplit(":", 1)[-1]
    return {
        "bootstrap_table": "Bootstrap & replicate statistics",
        "comparison_table": "Hypothesis tests & condition comparisons",
        "autocorrelation": "Autocorrelation",
    }.get(variant, "Other figures")


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
    """Cache one resolved snapshot of plan and figures for the current draft."""
    from mdscope.studio.model import planning_snapshot
    return planning_snapshot(yaml.safe_load(text), base)


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
        st.write(f"{display_metric(variable['metric'])} · {variable['label']}")
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



def _show_topology_progress(state, systems, groups):
    """Metadata remains visible even while newly edited groups are validating."""
    topology_future = state.get("studio:topology_future")
    group_future = state.get("studio:group_future")
    with st.expander("Topology identifiers and group matches", expanded=True):
        if topology_future is None or not topology_future.done():
            st.info("Reading topology in the background. You can continue editing Studio.")
        else:
            try:
                metadata = topology_future.result()
                for system_id, spec in systems.items():
                    info = metadata.get(system_id, {})
                    st.markdown(f"**{spec.get('label') or system_id}**")
                    if info.get("error"):
                        st.warning(info["error"])
                        continue
                    st.caption(f"{info['atoms']:,} atoms · {info['residues']:,} residues · {info['segments']:,} segments")
                    for name, values in info["identifiers"].items():
                        labels = ", ".join(values[:40])
                        if len(values) > 40:
                            labels += f" (+{len(values) - 40} more)"
                        st.write(f"{name}: {labels or 'Not provided by topology'}")
            except (OSError, ValueError, TypeError) as exc:
                st.error(f"Topology inspection failed: {exc}")

        if groups:
            if group_future is None or not group_future.done():
                st.caption("Checking updated group selections in the background…")
            else:
                try:
                    matches_by_system = group_future.result()
                    for system_id, matches in matches_by_system.items():
                        if "error" in matches:
                            continue
                        st.markdown(f"**Group matches · {system_id}**")
                        for group_id, match in matches.items():
                            if "error" in match:
                                st.warning(f"{group_id}: {match['error']}")
                            elif match["atoms"] == 0:
                                st.warning(f"{group_id}: matches zero atoms")
                            else:
                                st.caption(f"{group_id}: {match['atoms']:,} atoms, {match['residues']:,} residues")
                except (OSError, ValueError, TypeError) as exc:
                    st.error(f"Group inspection failed: {exc}")


@st.fragment(run_every="3s")
def _pending_topology_fragment(state, systems, groups):
    topology_future = state.get("studio:topology_future")
    group_future = state.get("studio:group_future")
    if (topology_future is not None and topology_future.done() and
            (not groups or group_future is not None and group_future.done())):
        st.rerun(scope="app")
    _show_topology_progress(state, systems, groups)


@st.fragment
def _completed_topology_fragment(state, systems, groups):
    _show_topology_progress(state, systems, groups)


def main():
    st.set_page_config(page_title="MDscope Studio", layout="wide")
    header_html = """
    <style>
      :root {
        --mdscope-cyan: #54B6E6;
        --mdscope-cyan-deep: #147EA9;
        --mdscope-violet: #A36CB0;
        --mdscope-slate: #92959C;
        --mdscope-ink: #172431;
      }
      .block-container { max-width: 1280px; padding-top: 3.5rem; padding-bottom: 3rem; }
      h1, h2, h3 { letter-spacing: -.035em; }
      h1 { font-weight: 760 !important; }
      [data-testid="stTabs"] [data-baseweb="tab-list"] {
        gap: .25rem; overflow-x: auto; scrollbar-width: thin;
        border-bottom: 1px solid rgba(128,145,165,.25);
      }
      [data-testid="stTabs"] button[role="tab"] {
        font-weight: 620; border-radius: .6rem .6rem 0 0;
        padding: .65rem 1rem; white-space: nowrap;
      }
      [data-testid="stTabs"] button[role="tab"][aria-selected="true"] {
        color: var(--mdscope-cyan-deep) !important;
        background: rgba(84,182,230,.085);
        border-bottom-color: var(--mdscope-cyan) !important;
      }
      [data-testid="stExpander"] {
        border-radius: .85rem; border-color: rgba(128,145,165,.22);
      }
      [data-testid="stAlert"] { border-radius: .75rem; }
      [data-testid="stButton"] button, [data-testid="stPopover"] button {
        border-radius: .65rem; transition: border-color .15s ease;
      }
      [data-testid="stButton"] button:hover { border-color: var(--mdscope-cyan-deep); }
      [data-testid="stMetric"] {
        border: 1px solid rgba(128,145,165,.18); border-radius: .75rem;
        padding: .65rem .9rem;
      }
      .mdscope-identity {
        display: flex; align-items: center; gap: 1rem;
        padding: 1.45rem 1.2rem; margin-top: 1rem; margin-bottom: .8rem;
        border: 1px solid rgba(128,145,165,.22); border-radius: 1rem;
        background: linear-gradient(115deg, rgba(84,182,230,.085), rgba(163,108,176,.06) 58%, transparent);
      }
      .mdscope-mark { flex: none; width: 54px; height: 54px; overflow: visible; }
      .mdscope-eyebrow { font-size: .72rem; letter-spacing: .15em;
        text-transform: uppercase; font-weight: 750; color: var(--mdscope-cyan-deep); }
      .mdscope-name { font-size: clamp(1.5rem, 3vw, 2.1rem); font-weight: 760;
        letter-spacing: -.045em; line-height: 1.2; }
      .mdscope-subtitle { font-size: .86rem; opacity: .7; margin-top: .15rem; }
      .mdscope-mark polygon, .mdscope-mark circle { vector-effect: non-scaling-stroke; }
      @media (max-width: 640px) { .mdscope-identity { gap: .7rem; padding: .9rem; }
        .mdscope-mark { width: 44px; height: 44px; } }
    </style>
    <div class="mdscope-identity">
      __MDSCOPE_LOGO__
      <div>
        <div class="mdscope-eyebrow">Brown Lab · Scientific workspace</div>
        <div class="mdscope-name">MDscope Studio</div>
        <div class="mdscope-subtitle">Configure analyses · Define statistics · Build figures · Validate and export</div>
      </div>
    </div>
    """
    logo_file = Path(__file__).parent / "assets" / "MM+Lab+Logo.png"
    logo_markup = _studio_logo_markup(logo_file)
    st.markdown(header_html.replace("__MDSCOPE_LOGO__", logo_markup), unsafe_allow_html=True)
    st.info("Use Studio to set up your systems, choose analyses, configure statistics and figures, and export a validated YAML configuration. To process trajectories, run MDscope separately using the exported file.")
    state = st.session_state
    state.setdefault("document", new_document())
    state.setdefault("choices", {})
    with st.expander("Workspace, existing results and YAML import", expanded=False):
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
            st.caption("Verification checks saved configuration and input paths, "
                       "completed stage markers, replica summaries, and registered artifact "
                       "schemas and summary time/checkpoint identity. Source-file content hashes "
                       "and chemistry/convergence are not established. Offline CLI "
                       "finalization rechecks inputs before publishing.")
            from mdscope.studio.model import results_coverage
            coverage = results_coverage(inventory)
            if coverage:
                st.write("**Verified coverage by analysis, run and condition**")
                st.dataframe(coverage, hide_index=True)
                if all(row["status"] == "verified_existing" for row in coverage):
                    st.caption("For saved-data-only execution, use "
                               "mdscope statistics analysis.yaml --results-dir PATH "
                               "or mdscope plotting analysis.yaml --results-dir PATH. "
                               "The CLI revalidates the files before finalization.")
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
    tabs = st.tabs(["1 · Project & systems", "2 · Analyses",
                    "3 · Statistics", "4 · Figures", "5 · Review & export"])
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
        except (ValueError, TypeError, OSError, yaml.YAMLError) as exc:
            errors.append(str(exc)); st.error(str(exc))
        st.subheader("Reusable measurement groups")
        st.caption("Define a selection once and reuse it across analyses. For example, Protein A can select `protein and segid PROA`. Groups have stable YAML IDs and editable display names.")
        try:
            doc.setdefault("selections", {}).setdefault("groups", {})
            groups = doc["selections"]["groups"]
            def next_group_id(existing, prefix="group"):
                number = 1
                while f"{prefix}_{number}" in existing:
                    number += 1
                return f"{prefix}_{number}"

            for group_id in list(groups):
                if isinstance(groups[group_id], str):
                    groups[group_id] = {"selection": groups[group_id]}
                spec = groups[group_id]
                with st.expander(spec.get("label") or group_id):
                    label = st.text_input("Group name", spec.get("label") or group_id,
                                          key=f"edit:group:{group_id}:label",
                                          help="Human-readable label for this reusable selection. Analyses reference its stable YAML group ID (shown below); renaming the label does not change which atoms are selected.")
                    if label and label != (spec.get("label") or group_id):
                        spec["label"] = label
                    st.caption(f"Reusable YAML group ID: `{group_id}`")
                    measurement = st.text_input("Atoms to measure · MDAnalysis selection",
                                                spec.get("selection") or "",
                                                key=f"edit:group:{group_id}:selection",
                                                help="MDAnalysis selection syntax, e.g. protein and segid PROA.")
                    spec["selection"] = measurement
                    if not measurement.strip():
                        st.warning("Choose a specific atom selection for this group.")
                    left, right = st.columns(2)
                    if left.button("Duplicate group", key=f"edit:group:{group_id}:Duplicate"):
                        target = next_group_id(groups)
                        change_group(doc, group_id, target, duplicate=True, choices=state.choices)
                        doc["selections"]["groups"][target]["label"] = f"{spec.get('label') or group_id} copy"
                        restart_editor(doc)
                    if right.button("Delete group", key=f"edit:group:{group_id}:Delete"):
                        change_group(doc, group_id, None, choices=state.choices)
                        restart_editor(doc)
            if st.button("Add group", key="edit:new_group:add"):
                from mdscope.studio.model import add_group
                target = next_group_id(groups)
                add_group(doc, target)
                doc["selections"]["groups"][target]["label"] = "New group"
                doc["selections"]["groups"][target]["selection"] = ""
                restart_editor(doc)
            import hashlib

            from mdscope.studio.topology_preview import (
                submit_groups,
                submit_topologies,
                topology_signature,
            )
            systems = doc.get("systems", {})
            signature = topology_signature(systems, base)
            groups_snapshot = deepcopy(doc["selections"]["groups"])
            if state.get("studio:topology_signature") != signature:
                state["studio:topology_signature"] = signature
                state["studio:topology_future"] = submit_topologies(
                    deepcopy(systems), base)
                state.pop("studio:group_signature", None)
            group_signature = hashlib.sha256(
                (signature + yaml.safe_dump(groups_snapshot, sort_keys=True)).encode()
            ).hexdigest()
            if state.get("studio:group_signature") != group_signature:
                state["studio:group_signature"] = group_signature
                state["studio:group_future"] = submit_groups(
                    deepcopy(systems), base, groups_snapshot)
            topology_future = state["studio:topology_future"]
            group_future = state.get("studio:group_future")
            if (topology_future.done() and
                    (not groups_snapshot or group_future is not None and group_future.done())):
                _completed_topology_fragment(state, systems, groups_snapshot)
            else:
                _pending_topology_fragment(state, systems, groups_snapshot)
        except (ValueError, TypeError, OSError, yaml.YAMLError) as exc:
            errors.append(str(exc)); st.error(str(exc))
    with tabs[1]:
        st.subheader("Global analysis selections")
        try:
            selections = doc.setdefault("selections", {})
            common = {key: value for key, value in selections.items() if key != "groups"}
            common = yaml_editor("Global selections", common, "edit:selections")
            from mdscope.studio.model import replace_global_selections
            replace_global_selections(doc, common)
        except (ValueError, TypeError, yaml.YAMLError) as exc:
            errors.append(str(exc)); st.error(str(exc))
        card_previews = []
        card_plan = None
        available_analyses = analysis_options()
        filter_text = st.text_input("Find an analysis", key="studio:analysis_search",
                                    placeholder="Search by analysis name")
        show_enabled = st.checkbox("Show enabled analyses only", key="studio:enabled_only")
        for category, names in grouped_analyses(available_analyses).items():
            visible = [name for name in names
                       if filter_text.lower() in name.lower() and
                       (not show_enabled or doc.get("analyses", {}).get(name, {}).get("enabled", False))]
            if not visible:
                continue
            with st.expander(f"{category} · {len(visible)} analyses", expanded=False):
                for name in visible:
                    schema = available_analyses[name]
                    spec = doc.setdefault("analyses", {}).setdefault(name, {"enabled": DEFAULTS["analyses"].get(name, {}).get("enabled", False)})
                    if isinstance(spec, bool):
                        spec = doc["analyses"][name] = {"enabled": spec}
                    enabled = st.checkbox(f"Enable {analysis_title(name)}", value=spec.get("enabled", False), key=f"enable:{name}")
                    spec["enabled"] = enabled
                    if not enabled:
                        continue
                    with st.expander(f"{analysis_title(name)} · settings", expanded=False):
                        try:
                            if name in {"rmsd", "radius_of_gyration"}:
                                default = spec.get("groups", list(doc["selections"].get("groups", {})))
                                spec["groups"] = st.multiselect("Named measurement groups", list(doc["selections"].get("groups", {})),
                                                                 default=default, key=f"edit:groups:{name}")
                                st.caption("These group IDs are defined under Project & systems. Changing group definitions changes the corresponding measurement atoms; global measurement is also saved. This list filters grouped rows only.")
                                if name == "rmsd":
                                    for group_id in spec["groups"]:
                                        group = doc["selections"]["groups"][group_id]
                                        if isinstance(group, str):
                                            group = doc["selections"]["groups"][group_id] = {"selection": group}
                                        fit = st.text_input(
                                            f"Alignment atoms · {group.get('label') or group_id}",
                                            value=group.get("alignment_selection") or "",
                                            key=f"edit:rmsd:alignment:{group_id}",
                                            help="Optional MDAnalysis selection used to fit this group before RMSD. Leave blank to use the normal alignment behavior.")
                                        if fit.strip():
                                            group["alignment_selection"] = fit
                                        else:
                                            group.pop("alignment_selection", None)
                            edited = analysis_settings(name, spec, doc)
                            doc["analyses"][name] = {"enabled": True, **edited,
                                                      **({"groups": spec["groups"]} if name in {"rmsd", "radius_of_gyration"} and "groups" in spec else {})}
                            with st.expander("Technical option names"):
                                st.caption("Supported keys: " + ", ".join(key for key in schema if key != "enabled"))
                            st.caption("Units: Å for distances, ps for time, degrees for angles. Named runs are available where supported.")
                            card_previews.append((name, st.empty()))
                        except (ValueError, TypeError, KeyError, yaml.YAMLError) as exc:
                            errors.append(str(exc)); st.error(str(exc))
        # Fill card placeholders after all settings have been edited, using
        # one coherent plan/catalog snapshot for this Streamlit rerun.
        if card_previews:
            try:
                card_plan, card_figures = preview_document(
                    yaml.safe_dump(doc, sort_keys=False), base)
                for name, placeholder in card_previews:
                    with placeholder.container():
                        show_outputs(card_plan, name, doc, base)
                        if not card_figures.get(name):
                            st.info("No analysis-only figures verified for this analysis; nothing selected automatically.")
            except (ValueError, TypeError, KeyError, yaml.YAMLError) as exc:
                errors.append(str(exc)); st.error(str(exc))
    with tabs[2]:
        stats = doc.setdefault("statistics", {})
        stats["enabled"] = st.checkbox("Enable statistics", value=stats.get("enabled", DEFAULTS["statistics"]["enabled"]), key="stats:master")
        st.caption("Independent units: simulation replicates. Eligibility requires compatible sampling/identities; CI availability is determined by the backend, never by frame count alone.")
        if stats["enabled"]:
            try:
                # Reuse the current draft plan built for the analysis cards.
                # Rebuild below only after statistic selector edits.
                plan = card_plan
                if plan is None:
                    plan, _ = preview_document(yaml.safe_dump(doc, sort_keys=False), base)
                import pandas as pd

                from mdscope.statistics.selection import validate_metric_selection
                from mdscope.studio.model import (
                    set_statistics_exact_family,
                    set_statistics_family_mode,
                    statistics_family_mode,
                    statistics_selector_fields,
                )
                families = tuple(studio_metric_selection_defaults())
                from mdscope.studio.model import statistics_source_status
                inventory_fresh = (
                    state.get("studio:results_inventory_path") == results_dir
                    and state.get("studio:results_inventory_config") ==
                    yaml.safe_dump(doc, sort_keys=True))
                source_status = statistics_source_status(
                    plan, state.get("studio:results_inventory") if inventory_fresh else None)
                if inventory_fresh:
                    saved_readiness_details(state.get("studio:results_inventory"))
                eligible_families = eligible_statistics_families(source_status, families)
                if not eligible_families:
                    st.info("No eligible statistics sources. Enable an analysis or inspect verified saved results.")
                st.caption("Measurements come from enabled analyses or fully verified saved outputs. Statistical methods are configured below.")
                def set_family(family, checked):
                    if family not in eligible_families:
                        return
                    set_statistics_family_mode(state.document, family, "all" if checked else "none")
                    # The widget values are reset when the family changes mode.
                    state.pop(f"stats:table:{family}", None)
                    state.pop(f"stats:tree:convert:{family}", None)
                def set_all_families():
                    checked = state["stats:tree:all"]
                    for family in eligible_families:
                        set_family(family, checked)
                        state[f"stats:tree:family:{family}"] = checked
                with st.expander("Advanced: select all eligible measurement families"):
                    st.caption("This includes dynamic measurement identities discovered during execution.")
                    st.checkbox("Include all eligible measurements",
                                value=all(statistics_family_mode(doc, family) == "all"
                                          for family in eligible_families) if eligible_families else False,
                                key="stats:tree:all", on_change=set_all_families,
                                disabled=not eligible_families)
                for family in eligible_families:
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
                        if not availability["scheduled_sources"] and not availability.get("existing_coverage"):
                            st.warning("No scheduled or verified saved producer for this family. "
                                       "Select a compatible analysis or inspect existing outputs.")
                        st.checkbox(f"Include all {family.replace('_', ' ')} measurements",
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
                            by_analysis = {}
                            for var in variables:
                                by_analysis.setdefault(var.get("analysis", "Other"), []).append(var)
                            for analysis_name, analysis_variables in by_analysis.items():
                                with st.expander(analysis_name.replace("_", " ").title(), expanded=False):
                                    for var in analysis_variables:
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
                                        if st.checkbox(display_metric(var["metric"]) + " · " + str(var.get("group", ""))
                                                       if family == "time_series" else var.get("label", str(identity)),
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
                st.markdown("**Statistical methods**")
                st.caption("Methods apply to eligible measurements; simulation replicas are independent units.")
                stats["bootstrap"] = st.checkbox("Bootstrap confidence intervals",
                    value=stats.get("bootstrap", True), key="stats:bootstrap",
                    help="Resample blocks and replicas according to the configured hierarchy.")
                stats["ci"] = st.number_input("Confidence level", min_value=0.50,
                    max_value=0.999, step=0.01, value=float(stats.get("ci", 0.95)),
                    key="stats:ci", disabled=not stats["bootstrap"])
                tests = stats.setdefault("hypothesis_tests", {})
                tests["enabled"] = st.checkbox("Hypothesis tests",
                    value=tests.get("enabled", True), key="stats:tests")
                correction = ["fdr_bh", "none"]
                tests["multiple_testing"] = st.selectbox("Multiple-testing correction",
                    options=correction, index=correction.index(tests.get("multiple_testing", "fdr_bh")),
                    format_func=lambda v: "Benjamini–Hochberg FDR" if v == "fdr_bh" else "None",
                    disabled=not tests["enabled"], key="stats:fdr")
                with st.expander("Advanced: resampling, test design and backend settings"):
                    controls = yaml_editor("Advanced statistics YAML",
                        {key: val for key, val in stats.items() if key not in
                         {"enabled", "metric_selection", "interaction_min_occupancy",
                          "bootstrap", "ci", "hypothesis_tests"}}, "stats:settings") or {}
                    if set(controls) & {"enabled", "metric_selection", "interaction_min_occupancy",
                                        "bootstrap", "ci", "hypothesis_tests"}:
                        raise ValueError("Use the statistical methods controls above")
                    stats.update(controls)
                    extra_tests = yaml_editor("Advanced hypothesis test options",
                        {key: val for key, val in tests.items()
                         if key not in {"enabled", "multiple_testing"}},
                        "stats:tests:advanced") or {}
                    if set(extra_tests) & {"enabled", "multiple_testing"}:
                        raise ValueError("Use the hypothesis test controls above")
                    tests.update(extra_tests)
                # Statistics-stage previews must reflect the current rerun's
                # edited metric selection, not the earlier analysis-card plan.
                plan, _ = preview_document(yaml.safe_dump(doc, sort_keys=False), base)
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
                        st.json({key: family[key] for key in
                                 ("family", "metric_selection", "files", "metric_contract")
                                 if key in family})
                from mdscope.planning import authoring_catalog
            except (ValueError, TypeError, KeyError, yaml.YAMLError) as exc:
                errors.append(str(exc)); st.error(str(exc))
    with tabs[3]:
        st.caption("Select figures from the current analysis and statistics configuration. "
                   "Figures are optional; this does not execute plotting.")
        plots = doc.get("plots", {})
        # One hierarchy for imported family switches and exact selections.
        # Imported configurations are not converted until a checkbox is changed.
        st.caption("Figures · select all, whole categories, or individual outputs.")
        from mdscope.config import figure_enabled, merge
        from mdscope.planning import authoring_catalog
        merged = merge(DEFAULTS, doc)
        catalog = authoring_catalog(merged, include_inactive=True)["figures"]
        active_catalog = authoring_catalog(merged)["figures"]
        inventory_fresh = (
            state.get("studio:results_inventory_path") == results_dir
            and state.get("studio:results_inventory_config") ==
            yaml.safe_dump(doc, sort_keys=True))
        if inventory_fresh:
            from mdscope.studio.model import saved_source_catalog
            saved_catalog = saved_source_catalog(doc, base, state.get("studio:results_inventory", []))["figures"]
            for owner, figures in saved_catalog.items():
                known = {fig["id"] for fig in active_catalog.get(owner, [])}
                active_catalog.setdefault(owner, []).extend(fig for fig in figures if fig["id"] not in known)
            for owner, figures in active_catalog.items():
                known = {fig["id"] for fig in catalog.get(owner, [])}
                catalog.setdefault(owner, []).extend(fig for fig in figures if fig["id"] not in known)
            st.caption("Saved-only figure sources require all configured replicas to pass the same input checks used by finalization.")
            saved_readiness_details(state.get("studio:results_inventory"))
        figures_by_family = {}
        for owner, figures in active_catalog.items():
            for fig in figures:
                family = fig["id"].split(":", 1)[0]
                figures_by_family.setdefault(family, []).append((owner, fig))
        if not figures_by_family:
            st.info("No eligible figures. Enable an analysis or statistics producer first.")
        else:
            st.caption("Only currently eligible outputs appear below. Imported selections are preserved until edited.")
        all_ids = list(dict.fromkeys(
            fig["id"] for entries in figures_by_family.values() for _, fig in entries))
        exact_import = plots.get("include") is not None
        resolved_figure_cfg = resolve_config(doc) if not exact_import else None
        initial = set(plots.get("include") or []) if exact_import else {
            identifier for identifier in all_ids
            if figure_enabled(resolved_figure_cfg, identifier)
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
        # Retain the established family-level selection controls for existing
        # workflows and imported YAML; the browsing hierarchy below is scientific.
        with st.expander("Select by output type (advanced)"):
            for family, entries in figures_by_family.items():
                ids = list(dict.fromkeys(fig["id"] for _, fig in entries))
                st.checkbox(f"Select all {family.replace('_', ' ')}",
                            value=all(selected(i) for i in ids),
                            key=f"figure:tree:family:{family}",
                            on_change=lambda ids=ids, family=family: set_figures(
                                ids, state[f"figure:tree:family:{family}"]))
        # Statistical outputs are grouped by inference type before source family.
        # Exact plot IDs stay authoritative; this is presentation only.
        by_category = {}
        for family, entries in figures_by_family.items():
            for owner, fig in entries:
                category = "Statistics" if owner == "statistics" else analysis_category(owner)
                method = figure_stat_type(fig["id"]) if owner == "statistics" else "Analysis figures"
                by_category.setdefault(category, {}).setdefault(method, {}).setdefault(owner, {}).setdefault(family, []).append(fig)

        stat_entries = by_category.get("Statistics", {})
        if stat_entries:
            st.markdown("**Select by statistical output type**")
            st.caption("Toggle an entire class of figures. Bootstrap summaries and hypothesis-test comparisons are independent; autocorrelation is a diagnostic, not a hypothesis test.")
            for method, owners in stat_entries.items():
                ids = list(dict.fromkeys(
                    fig["id"] for families in owners.values()
                    for figs in families.values() for fig in figs))
                key = f"figure:tree:method:{method}"
                # A previously checked type checkbox must update when individual
                # children are modified, without overwriting child choices.
                state[key] = bool(ids) and all(selected(i) for i in ids)
                st.checkbox(f"{method} · {len(ids)} figures",
                            key=key,
                            on_change=lambda ids=ids, key=key: set_figures(ids, state[key]))

        for category, methods in by_category.items():
            count = sum(len(figs) for owners in methods.values()
                        for families in owners.values() for figs in families.values())
            with st.expander(f"{category} · {count} figures"):
                for method, owners in methods.items():
                    method_context = (st.expander(f"{method} · " + str(sum(
                        len(figs) for families in owners.values() for figs in families.values())))
                        if category == "Statistics" else st.container())
                    with method_context:
                        for owner, families in owners.items():
                            with st.expander(analysis_title(owner)):
                                for family, figs in families.items():
                                    ids = list(dict.fromkeys(fig["id"] for fig in figs))
                                    with st.expander(f"{family.replace('_', ' ').title()} · {len(ids)} outputs"):
                                        st.checkbox(f"Select all {analysis_title(owner)} / {method} / {family.replace('_', ' ')}",
                                                    value=all(selected(i) for i in ids),
                                                    key=f"figure:tree:family:{family}:owner:{owner}:method:{method}",
                                                    on_change=lambda ids=ids, owner=owner, family=family, method=method: set_figures(
                                                        ids, state[f"figure:tree:family:{family}:owner:{owner}:method:{method}"]))
                                        for fig in figs:
                                            identifier = fig["id"]
                                            st.checkbox(figure_title(fig), value=identifier in initial,
                                                        key=f"figure:tree:{identifier}", help=identifier,
                                                        on_change=changed_figure)
                                            st.caption(fig["destination"])
        exact_figures = exact_import or bool(state.get("figure:tree:edited"))
        if exact_import and not state.get("figure:tree:edited") and not state.choices:
            state.choices = {
                owner: [fig["id"] for fig in figures if fig["id"] in initial]
                for owner, figures in catalog.items()
            }
        if state.get("figure:tree:edited"):
            from mdscope.studio.model import set_figure_selection
            chosen = {identifier for identifier in all_ids if selected(identifier)}
            state.choices = set_figure_selection(doc, catalog, chosen)
            if not exact_import:
                st.warning("Editing individual figures converts category switches to exact IDs. "
                           "Figures whose identities depend on unavailable topology/data "
                           "may not be listed; review the exported figure list.")
    
    with tabs[4]:
        try:
            missing_groups = [group_id for group_id, spec in doc.get("selections", {}).get("groups", {}).items()
                              if isinstance(spec, dict) and not (spec.get("selection") or "").strip()]
            if missing_groups:
                st.error("Choose a measurement selection for: " + ", ".join(missing_groups))
            text, plan = review(doc, base, state.choices if exact_figures else None, inventory=(state.get("studio:results_inventory") if inventory_fresh else None))
            st.warning("Schema-only validation: topology, chemistry, trajectory alignment and sampling eligibility remain unverified.")
            from mdscope.studio.model import unavailable_exact_statistics
            missing_stats = unavailable_exact_statistics(
                doc, plan,
                state.get("studio:results_inventory") if inventory_fresh else None)
            if missing_stats:
                st.error("Exact statistics requested without scheduled or verified saved sources: " +
                         ", ".join(missing_stats) +
                         ". Enable a compatible analysis or inspect existing results "
                         "for the current configuration and path.")

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
            valid = not errors and not missing_groups and not missing_stats and (not require_topology or bool(preflight and preflight[0] == signature))
            st.download_button("Download analysis.yaml", text, file_name="analysis.yaml", mime="text/yaml", disabled=not valid)
            st.info("Export saves a YAML configuration; it does not overwrite analysis results.")
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
