"""Guided project inputs; configuration semantics stay in the backend."""
import streamlit as st
from pathlib import Path

from mdscope.config import DEFAULTS
from mdscope.studio.model import add_replica, add_system, delete_system, remove_replica


def path_picker(label, value, key, *, directory=False, optional=False):
    """Browse paths on the HPC host, not the user's browser machine."""
    value = st.text_input(label, value, key=key,
                          help="Path on the machine running MDscope Studio. Click Browse to choose an existing file.")
    if value:
        path = Path(value).expanduser()
        exists = path.is_dir() if directory else path.is_file()
        if exists:
            st.caption("✓ Found on this machine")
        elif path.exists():
            st.warning("Path exists but is not the expected file/directory type.")
        else:
            st.caption("Not yet verified on this machine" if optional else "⚠ Path not found on this machine")
    with st.popover("Browse files" if not directory else "Browse folders"):
        location = st.text_input("Directory to browse", str(Path(value).expanduser().parent) if value else str(Path.cwd()),
                                 key=f"{key}:browse_dir")
        folder = Path(location).expanduser()
        if not folder.is_dir():
            st.warning("Directory not found or inaccessible.")
        else:
            try:
                entries = sorted(folder.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower()))[:200]
                options = [str(entry) for entry in entries
                           if entry.is_dir() or directory or entry.is_file()]
                if options:
                    chosen = st.selectbox("Choose a folder or file", options, key=f"{key}:browse_choice")
                    if st.button("Use selected path", key=f"{key}:browse_use"):
                        if directory or Path(chosen).is_file():
                            st.session_state[key] = chosen
                            st.rerun()
                        else:
                            st.session_state[f"{key}:browse_dir"] = chosen
                            st.rerun()
                else:
                    st.info("No entries in this directory.")
                if len(entries) == 200:
                    st.caption("Showing at most 200 entries. Enter a more specific directory above.")
            except OSError as exc:
                st.error(f"Cannot list directory: {exc}")
    return value


def project_settings(document, restart):
    project = document.setdefault("project", {})
    for field, label in (("name", "Project name"), ("output", "Output directory")):
        project[field] = (path_picker(label, project.get(field, DEFAULTS["project"][field]), f"edit:project:{field}", directory=True)
                          if field == "output" else st.text_input(label, project.get(field, DEFAULTS["project"][field]), key=f"edit:project:{field}"))
    systems = document.setdefault("systems", {})
    for system_id, spec in list(systems.items()):
        with st.expander(f"Condition · {system_id}", expanded=(system_id == next(iter(systems)))):
            for field, label in (("label", "Display label"), ("topology", "Topology path")):
                value = (path_picker(label, spec.get(field, ""), f"edit:system:{system_id}:{field}")
                         if field == "topology" else st.text_input(label, spec.get(field, ""), key=f"edit:system:{system_id}:{field}"))
                if value or field in spec:
                    spec[field] = value
            st.caption("Paths are configuration inputs; files are opened only by an explicit preflight.")
            replicas = spec.setdefault("replicates", [])
            for index, rep in enumerate(replicas):
                # CLI shorthand remains accepted; a guided edit preserves its path.
                if isinstance(rep, str):
                    rep = replicas[index] = {"trajectory": rep}
                st.markdown(f"**Replica {index + 1}**")
                left, right = st.columns(2)
                for column, field, label in (
                    (left, "label", "Replica ID"),
                    (right, "plot_label", "Replica display label"),
                    (left, "trajectory", "Trajectory path"),
                    (right, "msd_trajectory", "Optional source-frame MSD trajectory"),
                ):
                    with column:
                        value = (path_picker(label, rep.get(field) or "", f"edit:rep:{system_id}:{index}:{field}", optional=field == "msd_trajectory")
                                 if field in {"trajectory", "msd_trajectory"} else
                                 st.text_input(label, rep.get(field) or "", key=f"edit:rep:{system_id}:{index}:{field}"))
                    if value:
                        rep[field] = value
                    else:
                        rep.pop(field, None)
                if st.button("Remove replica", key=f"edit:rep:{system_id}:{index}:remove", disabled=len(replicas) <= 1):
                    remove_replica(document, system_id, index)
                    restart(document)
            if st.button("Add replica", key=f"edit:system:{system_id}:add_rep"):
                add_replica(document, system_id)
                restart(document)
            if st.button("Remove condition", key=f"edit:system:{system_id}:remove", disabled=len(systems) <= 1):
                delete_system(document, system_id, choices=st.session_state.choices)
                restart(document)
    new_id = st.text_input("New condition ID", key="edit:system:new")
    if st.button("Add condition", key="edit:system:add"):
        add_system(document, new_id)
        restart(document)
    trajectory = document.setdefault("trajectory", {})
    left, right = st.columns(2)
    for column, field in zip((left, right), ("stride", "chunk_frames")):
        with column:
            trajectory[field] = st.number_input(field.replace("_", " "), min_value=1,
                value=trajectory.get(field, DEFAULTS["trajectory"][field]), step=1, key=f"edit:trajectory:{field}")
    st.caption("Frame start is inclusive; stop is exclusive. Negative indices follow Python slicing. Blank bounds use the trajectory extent.")
    for column, field in zip((left, right), ("start", "stop")):
        with column:
            value = st.text_input(f"Frame {field} (optional)", str(trajectory.get(field)) if trajectory.get(field) is not None else "", key=f"edit:trajectory:{field}")
        trajectory[field] = int(value) if value.strip() else None
    time = document.setdefault("time", {})
    value = st.text_input("Override frame spacing (ps; blank uses trajectory time)", str(time.get("dt_ps")) if time.get("dt_ps") is not None else "", key="edit:time:dt_ps")
    time["dt_ps"] = float(value) if value.strip() else None
    time["origin_ps"] = st.number_input("Time origin (ps)", value=float(time.get("origin_ps", 0.)), key="edit:time:origin_ps")
