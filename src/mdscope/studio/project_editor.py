"""Guided project inputs; configuration semantics stay in the backend."""
import streamlit as st

from mdscope.config import DEFAULTS
from mdscope.studio.model import add_system, delete_system


def project_settings(document, restart):
    project = document.setdefault("project", {})
    for field, label in (("name", "Project name"), ("output", "Output directory")):
        project[field] = st.text_input(label, project.get(field, DEFAULTS["project"][field]), key=f"edit:project:{field}")
    systems = document.setdefault("systems", {})
    for system_id, spec in list(systems.items()):
        with st.expander(f"Condition · {system_id}", expanded=(system_id == next(iter(systems)))):
            for field, label in (("label", "Display label"), ("topology", "Topology path")):
                value = st.text_input(label, spec.get(field, ""), key=f"edit:system:{system_id}:{field}")
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
                        value = st.text_input(label, rep.get(field) or "", key=f"edit:rep:{system_id}:{index}:{field}")
                    if value:
                        rep[field] = value
                    else:
                        rep.pop(field, None)
                if st.button("Remove replica", key=f"edit:rep:{system_id}:{index}:remove", disabled=len(replicas) <= 1):
                    replicas.pop(index)
                    restart(document)
            if st.button("Add replica", key=f"edit:system:{system_id}:add_rep"):
                used = {(rep.get("label") or f"rep{i}") if isinstance(rep, dict) else f"rep{i}"
                        for i, rep in enumerate(replicas, 1)}
                number = 1
                while f"rep{number}" in used:
                    number += 1
                replicas.append({"label": f"rep{number}", "trajectory": f"rep{number}.xtc"})
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
