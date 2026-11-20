"""Guided project inputs; configuration semantics stay in the backend."""
import streamlit as st
from pathlib import Path

from mdscope.config import DEFAULTS
from mdscope.studio.model import add_replica, add_system, delete_system, remove_replica


def path_picker(label, value, key, *, directory=False, optional=False, create_on_run=False, help_text=None):
    """Browse paths on the host running Studio, not the browser machine."""
    pending = st.session_state.pop(f"{key}:pending_path", None)
    if pending is not None:
        st.session_state[key] = pending
        value = pending
    value = st.text_input(label, value, key=key,
                          help=help_text or ("MDscope creates this directory when analyses run." if create_on_run else None))
    if value:
        path = Path(value).expanduser()
        exists = path.is_dir() if directory else path.is_file()
        if exists:
            st.success("Path found on this machine.")
        elif path.exists():
            st.warning("Path exists but is not the expected file/directory type.")
        else:
            st.info("Output directory will be created when MDscope runs.") if create_on_run else (
                st.info("Optional path not present on this machine.") if optional else
                st.error("Required input path not found on this machine.")
            )
    # Keep browser navigation independent from the configured input path.
    # All browser locations are absolute so '..' and parent navigation are stable.
    browse_key = f"{key}:browse_dir"
    pending_dir = st.session_state.pop(f"{key}:pending_dir", None)
    if pending_dir is not None:
        st.session_state[browse_key] = str(Path(pending_dir).expanduser().resolve())
    elif browse_key not in st.session_state:
        candidate = Path(value).expanduser() if value else Path.cwd()
        st.session_state[browse_key] = str(
            (candidate if candidate.is_dir() else candidate.parent).resolve()
        )
    with st.popover("Browse folders" if directory else "Browse files"):
        location = st.text_input(
            "Directory to browse", key=browse_key,
            help="Enter an absolute directory path (e.g., /home/user/data), or navigate using the buttons below.",
        )
        folder = Path(location).expanduser().resolve()
        if str(folder) != st.session_state.get(f"{key}:browse_listing_dir"):
            st.session_state[f"{key}:browse_listing_dir"] = str(folder)
            st.session_state.pop(f"{key}:browse_choice", None)
        if st.button("Go to parent folder", key=f"{key}:browse_parent"):
            st.session_state[f"{key}:pending_dir"] = str(folder.parent)
            st.rerun()
        if not folder.is_dir():
            st.warning("Directory not found or inaccessible. Enter an existing directory path.")
        else:
            if directory and st.button("Use this folder", key=f"{key}:browse_current"):
                st.session_state[f"{key}:pending_path"] = str(folder)
                st.rerun()
            try:
                entries = sorted(
                    folder.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower())
                )[:200]
                options = [str(entry.resolve()) for entry in entries
                           if entry.is_dir() or (not directory and entry.is_file())]
                if options:
                    chosen = st.selectbox(
                        "Choose a folder or file", options,
                        format_func=lambda path: Path(path).name,
                        key=f"{key}:browse_choice",
                    )
                    if Path(chosen).is_dir():
                        if st.button("Open folder", key=f"{key}:browse_open"):
                            st.session_state[f"{key}:pending_dir"] = chosen
                            st.rerun()
                    if directory or Path(chosen).is_file():
                        if st.button("Use selected path", key=f"{key}:browse_use"):
                            st.session_state[f"{key}:pending_path"] = chosen
                            st.rerun()
                else:
                    st.info("No entries in this directory.")
                if len(entries) == 200:
                    st.caption("Showing at most 200 entries. Enter a more specific directory above.")
            except OSError as exc:
                st.error(f"Cannot list directory: {exc}")
    return value


def project_settings(document, restart):
    st.subheader("Project")
    project = document.setdefault("project", {})
    for field, label in (("name", "Project name"), ("output", "Output directory")):
        project[field] = (path_picker(label, project.get(field, DEFAULTS["project"][field]), f"edit:project:{field}", directory=True, create_on_run=True)
                          if field == "output" else st.text_input(label, project.get(field, DEFAULTS["project"][field]), key=f"edit:project:{field}"))
    st.subheader("Simulation conditions")
    st.markdown("**Compare simulation conditions** — each condition uses one topology and one or more independent trajectories.")
    systems = document.setdefault("systems", {})
    for system_id, spec in list(systems.items()):
        with st.expander(f"{spec.get('label') or 'Condition ' + str(list(systems).index(system_id) + 1)} · {len(spec.get('replicates', []))} replica(s)", expanded=(system_id == next(iter(systems)))):
            for field, label in (("label", "Display label"), ("topology", "Topology path")):
                value = (path_picker(label, spec.get(field, ""), f"edit:system:{system_id}:{field}", help_text="Topology used to define the atom identities, order, and structure for every replica in this condition. Examples: topology.tpr, topology.pdb, topology.gro, topology.psf, topology.prmtop. Format support depends on the MDAnalysis readers installed; verify the selected file with topology preflight.")
                         if field == "topology" else st.text_input(label, spec.get(field, ""), key=f"edit:system:{system_id}:{field}", placeholder="e.g., Wild Type, Apo, Ligand Bound", help="Shown in condition headings, figure legends, and reports. Leave blank to use the default condition name."))
                if value or field in spec:
                    spec[field] = value
            replicas = spec.setdefault("replicates", [])
            for index, rep in enumerate(replicas):
                # CLI shorthand remains accepted; a guided edit preserves its path.
                if isinstance(rep, str):
                    rep = replicas[index] = {"trajectory": rep}
                st.markdown(f"#### {rep.get('plot_label') or 'Rep ' + str(index + 1)}")
                left, right = st.columns(2)
                for column, field, label in (
                    (left, "trajectory", "Trajectory path"),
                    (right, "msd_trajectory", "Optional source-frame MSD trajectory"),
                ):
                    with column:
                        value = (path_picker(label, rep.get(field) or "", f"edit:rep:{system_id}:{index}:{field}", optional=field == "msd_trajectory", help_text=("Primary production trajectory. Frames must contain atoms in the same order as the topology. Examples: rep1.xtc, rep1.trr, rep1.dcd, rep1.nc. Supported formats depend on installed MDAnalysis readers." if field == "trajectory" else "Optional trajectory with original, unmodified coordinates for membrane MSD. Must represent the same replica, with compatible atoms and corresponding frames as the primary trajectory. Examples: rep1_unmodified.xtc, rep1_unmodified.trr. Leave blank unless this analysis needs it."))
                                 if field in {"trajectory", "msd_trajectory"} else
                                 st.text_input(label, rep.get(field) or "", key=f"edit:rep:{system_id}:{index}:{field}", placeholder=(f"Rep {index + 1}" if field == "plot_label" else f"rep{index + 1}")))
                    if value:
                        rep[field] = value
                    else:
                        rep.pop(field, None)
                with st.expander("Replica naming (optional)"):
                    for field, label in (("plot_label", "Figure and report label"), ("label", "Advanced: internal replica ID")):
                        value = st.text_input(label, rep.get(field) or "", key=f"edit:rep:{system_id}:{index}:{field}", help=("Appears in figure legends and reports; default is Rep " + str(index + 1) + "." if field == "plot_label" else "Filesystem-safe identifier for output directories and replica matching. Leave unchanged unless necessary."))
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
    if st.button("Add condition", key="edit:system:add"):
        index = 1
        while f"condition{index}" in systems:
            index += 1
        add_system(document, f"condition{index}")
        restart(document)
    with st.expander("Advanced sampling & time settings"):
        st.caption("Choose which trajectory frames to process. Leave start and stop blank to use the entire trajectory.")
        trajectory = document.setdefault("trajectory", {})
        left, right = st.columns(2)
        for column, field in zip((left, right), ("stride", "chunk_frames")):
            with column:
                trajectory[field] = st.number_input({"stride": "Frame stride (frames)", "chunk_frames": "Processing chunk size (frames)"}[field], min_value=1,
                    value=trajectory.get(field, DEFAULTS["trajectory"][field]), step=1, key=f"edit:trajectory:{field}")
        st.caption("Frame indices start at 0. For example: start 100, stop 1000, stride 10 processes frames 100, 110, …, 990.")
        for column, field in zip((left, right), ("start", "stop")):
            with column:
                value = st.text_input({"start": "Start frame index (inclusive; optional)", "stop": "Stop frame index (exclusive; optional)"}[field], str(trajectory.get(field)) if trajectory.get(field) is not None else "", key=f"edit:trajectory:{field}")
            trajectory[field] = int(value) if value.strip() else None
        time = document.setdefault("time", {})
        value = st.text_input("Override time between frames (ps; optional)", str(time.get("dt_ps")) if time.get("dt_ps") is not None else "", key="edit:time:dt_ps")
        time["dt_ps"] = float(value) if value.strip() else None
        time["origin_ps"] = st.number_input("Initial simulation time (ps)", value=float(time.get("origin_ps", 0.)), key="edit:time:origin_ps")
