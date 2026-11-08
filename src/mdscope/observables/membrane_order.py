"""Lipid C-H order parameters and acyl-tail order profiles."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.lib.distances import capped_distance, minimize_vectors

from ..config import atomic_parquet, ensure_dir, time_ps
from ..core.lipid_geometry import classify_leaflets, planar_box
from ..pipeline.progress import progress


AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


@dataclass(frozen=True)
class OrderSite:
    atom_index: int
    residue_index: int
    resname: str
    tail: str
    position: int
    atom_name: str


@dataclass(frozen=True)
class OrderProfile:
    profile_id: str
    positions: np.ndarray
    labels: np.ndarray
    site_bins: np.ndarray
    bin_counts: np.ndarray


def _is_hydrogen(atom) -> bool:
    name = str(atom.name).upper()
    if name.startswith("H") or (len(name) > 1 and name[0].isdigit() and name[1] == "H"):
        return True
    try:
        return str(atom.element).upper() == "H"
    except (AttributeError, NoDataError):
        return False


def _valid_box(universe):
    dims = np.asarray(universe.dimensions, float) if universe.dimensions is not None else None
    if dims is None or len(dims) < 3 or not np.all(np.isfinite(dims[:3])):
        return None
    return universe.dimensions if np.all(dims[:3] > 0) else None


def _tail_atoms(residue, spec: dict) -> list[tuple[object, int]]:
    names = spec.get("carbon_names")
    if names:
        positions = spec.get("positions")
        if positions is not None and len(positions) != len(names):
            raise ValueError("membrane_order tail positions must match carbon_names length")
        lookup = {str(name): i for i, name in enumerate(names)}
        found = {str(atom.name): atom for atom in residue.atoms if str(atom.name) in lookup}
        return [
            (found[str(name)], int(positions[i]) if positions is not None else i + 1)
            for i, name in enumerate(names)
            if str(name) in found
        ]

    pattern = re.compile(str(spec["carbon_pattern"]))
    found = []
    for atom in residue.atoms:
        match = pattern.fullmatch(str(atom.name))
        if not match:
            continue
        if match.lastindex:
            try:
                position = int(match.group(1))
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "membrane_order carbon_pattern capture group 1 must be an integer"
                ) from exc
        else:
            digits = re.findall(r"\d+", str(atom.name))
            if not digits:
                raise ValueError(
                    "membrane_order carbon_pattern needs an integer capture group "
                    "or numeric atom names"
                )
            position = int(digits[-1])
        found.append((atom, position))
    return sorted(found, key=lambda item: item[1])


def _tail_spec(spec: dict, resname: str, tail: str) -> dict:
    base = dict((spec.get("tails", {}) or {}).get(tail, {}) or {})
    override = (
        (spec.get("lipid_overrides", {}) or {})
        .get(str(resname), {})
        .get(tail, {})
        or {}
    )
    base.update(override)
    return base


def _candidate_sites(universe, spec: dict) -> list[OrderSite]:
    selection = spec.get("lipid_selection")
    atoms = universe.select_atoms(str(selection)) if selection else universe.atoms
    excluded = {str(name) for name in spec.get("exclude_resnames", [])}
    tails = list((spec.get("tails", {}) or {}).keys())
    minimum = int(spec.get("min_carbons_per_tail", 4))
    sites = []
    # Iterate selected residues only: full hydrated TPRs contain hundreds of
    # thousands of water residues and scanning all of them is unnecessary.
    for residue in atoms.residues:
        if str(residue.resname) in excluded:
            continue
        for tail in tails:
            tail_spec = _tail_spec(spec, str(residue.resname), tail)
            if not tail_spec.get("carbon_pattern") and not tail_spec.get("carbon_names"):
                continue
            found = _tail_atoms(residue, tail_spec)
            if len(found) < minimum:
                continue
            sites.extend(
                OrderSite(
                    int(atom.index),
                    int(residue.ix),
                    str(residue.resname),
                    str(tail),
                    int(position),
                    str(atom.name),
                )
                for atom, position in found
            )
    return sites


def _fixed_ch_pairs(universe, sites: list[OrderSite], cutoff: float):
    carbon_atoms = universe.atoms[[site.atom_index for site in sites]]
    residue_indices = sorted({site.residue_index for site in sites})
    residue_mask = np.isin(universe.atoms.resindices, residue_indices)
    hydrogen_indices = [
        int(atom.index)
        for atom in universe.atoms[residue_mask]
        if _is_hydrogen(atom)
    ]
    if not hydrogen_indices:
        raise ValueError(
            "membrane_order found no explicit lipid hydrogens; S_CH requires an all-atom topology"
        )
    hydrogens = universe.atoms[hydrogen_indices]
    box = _valid_box(universe)
    # Prefer chemically exact C-H connectivity when topology bonds exist.
    # A radius search can otherwise spuriously associate nearby non-bonded H atoms.
    try:
        topology_bonds=np.asarray(universe.bonds.indices,dtype=np.int64)
    except (NoDataError,AttributeError):
        topology_bonds=None
    if topology_bonds is not None and len(topology_bonds):
        ci=np.full(len(universe.atoms),-1,dtype=np.int32)
        hi=np.full(len(universe.atoms),-1,dtype=np.int32)
        ci[carbon_atoms.indices]=np.arange(len(carbon_atoms),dtype=np.int32)
        hi[hydrogens.indices]=np.arange(len(hydrogens),dtype=np.int32)
        forward=np.stack([ci[topology_bonds[:,0]],hi[topology_bonds[:,1]]],axis=1)
        reverse=np.stack([ci[topology_bonds[:,1]],hi[topology_bonds[:,0]]],axis=1)
        pairs=np.concatenate([forward,reverse])
        pairs=pairs[np.all(pairs>=0,axis=1)]
        if not len(pairs):
            raise ValueError("No bonded C-H pairs in selected lipid tails")
    else:
        # Explicit fallback for PDBs without bond connectivity.
        pairs=capped_distance(
            carbon_atoms.positions,hydrogens.positions,max_cutoff=float(cutoff),
            box=box,return_distances=False)
        if not len(pairs):
            raise ValueError("membrane_order found no C-H pairs within hydrogen_cutoff_A")
    same_residue=carbon_atoms.resindices[pairs[:,0]] == hydrogens.resindices[pairs[:,1]]
    pairs=np.asarray(pairs[same_residue],dtype=np.int32)
    if not len(pairs):
        raise ValueError("membrane_order found no same-residue C-H pairs")

    bonded_counts = np.bincount(pairs[:, 0], minlength=len(sites))
    keep = bonded_counts > 0
    if not np.any(keep):
        raise ValueError("membrane_order found no tail carbons carrying hydrogens")
    if not np.all(keep):
        old_to_new = np.full(len(sites), -1, dtype=np.int32)
        old_to_new[np.flatnonzero(keep)] = np.arange(np.count_nonzero(keep), dtype=np.int32)
        pair_keep = keep[pairs[:, 0]]
        pairs = pairs[pair_keep]
        pairs[:, 0] = old_to_new[pairs[:, 0]]
        sites = [site for site, use in zip(sites, keep) if use]
        carbon_atoms = universe.atoms[[site.atom_index for site in sites]]
    return sites, carbon_atoms, hydrogens, pairs


def _profile_definition(
    sites: list[OrderSite],
    tail: str,
    resname: str | None = None,
    leaflet: str | None = None,
    leaflet_by_residue: dict[int, str] | None = None,
) -> OrderProfile | None:
    selected = np.array(
        [
            site.tail == tail and (resname is None or site.resname == resname)
            and (leaflet is None or (leaflet_by_residue or {}).get(site.residue_index)==leaflet)
            for site in sites
        ],
        dtype=bool,
    )
    if not np.any(selected):
        return None
    positions = np.array(sorted({site.position for site, use in zip(sites, selected) if use}), int)
    position_to_bin = {int(position): i for i, position in enumerate(positions)}
    site_bins = np.full(len(sites), -1, dtype=np.int32)
    for i, (site, use) in enumerate(zip(sites, selected)):
        if use:
            site_bins[i] = position_to_bin[site.position]
    valid_bins = site_bins[site_bins >= 0]
    counts = np.bincount(valid_bins, minlength=len(positions)).astype(float)
    parts = ["lipid_order"] + ([resname] if resname is not None else []) + [tail] + ([leaflet] if leaflet is not None else [])
    profile_id = "__".join(parts)
    return OrderProfile(
        profile_id,
        positions.astype(float),
        positions.astype(str),
        site_bins,
        counts,
    )


def _profiles(sites: list[OrderSite], spec: dict, leaflet_by_residue=None) -> list[OrderProfile]:
    tails = sorted({site.tail for site in sites})
    profiles = [
        profile
        for tail in tails
        if (profile := _profile_definition(sites, tail)) is not None
    ]
    if spec.get("per_lipid_type", True):
        resnames = sorted({site.resname for site in sites})
        profiles.extend(
            profile
            for resname in resnames
            for tail in tails
            if (profile := _profile_definition(sites, tail, resname)) is not None
        )
    if spec.get("per_leaflet", True) and leaflet_by_residue is not None:
        for leaf in ("upper", "lower"):
            profiles.extend(
                profile for tail in tails
                if (profile := _profile_definition(sites, tail, leaflet=leaf,
                            leaflet_by_residue=leaflet_by_residue)) is not None
            )
            if spec.get("per_lipid_type", True):
                profiles.extend(
                    profile for name in sorted({site.resname for site in sites})
                    for tail in tails
                    if (profile := _profile_definition(sites, tail, name, leaf,
                                  leaflet_by_residue)) is not None
                )
    return profiles


def _profile_values(site_values: np.ndarray, profile: OrderProfile) -> np.ndarray:
    valid = (profile.site_bins >= 0) & np.isfinite(site_values)
    sums = np.bincount(
        profile.site_bins[valid],
        weights=site_values[valid],
        minlength=len(profile.positions),
    )
    counts = np.bincount(
        profile.site_bins[valid],
        minlength=len(profile.positions),
    ).astype(float)
    return np.divide(
        sums,
        counts,
        out=np.full(len(profile.positions), np.nan, dtype=float),
        where=counts > 0,
    )


def _atomic_npz(path: Path, **arrays) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as handle:
        np.savez_compressed(handle, **arrays)
    tmp.replace(path)


def _metric_name(profile_id: str) -> str:
    parts = profile_id.split("__")
    return "lipid_order_" + "_".join(parts[1:]) + "_mean"


def run_membrane_order(rep, cfg: dict, repout: Path) -> dict:
    """Calculate conventional -S_CH from explicit C-H bond orientations."""
    spec = cfg["analyses"]["membrane_order"]
    trajectory = cfg["trajectory"]
    axis = str(spec.get("axis", "z")).lower()
    axis_index = AXIS_INDEX[axis]
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    first_frame = trajectory.get("start")
    universe.trajectory[0 if first_frame is None else int(first_frame)]
    sites = _candidate_sites(universe, spec)
    if not sites:
        raise ValueError(
            "membrane_order found no tail carbons; check lipid_selection, tail patterns, "
            "and lipid_overrides"
        )
    sites, carbons, hydrogens, pairs = _fixed_ch_pairs(
        universe,
        sites,
        float(spec.get("hydrogen_cutoff_A", 1.25)),
    )
    leaflet_by_residue = None
    if spec.get("per_leaflet", True):
        representative = {}
        for site in sorted(sites, key=lambda item: (item.residue_index, item.position)):
            representative.setdefault(site.residue_index, site.atom_index)
        indices = np.asarray(list(representative), int)
        coords = universe.atoms[list(representative.values())].positions.copy()
        anchor_selection = spec.get("leaflet_anchor_selection")
        if anchor_selection:
            selected = universe.select_atoms(str(anchor_selection))
            for i, index in enumerate(indices):
                matching = selected[selected.resindices == index]
                if len(matching):
                    coords[i] = matching.center_of_geometry()
        box = planar_box(universe.dimensions, axis_index)
        leaf = classify_leaflets(coords[:, axis_index], float(box[axis_index]))
        leaflet_by_residue = {int(ix): "upper" if flag else "lower"
                              for ix, flag in zip(indices, leaf)}
    profiles = _profiles(sites, spec, leaflet_by_residue)
    if not profiles:
        raise ValueError("membrane_order found no usable sn-1/sn-2 order profiles")

    profile_values = {profile.profile_id: [] for profile in profiles}
    times = []
    frame_metrics = []
    frames = universe.trajectory[
        trajectory.get("start"):trajectory.get("stop"):int(trajectory["stride"])
    ]
    with progress(frames, total=len(frames), desc=f"{rep.system}/{rep.label}: lipid order",
                  enabled=cfg.get("progress", {}).get("enabled", True)) as iterator:
        for frame_index, ts in enumerate(iterator):
            vectors = hydrogens.positions[pairs[:, 1]] - carbons.positions[pairs[:, 0]]
            box = _valid_box(universe)
            if box is not None:
                vectors = minimize_vectors(vectors, box)
            norms = np.linalg.norm(vectors, axis=1)
            valid = norms > 0
            minus_s_ch = np.full(len(pairs), np.nan, dtype=float)
            cosine = vectors[valid, axis_index] / norms[valid]
            s_ch = 0.5 * (3.0 * cosine * cosine - 1.0)
            minus_s_ch[valid] = -s_ch
            site_sums = np.bincount(
                pairs[:, 0][valid],
                weights=minus_s_ch[valid],
                minlength=len(sites),
            )
            valid_counts = np.bincount(pairs[:, 0][valid], minlength=len(sites)).astype(float)
            site_values = np.divide(
                site_sums,
                valid_counts,
                out=np.full(len(sites), np.nan, dtype=float),
                where=valid_counts > 0,
            )
            row = {"time_ps": time_ps(ts, frame_index, cfg)}
            times.append(row["time_ps"])
            for profile in profiles:
                values = _profile_values(site_values, profile)
                profile_values[profile.profile_id].append(values.astype(np.float32))
                row[_metric_name(profile.profile_id)] = float(np.nanmean(values)) if np.any(np.isfinite(values)) else np.nan
                if profile.profile_id in {"lipid_order__sn1", "lipid_order__sn2"}:
                    position3 = np.flatnonzero(profile.positions == 3)
                    if len(position3):
                        row[f"lipid_order_{profile.profile_id.split('__')[1]}_position3_mean"] = float(values[position3[0]])
            tail_metrics = [
                row[_metric_name(profile.profile_id)]
                for profile in profiles
                if profile.profile_id.count("__") == 1
            ]
            if tail_metrics:
                row["lipid_order_mean"] = float(np.nanmean(tail_metrics))
            frame_metrics.append(row)

    out = ensure_dir(Path(repout) / "profiles" / "membrane_order")
    atomic_parquet(out / "timeseries.parquet", pd.DataFrame(frame_metrics))
    summary = pd.DataFrame(
        [
            {
                "resname": resname,
                "tail": tail,
                "carbon_position": position,
                "atom_name": atom_name,
                "n_lipids": int(
                    sum(
                        site.resname == resname
                        and site.tail == tail
                        and site.position == position
                        for site in sites
                    )
                ),
            }
            for resname, tail, position, atom_name in sorted(
                {
                    (site.resname, site.tail, site.position, site.atom_name)
                    for site in sites
                }
            )
        ]
    )
    atomic_parquet(out / "site_mapping.parquet", summary)
    index_rows = []
    for i, profile in enumerate(profiles):
        filename = f"profile-{i:04d}.npz"
        _atomic_npz(
            out / filename,
            time_ps=np.asarray(times, float),
            values=np.asarray(profile_values[profile.profile_id], np.float32),
            positions=profile.positions,
            labels=profile.labels,
        )
        index_rows.append({"profile": profile.profile_id, "file": filename})
    atomic_parquet(out / "profile_index.parquet", pd.DataFrame(index_rows))
    return {
        "frames": len(times),
        "carbon_sites": len(sites),
        "ch_bonds": len(pairs),
        "profiles": len(profiles),
    }


def membrane_order_timeseries(repout: Path) -> list[pd.DataFrame]:
    path = Path(repout) / "profiles" / "membrane_order" / "timeseries.parquet"
    return [pd.read_parquet(path)] if path.exists() else []


def temporal_membrane_order_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    del cfg
    out = {}
    root = Path(repout) / "profiles" / "membrane_order"
    index_path = root / "profile_index.parquet"
    if not index_path.exists():
        return out
    for row in pd.read_parquet(index_path).itertuples(index=False):
        path = root / str(row.file)
        if not path.exists():
            continue
        with np.load(path) as data:
            out[str(row.profile)] = {
                "times": data["time_ps"].astype(float),
                "values": data["values"].astype(float),
                "positions": data["positions"].astype(float),
                "labels": data["labels"].astype(str),
                "reducer": "mean",
            }
    return out
