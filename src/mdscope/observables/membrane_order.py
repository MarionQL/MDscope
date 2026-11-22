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

from ..config import time_ps
from ..core.persistence import atomic_parquet, ensure_dir
from ..core.lipid_geometry import classify_leaflets, planar_box
from ..core.progress import progress

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


def _tail_atoms(residue, spec: dict) -> list[tuple[int, str, int]]:
    """``(atom index, atom name, carbon position)`` per mapped tail carbon."""
    # Name/index arrays avoid building one Atom object per lipid atom.
    atoms = list(zip(residue.atoms.indices.tolist(), map(str, residue.atoms.names)))
    names = spec.get("carbon_names")
    if names:
        positions = spec.get("positions")
        if positions is not None and len(positions) != len(names):
            raise ValueError("membrane_order tail positions must match carbon_names length")
        lookup = {str(name): i for i, name in enumerate(names)}
        found = {name: index for index, name in atoms if name in lookup}
        return [
            (found[str(name)], str(name), int(positions[i]) if positions is not None else i + 1)
            for i, name in enumerate(names)
            if str(name) in found
        ]

    pattern = re.compile(str(spec["carbon_pattern"]))
    found = []
    for index, name in atoms:
        match = pattern.fullmatch(name)
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
            digits = re.findall(r"\d+", name)
            if not digits:
                raise ValueError(
                    "membrane_order carbon_pattern needs an integer capture group "
                    "or numeric atom names"
                )
            position = int(digits[-1])
        found.append((index, name, position))
    return sorted(found, key=lambda item: item[2])


NAMING_KEYS = ("carbon_pattern", "carbon_names", "positions")


def _tail_spec(spec: dict, resname: str, tail: str) -> dict:
    """Shared tail naming, replaced per lipid by an explicit naming override."""
    base = dict((spec.get("tails", {}) or {}).get(tail, {}) or {})
    override = (
        (spec.get("lipid_overrides", {}) or {})
        .get(str(resname), {})
        .get(tail, {})
        or {}
    )
    if override.get("carbon_pattern") or override.get("carbon_names"):
        base = {key: value for key, value in base.items() if key not in NAMING_KEYS}
    base.update(override)
    return base


def _scan_sites(universe, spec: dict, bonds):
    """Map tail carbons per selected residue and verify them against the topology.

    Returns ``(sites, matched, excluded, chains)``: matched-carbon counts per
    residue for each ``(resname, tail)``, residue counts skipped through
    ``exclude_resnames``, and the bond-graph verification summary.
    """
    selection = spec.get("lipid_selection")
    atoms = universe.select_atoms(str(selection)) if selection else universe.atoms
    excluded_names = {str(name) for name in spec.get("exclude_resnames", [])}
    tails = list((spec.get("tails", {}) or {}).keys())
    minimum = int(spec.get("min_carbons_per_tail", 4))
    sites, matched, excluded = [], {}, {}
    # Iterate selected residues only: full hydrated TPRs contain hundreds of
    # thousands of water residues and scanning all of them is unnecessary.
    for residue in atoms.residues:
        resname = str(residue.resname)
        if resname in excluded_names:
            excluded[resname] = excluded.get(resname, 0) + 1
            continue
        for tail in tails:
            tail_spec = _tail_spec(spec, resname, tail)
            found = (_tail_atoms(residue, tail_spec)
                     if tail_spec.get("carbon_pattern") or tail_spec.get("carbon_names") else [])
            matched.setdefault((resname, str(tail)), []).append(len(found))
            if len(found) < minimum:
                continue
            sites.extend(
                OrderSite(int(index), int(residue.ix), resname, str(tail), int(position), name)
                for index, name, position in found
            )
    return sites, matched, excluded, _verify_tail_chains(universe, sites, bonds)


def _site_label(site: OrderSite) -> str:
    return f"{site.resname} tail {site.tail} atom {site.atom_name} (position {site.position})"


def _check_tail_atoms(universe, sites: list[OrderSite], group_ids, positions) -> None:
    """Each mapped atom is one carbon at one position of one residue tail."""
    owner = {}
    for site in sites:
        other = owner.setdefault(site.atom_index, site)
        if other is not site:
            raise ValueError(
                f"membrane_order maps {site.resname} atom {site.atom_name} to both tails "
                f"{other.tail} and {site.tail}; tail naming must not overlap"
            )
    try:
        elements = universe.atoms[[site.atom_index for site in sites]].elements
    except (NoDataError, AttributeError):
        elements = ()
    for site, element in zip(sites, elements):
        element = str(element).strip().upper()
        if element and element != "C":
            raise ValueError(
                f"membrane_order maps non-carbon atom {_site_label(site)} with "
                f"element {element}; check tails and lipid_overrides"
            )
    seen = set()
    for site, group, position in zip(sites, group_ids, positions):
        if (group, position) in seen:
            raise ValueError(f"membrane_order maps duplicate carbon positions in {_site_label(site)}")
        seen.add((group, position))


def _verify_tail_chains(universe, sites: list[OrderSite], bonds) -> dict:
    """Reject tail mappings contradicted by topology elements or bonds.

    With bonds, carbons of one residue tail must form one chain whose bonded
    neighbours differ by exactly one position, and no tail carbon may bond to
    another tail or residue. Residue tails without any bonded tail-carbon pair
    (C-H-only bond graphs) are counted as unverified instead.
    """
    summary = {"bond_graph": False, "residue_tails": 0, "verified_residue_tails": 0}
    if not sites:
        return summary
    group_of = {}
    group_ids = np.asarray([group_of.setdefault((site.residue_index, site.tail), len(group_of))
                            for site in sites], dtype=np.int64)
    positions = np.asarray([site.position for site in sites], dtype=np.int64)
    summary["residue_tails"] = n_groups = len(group_of)
    _check_tail_atoms(universe, sites, group_ids, positions)
    if bonds is None:
        return summary
    summary["bond_graph"] = True
    site_of = np.full(len(universe.atoms), -1, dtype=np.int64)
    site_of[[site.atom_index for site in sites]] = np.arange(len(sites))
    left, right = site_of[bonds[:, 0]], site_of[bonds[:, 1]]
    keep = (left >= 0) & (right >= 0)
    pairs = np.unique(np.sort(np.stack([left[keep], right[keep]], axis=1), axis=1), axis=0)
    left, right = pairs[:, 0], pairs[:, 1]
    invalid = np.flatnonzero((group_ids[left] != group_ids[right])
                             | (np.abs(positions[left] - positions[right]) != 1))
    if len(invalid):
        a, b = int(left[invalid[0]]), int(right[invalid[0]])
        raise ValueError(
            f"membrane_order tail mapping contradicts topology bonds: "
            f"{_site_label(sites[a])} is bonded to {_site_label(sites[b])}; "
            "bonded tail carbons must be consecutive positions of one tail"
        )
    edges = np.bincount(group_ids[left], minlength=n_groups)
    nodes = np.bincount(group_ids, minlength=n_groups)
    checked = edges > 0
    split = np.flatnonzero(checked & (nodes - edges > 1))
    if len(split):
        group = int(split[0])
        members = [site for site, gid in zip(sites, group_ids) if gid == group]
        raise ValueError(
            f"membrane_order tail {members[0].tail} of {members[0].resname} maps "
            f"{int(nodes[group] - edges[group])} disconnected bonded segments "
            f"({', '.join(site.atom_name for site in members)}); a tail must be one bonded chain"
        )
    summary["verified_residue_tails"] = int(np.count_nonzero(checked))
    return summary


def _fixed_ch_pairs(universe, sites: list[OrderSite], cutoff: float, topology_bonds):
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
    if topology_bonds is not None:
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


def _order_sites(universe, spec: dict):
    """Verified tail sites and fixed C-H pairs; one bond-graph read for both."""
    try:
        bonds = np.asarray(universe.bonds.indices, dtype=np.int64)
    except (NoDataError, AttributeError):
        bonds = None
    if bonds is not None and not len(bonds):
        bonds = None
    scan = _scan_sites(universe, spec, bonds)
    if not scan[0]:
        raise ValueError(
            "membrane_order found no tail carbons; check lipid_selection, tails, "
            "tail_preset and lipid_overrides"
        )
    pairs = _fixed_ch_pairs(universe, scan[0], float(spec.get("hydrogen_cutoff_A", 1.25)), bonds)
    return scan, pairs, bonds is not None


def _initial_leaflets(universe, sites: list[OrderSite], spec: dict, axis_index: int) -> dict[int, str]:
    """Fixed leaflet per order lipid from the current frame.

    With ``leaflet_anchor_selection`` every order lipid needs an anchor atom;
    otherwise each lipid uses its lowest-position tail carbon.
    """
    representative = {}
    for site in sorted(sites, key=lambda item: (item.residue_index, item.position)):
        representative.setdefault(site.residue_index, site.atom_index)
    indices = np.asarray(list(representative), int)
    coords = universe.atoms[list(representative.values())].positions.copy()
    anchor_selection = spec.get("leaflet_anchor_selection")
    if anchor_selection:
        selected = universe.select_atoms(str(anchor_selection))
        missing = []
        for i, index in enumerate(indices):
            matching = selected[selected.resindices == index]
            if len(matching):
                coords[i] = matching.center_of_geometry()
            else:
                missing.append(universe.residues[int(index)])
        if missing:
            examples = ", ".join(f"{res.resname} {res.resid}" for res in missing[:3])
            raise ValueError(
                f"membrane_order.leaflet_anchor_selection matched no atom in {len(missing)} "
                f"order lipids (e.g. {examples}); leaflets need one comparable anchor per lipid"
            )
    box = planar_box(universe.dimensions, axis_index)
    leaf = classify_leaflets(coords[:, axis_index], float(box[axis_index]))
    return {int(ix): "upper" if flag else "lower" for ix, flag in zip(indices, leaf)}


def order_setup(universe, spec: dict) -> dict:
    """Fixed C-H pairs, initial leaflets and profiles from the current frame."""
    axis = AXIS_INDEX[str(spec.get("axis", "z")).lower()]
    _, (sites, carbons, hydrogens, pairs), _ = _order_sites(universe, spec)
    leaflet_by_residue = (_initial_leaflets(universe, sites, spec, axis)
                          if spec.get("per_leaflet", True) else None)
    profiles = _profiles(sites, spec, leaflet_by_residue)
    if not profiles:
        raise ValueError("membrane_order found no usable sn-1/sn-2 order profiles")
    return {"sites": sites, "carbons": carbons, "hydrogens": hydrogens,
            "pairs": pairs, "profiles": profiles, "axis": axis}


def order_preflight(universe, spec: dict, *, leaflets: bool = True,
                    anchor_residues: set[int] | None = None) -> dict:
    """Topology-only audit of tail mapping, C-H pairs and initial leaflets.

    Raises on contradictory mappings or unavailable S_CH (no bonded/near
    hydrogens). Species/tails without enough mapped carbons are reported as
    exclusions, matching what the observer will compute. ``anchor_residues``
    counts order lipids outside a population defined by other anchors.
    """
    (_, matched, excluded, chains), (sites, _, _, pairs), has_bonds = _order_sites(universe, spec)
    minimum = int(spec.get("min_carbons_per_tail", 4))
    ch_bonds = np.bincount(pairs[:, 0], minlength=len(sites))
    rows = []
    for (resname, tail), found in sorted(matched.items()):
        found = np.asarray(found)
        used = [(site, int(n)) for site, n in zip(sites, ch_bonds)
                if site.resname == resname and site.tail == tail]
        positions = sorted({site.position for site, _ in used})
        mapped = len({site.residue_index for site, _ in used})
        rows.append({
            "resname": resname, "tail": tail, "residues": int(len(found)),
            "mapped_residues": mapped,
            "no_ch_residues": int(np.count_nonzero(found >= minimum)) - mapped,
            "below_minimum_residues": int(np.count_nonzero((found > 0) & (found < minimum))),
            "unmatched_residues": int(np.count_nonzero(found == 0)),
            "ch_sites": len(used), "ch_bonds": int(sum(n for _, n in used)),
            "positions": f"{positions[0]}-{positions[-1]}" if positions else None,
        })
    audit = {
        "tail_preset": spec.get("tail_preset"), "min_carbons_per_tail": minimum,
        "ch_pairs": "topology bonds" if has_bonds else "distance fallback",
        "chains": chains, "excluded_resnames": excluded, "species": rows,
        "lipids": len({site.residue_index for site in sites}),
        "unused_overrides": sorted(set(spec.get("lipid_overrides") or {})
                                   - {resname for resname, _ in matched}),
    }
    residues = sorted({site.residue_index for site in sites})
    if anchor_residues is not None:
        audit["lipids_without_anchor"] = len(set(residues) - set(anchor_residues))
    if has_bonds:
        fragments = universe.residues[residues].atoms.fragindices
        owner = universe.residues[residues].atoms.resindices
        per_fragment = {}
        for fragment, residue in zip(fragments, owner):
            per_fragment.setdefault(int(fragment), set()).add(int(residue))
        audit["multi_residue_molecules"] = sum(len(v) > 1 for v in per_fragment.values())
    if leaflets and spec.get("per_leaflet", True):
        if universe.dimensions is None:
            audit["leaflets"] = None  # Checked on the first analyzed frame instead.
        else:
            axis = AXIS_INDEX[str(spec.get("axis", "z")).lower()]
            leaf = _initial_leaflets(universe, sites, spec, axis)
            values = list(leaf.values())
            audit["leaflets"] = {"upper": values.count("upper"), "lower": values.count("lower")}
    return audit


def run_membrane_order(rep, cfg: dict, repout: Path) -> dict:
    """Calculate conventional -S_CH from explicit C-H bond orientations."""
    spec = cfg["analyses"]["membrane_order"]
    trajectory = cfg["trajectory"]
    axis = str(spec.get("axis", "z")).lower()
    axis_index = AXIS_INDEX[axis]
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    first_frame = trajectory.get("start")
    universe.trajectory[0 if first_frame is None else int(first_frame)]
    setup = order_setup(universe, spec)
    sites, carbons, hydrogens = setup["sites"], setup["carbons"], setup["hydrogens"]
    pairs, profiles = setup["pairs"], setup["profiles"]

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
