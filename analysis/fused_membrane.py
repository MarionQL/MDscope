"""One raw-coordinate trajectory pass for lipid order, normals and enrichment.

This observer is called from the existing streaming loop *before* any alignment.
It writes chunked intermediate results with the same part numbers as the
streaming checkpoint, so an interrupted simulation can resume without
replaying committed frames. Standalone execution is supported when the
regular protein streaming metrics are disabled.
"""
from __future__ import annotations

import json
import os
import re
import time
from collections import defaultdict
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import minimize_vectors

from .config import atomic_parquet, ensure_dir, time_ps
from .lipid_geometry import anchor_groups, aligned_leaflets, centers, periodic_center
from .lipid_enrichment import enrichment_frame
from .bilayer_normal import fit_local_normals
from .membrane_order import (
    AXIS_INDEX, _atomic_npz, _candidate_sites, _fixed_ch_pairs,
    _metric_name, _profile_values, _profiles, _valid_box,
)
from .lipid_geometry import classify_leaflets
from .progress import track


FUSED_ANALYSES = ("membrane_order", "bilayer_normal", "lipid_enrichment")


def enabled_names(cfg):
    return tuple(
        name for name in FUSED_ANALYSES
        if cfg.get("analyses", {}).get(name, {}).get("enabled", False)
    )


def _parts(path, suffix):
    return sorted(path.glob("part-*" + suffix))


def _concat(parts):
    return pd.concat((pd.read_parquet(p) for p in parts), ignore_index=True) if parts else pd.DataFrame()


class FusedMembrane:
    """Stateful, chunk-flushed lipid observers sharing a single Universe."""

    def __init__(self, universe, cfg, repout):
        self.u = universe
        self.cfg = cfg
        self.repout = Path(repout)
        self.names = enabled_names(cfg)
        self.cache = ensure_dir(self.repout / "cache" / "fused_membrane")
        self.rows = defaultdict(list)
        self.profile_frames = defaultdict(list)
        self.profile_times = []
        self.timings = defaultdict(float)
        self.frame_count = 0
        self.anchor_states = {}
        self.geometry = {}
        self.order = None
        if "bilayer_normal" in self.names:
            self.geometry["bilayer_normal"] = self._geometry_setup("bilayer_normal")
        if "lipid_enrichment" in self.names:
            self.geometry["lipid_enrichment"] = self._geometry_setup("lipid_enrichment")
            spec = cfg["analyses"]["lipid_enrichment"]
            selection = str(spec.get("protein_selection") or cfg["selections"]["protein"])
            self.protein = universe.select_atoms(selection)
            if not len(self.protein):
                raise ValueError("lipid_enrichment.protein_selection is empty")
            self.inner = float(spec.get("inner_radius_A", 0.))
            self.outer = float(spec.get("outer_radius_A", 15.))
            if not 0 <= self.inner < self.outer:
                raise ValueError("Require 0 <= inner_radius_A < outer_radius_A")
        if "membrane_order" in self.names:
            self.order = self._order_setup()

    def _geometry_setup(self, name):
        spec = self.cfg["analyses"][name]
        axis = AXIS_INDEX[str(spec.get("axis", "z")).lower()]
        # Identical selections share both anchor mapping and per-frame
        # leaflet classification, even if the two analyses differ otherwise.
        key = json.dumps({
            "axis": axis,
            "selections": spec.get("species_selections") or {},
            "allow_absent": bool(spec.get("allow_absent_species", False)),
        }, sort_keys=True)
        if key not in self.anchor_states:
            species, indices, tracker, initial = anchor_groups(self.u, spec)
            self.anchor_states[key] = {
                "species": species, "indices": indices,
                "tracker": tracker, "initial": initial, "axis": axis,
            }
        state = self.anchor_states[key]
        if name == "bilayer_normal":
            radius = float(spec.get("neighbor_radius_A", 15.))
            minimum = int(spec.get("min_neighbors", 5))
            if radius <= 0 or minimum < 3:
                raise ValueError("neighbor_radius_A > 0 and min_neighbors >= 3")
        return {"key": key, "state": state, "spec": spec}

    def _order_setup(self):
        spec = self.cfg["analyses"]["membrane_order"]
        axis = AXIS_INDEX[str(spec.get("axis", "z")).lower()]
        sites = _candidate_sites(self.u, spec)
        if not sites:
            raise ValueError("membrane_order found no tail carbons; check lipid_selection and tail patterns")
        sites, carbons, hydrogens, pairs = _fixed_ch_pairs(
            self.u, sites, float(spec.get("hydrogen_cutoff_A", 1.25))
        )
        leaflet_by_residue = None
        if spec.get("per_leaflet", True):
            representative = {}
            for site in sorted(sites, key=lambda item: (item.residue_index, item.position)):
                representative.setdefault(site.residue_index, site.atom_index)
            indices = np.asarray(list(representative), int)
            coords = self.u.atoms[list(representative.values())].positions.copy()
            anchor_selection = spec.get("leaflet_anchor_selection")
            if anchor_selection:
                selected = self.u.select_atoms(str(anchor_selection))
                for i, index in enumerate(indices):
                    matching = selected[selected.resindices == index]
                    if len(matching):
                        coords[i] = matching.center_of_geometry()
            leaf = classify_leaflets(coords[:, axis], float(self.u.dimensions[axis]))
            leaflet_by_residue = {
                int(ix): "upper" if flag else "lower"
                for ix, flag in zip(indices, leaf)
            }
        profiles = _profiles(sites, spec, leaflet_by_residue)
        if not profiles:
            raise ValueError("membrane_order found no usable sn-1/sn-2 order profiles")
        return {
            "sites": sites, "carbons": carbons, "hydrogens": hydrogens,
            "pairs": pairs, "profiles": profiles, "axis": axis,
        }

    def update(self, ts, frame_index):
        """Consume one raw frame; never changes the Universe coordinates."""
        if not self.names:
            return
        t = time_ps(ts, frame_index, self.cfg)
        box = self.u.dimensions
        shared = {}
        for name in ("bilayer_normal", "lipid_enrichment"):
            if name not in self.geometry:
                continue
            start = time.perf_counter()
            entry = self.geometry[name]
            state = entry["state"]
            key = entry["key"]
            if key not in shared:
                xyz = centers(state["tracker"])
                current, changes = aligned_leaflets(
                    xyz[:, state["axis"]], float(box[state["axis"]]), state["initial"]
                )
                shared[key] = (xyz, current, changes)
            xyz, current, changes = shared[key]
            species, indices = state["species"], state["indices"]
            axis = state["axis"]
            spec = entry["spec"]
            if name == "bilayer_normal":
                angles, neighbors = fit_local_normals(
                    xyz, current, box,
                    float(spec.get("neighbor_radius_A", 15.)),
                    int(spec.get("min_neighbors", 5)), axis,
                )
                row = {
                    "time_ps": t,
                    "bilayer_normal_valid_fraction": float(np.mean(np.isfinite(angles))),
                    "bilayer_normal_leaflet_changes": changes,
                }
                for leaf, mask in (("upper", current), ("lower", ~current)):
                    values = angles[mask]
                    valid = values[np.isfinite(values)]
                    row[f"bilayer_normal_{leaf}_mean_deg"] = float(np.mean(valid)) if len(valid) else np.nan
                    row[f"bilayer_normal_{leaf}_sd_deg"] = float(np.std(valid)) if len(valid) > 1 else np.nan
                    for lipid in np.unique(species):
                        values = angles[mask & (species == lipid)]
                        valid = values[np.isfinite(values)]
                        row[f"bilayer_normal_{lipid}_{leaf}_mean_deg"] = float(np.mean(valid)) if len(valid) else np.nan
                    if spec.get("save_lipid_angles", False):
                        self.rows["angles"].extend(
                            {"time_ps": t, "resindex": int(indices[i]),
                             "species": species[i], "leaflet": leaf,
                             "normal_deviation_deg": float(angles[i]),
                             "n_neighbors": int(neighbors[i])}
                            for i in np.flatnonzero(mask)
                        )
                self.rows[name].append(row)
            else:
                center = periodic_center(self.protein.positions, box)
                vals = enrichment_frame(
                    species, current, xyz, center, box, self.inner, self.outer, axis
                )
                row = {"time_ps": t, "lipid_enrichment_leaflet_changes": changes}
                for value in vals:
                    prefix = f"lipid_enrichment_{value['species']}_{value['leaflet']}"
                    row[prefix + "_fold"] = value["enrichment_fold"]
                    row[prefix + "_local_fraction"] = value["local_fraction"]
                    row[prefix + "_contact_count"] = value["count_local"]
                    self.rows["enrichment_counts"].append({"time_ps": t, **value})
                self.rows[name].append(row)
            self.timings[name] += time.perf_counter() - start

        if self.order is not None:
            start = time.perf_counter()
            o = self.order
            pairs = o["pairs"]
            vectors = o["hydrogens"].positions[pairs[:, 1]] - o["carbons"].positions[pairs[:, 0]]
            valid_box = _valid_box(self.u)
            if valid_box is not None:
                vectors = minimize_vectors(vectors, valid_box)
            norms = np.linalg.norm(vectors, axis=1)
            valid = norms > 0
            minus_s_ch = np.full(len(pairs), np.nan)
            cosine = vectors[valid, o["axis"]] / norms[valid]
            minus_s_ch[valid] = -0.5 * (3. * cosine * cosine - 1.)
            site_sums = np.bincount(
                pairs[:, 0][valid], weights=minus_s_ch[valid],
                minlength=len(o["sites"])
            )
            counts = np.bincount(pairs[:, 0][valid], minlength=len(o["sites"])).astype(float)
            site_values = np.divide(
                site_sums, counts,
                out=np.full(len(o["sites"]), np.nan), where=counts > 0
            )
            row = {"time_ps": t}
            self.profile_times.append(t)
            for index, profile in enumerate(o["profiles"]):
                values = _profile_values(site_values, profile)
                self.profile_frames[index].append(values.astype(np.float32))
                row[_metric_name(profile.profile_id)] = (
                    float(np.nanmean(values)) if np.any(np.isfinite(values)) else np.nan
                )
                if profile.profile_id in {"lipid_order__sn1", "lipid_order__sn2"}:
                    position3 = np.flatnonzero(profile.positions == 3)
                    if len(position3):
                        row[f"lipid_order_{profile.profile_id.split('__')[1]}_position3_mean"] = float(values[position3[0]])
            tail_metrics = [
                row[_metric_name(profile.profile_id)]
                for profile in o["profiles"] if profile.profile_id.count("__") == 1
            ]
            if tail_metrics:
                row["lipid_order_mean"] = float(np.nanmean(tail_metrics))
            self.rows["membrane_order"].append(row)
            self.timings["membrane_order"] += time.perf_counter() - start
        self.frame_count += 1

    def flush(self, part):
        """Persist a complete chunk *before* the streaming checkpoint commits."""
        if not self.names:
            return
        start = time.perf_counter()
        for name in ("bilayer_normal", "lipid_enrichment", "membrane_order"):
            if self.rows[name]:
                atomic_parquet(
                    ensure_dir(self.cache / name) / f"part-{part:06d}.parquet",
                    pd.DataFrame(self.rows[name]),
                )
                self.rows[name].clear()
        for name, key in (("lipid_enrichment", "enrichment_counts"),
                          ("bilayer_normal", "angles")):
            if self.rows[key]:
                atomic_parquet(
                    ensure_dir(self.cache / key) / f"part-{part:06d}.parquet",
                    pd.DataFrame(self.rows[key]),
                )
                self.rows[key].clear()
        if self.profile_times:
            arrays = {"time_ps": np.asarray(self.profile_times, float)}
            for index, profile in enumerate(self.order["profiles"]):
                arrays[f"profile_{index:04d}"] = np.asarray(
                    self.profile_frames[index], np.float32
                ).reshape((-1, len(profile.positions)))
                self.profile_frames[index].clear()
            _atomic_npz(
                ensure_dir(self.cache / "order_profiles") / f"part-{part:06d}.npz",
                **arrays,
            )
            self.profile_times.clear()
        self.timings["chunk_write"] += time.perf_counter() - start

    def prune_from(self, next_part):
        """Drop uncommitted part files left by an interrupted frame chunk."""
        for path in self.cache.rglob("part-*"):
            match = re.search(r"part-(\d+)", path.name)
            if match and int(match.group(1)) >= next_part and path.is_file():
                path.unlink()

def _profile_batches(profiles, frame_count, budget_bytes):
    """Group profile exports by estimated float32 resident bytes.

    A single very large profile is emitted alone; the budget is advisory
    because concatenation and compressed NPZ writing need temporary arrays.
    """
    if budget_bytes <= 0:
        raise ValueError("profile export memory budget must be positive")
    batches, batch, resident = [], [], 0
    for index, profile in enumerate(profiles):
        estimate = int(frame_count) * len(profile.positions) * 4
        if batch and resident + estimate > budget_bytes:
            batches.append(batch)
            batch, resident = [], 0
        batch.append(index)
        resident += estimate
    if batch:
        batches.append(batch)
    return batches


    def finalize(self):
        """Restore original output contracts from committed small chunk files."""
        if not self.names:
            return {}
        start = time.perf_counter()
        result = {}
        if "bilayer_normal" in self.names:
            state = self.geometry["bilayer_normal"]["state"]
            out = ensure_dir(self.repout / "extensions" / "bilayer_normal")
            series = _concat(_parts(self.cache / "bilayer_normal", ".parquet"))
            atomic_parquet(out / "timeseries.parquet", series)
            if self.geometry["bilayer_normal"]["spec"].get("save_lipid_angles", False):
                detail_dir = ensure_dir(out / "per_lipid")
                for path in _parts(self.cache / "angles", ".parquet"):
                    atomic_parquet(detail_dir / ("angles_" + path.name), pd.read_parquet(path))
            result["bilayer_normal"] = {
                "frames": len(series), "n_lipids": len(state["species"])
            }
        if "lipid_enrichment" in self.names:
            state = self.geometry["lipid_enrichment"]["state"]
            out = ensure_dir(self.repout / "extensions" / "lipid_enrichment")
            series = _concat(_parts(self.cache / "lipid_enrichment", ".parquet"))
            atomic_parquet(out / "timeseries.parquet", series)
            atomic_parquet(
                out / "counts.parquet",
                _concat(_parts(self.cache / "enrichment_counts", ".parquet"))
            )
            atomic_parquet(out / "selection_audit.parquet", pd.DataFrame({
                "species": state["species"], "resindex": state["indices"],
                "leaflet": np.where(state["initial"], "upper", "lower"),
            }))
            result["lipid_enrichment"] = {
                "frames": len(series), "n_lipids": len(state["species"]),
                "n_species": len(np.unique(state["species"])),
            }
        if self.order is not None:
            o = self.order
            out = ensure_dir(self.repout / "profiles" / "membrane_order")
            series = _concat(_parts(self.cache / "membrane_order", ".parquet"))
            atomic_parquet(out / "timeseries.parquet", series)
            sites = o["sites"]
            mapping = pd.DataFrame([
                {
                    "resname": resname, "tail": tail, "carbon_position": position,
                    "atom_name": atom_name,
                    "n_lipids": int(sum(
                        site.resname == resname and site.tail == tail
                        and site.position == position for site in sites
                    )),
                }
                for resname, tail, position, atom_name in sorted({
                    (site.resname, site.tail, site.position, site.atom_name) for site in sites
                })
            ])
            atomic_parquet(out / "site_mapping.parquet", mapping)
            parts = _parts(self.cache / "order_profiles", ".npz")
            # Bound peak RAM: for large trajectory exports, loading every
            # profile for every frame at once can swap or exhaust HPC nodes.
            # On ordinary trajectories this forms one batch, retaining the
            # faster single-pass chunk access from the previous optimization.
            limit_mb = max(1, int(os.environ.get(
                "MDSCOPE_ORDER_FINALIZE_MEMORY_MB", "256"
            )))
            batches = _profile_batches(
                o["profiles"], len(series), limit_mb * 1024 * 1024
            )
            index_rows = []
            all_times = np.array([], float)
            for batch_id, indices in enumerate(batches):
                profile_chunks = {i: [] for i in indices}
                time_chunks = []
                for path in parts:
                    with np.load(path, allow_pickle=False) as data:
                        if batch_id == 0:
                            time_chunks.append(data["time_ps"].astype(float))
                        for i in indices:
                            profile_chunks[i].append(
                                data[f"profile_{i:04d}"].astype(np.float32)
                            )
                if batch_id == 0:
                    all_times = (
                        np.concatenate(time_chunks) if time_chunks
                        else np.array([], float)
                    )
                for i in indices:
                    profile = o["profiles"][i]
                    arrays = profile_chunks.pop(i)
                    filename = f"profile-{i:04d}.npz"
                    _atomic_npz(
                        out / filename,
                        time_ps=all_times,
                        values=np.concatenate(arrays, axis=0) if arrays else np.empty(
                            (0, len(profile.positions)), np.float32
                        ),
                        positions=profile.positions,
                        labels=profile.labels,
                    )
                    index_rows.append({"profile": profile.profile_id, "file": filename})
            atomic_parquet(out / "profile_index.parquet", pd.DataFrame(index_rows))
            result["membrane_order"] = {
                "frames": len(series), "carbon_sites": len(sites),
                "ch_bonds": len(o["pairs"]), "profiles": len(o["profiles"]),
            }
        self.timings["finalize"] += time.perf_counter() - start
        from .config import atomic_json
        atomic_json(self.repout / "fused_membrane_timings.json", {
            "frames_processed_this_invocation": self.frame_count,
            "compute_seconds_by_analysis": dict(self.timings),
            "shared_anchor_mappings": len(self.anchor_states),
            "enabled": list(self.names),
        })
        return result


def run_fused_membrane(rep, cfg, repout):
    """One-pass fallback for geometry-only systems without protein streaming."""
    if not enabled_names(cfg):
        return {}
    tr = cfg["trajectory"]
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    start = 0 if tr.get("start") is None else int(tr["start"])
    u.trajectory[start]
    observer = FusedMembrane(u, cfg, repout)
    observer.prune_from(0)
    frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
    chunk = max(1, int(tr.get("chunk_frames", 1000)))
    for k, ts in enumerate(track(
        frames, total=len(frames), desc=f"{rep.system}/{rep.label}: fused membrane",
        enabled=cfg.get("progress", {}).get("enabled", True)
    )):
        observer.update(ts, k)
        if (k + 1) % chunk == 0:
            observer.flush(k // chunk)
    if observer.frame_count % chunk:
        observer.flush(observer.frame_count // chunk)
    return observer.finalize()
