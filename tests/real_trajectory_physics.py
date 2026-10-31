"""Full-length reference-trajectory physical audit and per-stage timing.

Uses every available frame (200/201) for leaflet/PBC/geometry checks. The
glycolipid H-bond subtest uses 20 frames because it is considerably costlier.
This is NOT a claim of trajectory convergence or independent replicates.
"""
from __future__ import annotations

import gc
import json
import os
import resource
import time
import traceback
from collections import Counter, defaultdict
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import minimize_vectors

from mdscope.observables.bilayer_normal import fit_local_normals
from mdscope.config import Replicate
from mdscope.observables.lipid_enrichment import enrichment_frame
from mdscope.core.lipid_geometry import (
    aligned_leaflets, anchor_groups, centers, periodic_center,
)
from mdscope.interactions.lipid_hydrogen_bonds import (
    lipid_molecule_map, run_lipid_hydrogen_bonds,
)

DATA = Path(os.environ.get("MDSCOPE_REAL_DATA", "real_test_data")).resolve()
ROOT = Path(os.environ.get("MDSCOPE_REAL_OUTPUT", "real_test_results")).resolve()
ROOT.mkdir(parents=True, exist_ok=True)
SYSTEMS = {
    "carb_sdc4": ("carb-sdc4-step7_production.tpr", "carb-sdc4-rep1_test.xtc"),
    "nocarb_sdc4": ("nocarb-sdc4-step7_production.tpr", "nocarb-sdc4-rep1_test.xtc"),
}
SPECIES = {
    name: f"resname {name} and name {anchor}"
    for name, anchor in [
        ("POPC", "P"), ("POPE", "P"), ("POPI", "P"), ("POPS", "P"),
        ("DOPC", "P"), ("DPPC", "P"), ("PSM", "P"), ("CHL1", "O3"),
    ]
}
SPECIES.update({
    name: f"chainID {name} and resname CER160 and name NF"
    for name in ("GM1a", "BGAL", "CER160")
})


def seconds(timer, name, start):
    timer[name] += time.perf_counter() - start


def geometry_audit(system, topology, trajectory):
    begin = time.perf_counter()
    timer = defaultdict(float)
    u = mda.Universe(str(topology), str(trajectory))
    u.trajectory[0]
    spec = {"axis": "z", "species_selections": SPECIES, "allow_absent_species": True}
    start = time.perf_counter()
    species, indices, tracker, initial = anchor_groups(u, spec)
    seconds(timer, "anchor_topology_setup_s", start)
    protein = u.select_atoms("protein and name CA")
    assert len(protein) > 0
    frame_records = []
    max_wrap_enrichment_diff = 0
    max_wrap_normal_diff = 0.
    max_wrap_center_diff = 0.
    changed_counts = []
    changed_lipids = {}
    transition_counts = Counter()
    previous_leaf = None
    previous_time = None
    sampled = {0, 1, 2, 5, 10, 20, 40, 80, 120, 160, len(u.trajectory)-1}
    naive_center_disagreement = 0
    for i, ts in enumerate(u.trajectory):
        t0 = time.perf_counter()
        xyz = centers(tracker)
        seconds(timer, "anchor_coordinates_s", t0)
        t0 = time.perf_counter()
        current, changes = aligned_leaflets(xyz[:, 2], float(ts.dimensions[2]), initial)
        seconds(timer, "leaflet_assignment_s", t0)
        changed_counts.append(changes)
        for j in np.flatnonzero(current != initial):
            key = str(int(indices[j]))
            entry = changed_lipids.setdefault(key, {
                "resindex": int(indices[j]), "species": str(species[j]),
                "n_frames_changed": 0, "first_changed_frame": i,
                "first_changed_time_ps": float(ts.time),
                "first_changed_z_A": float(xyz[j, 2]),
            })
            entry["n_frames_changed"] += 1
            entry["last_changed_frame"] = i
            entry["last_changed_time_ps"] = float(ts.time)
            entry["last_changed_z_A"] = float(xyz[j, 2])
        if previous_leaf is not None:
            for j in np.flatnonzero(current != previous_leaf):
                transition_counts[str(int(indices[j]))] += 1
        previous_leaf = current.copy()
        t0 = time.perf_counter()
        center = periodic_center(protein.positions, ts.dimensions)
        seconds(timer, "periodic_protein_center_s", t0)
        naive = protein.center_of_geometry()
        center_delta = float(np.linalg.norm(minimize_vectors(
            (center-naive)[None, :], ts.dimensions)[0]))
        t0 = time.perf_counter()
        enrichment = enrichment_frame(
            species, current, xyz, center, ts.dimensions, 0., 25., 2)
        seconds(timer, "enrichment_s", t0)
        t0 = time.perf_counter()
        angles, neighbors = fit_local_normals(
            xyz, current, ts.dimensions, radius=15., min_neighbors=4, normal=2)
        seconds(timer, "bilayer_normals_s", t0)
        assert len(angles) == len(species)
        valid = np.isfinite(angles)
        assert valid.mean() > 0.75, f"{system} frame {i}: insufficient local normals"
        assert np.all((angles[valid] >= 0) & (angles[valid] <= 90.))
        assert np.all(np.isfinite([r["n_bulk"] for r in enrichment]))
        assert all(r["n_bulk"] > 0 for r in enrichment)
        if previous_time is not None:
            assert ts.time > previous_time
        previous_time = ts.time
        if i in sampled:
            # Independently shift selected lipids/protein atoms by entire box
            # vectors. All PBC-aware observables must remain invariant.
            shifted_xyz = xyz.copy()
            shifted_xyz[::4, 0] += ts.dimensions[0]
            shifted_xyz[1::4, 1] -= ts.dimensions[1]
            shifted_xyz[2::4, 2] += ts.dimensions[2]
            shifted_protein = protein.positions.astype(float).copy()
            shifted_protein[::3, 0] -= ts.dimensions[0]
            shifted_protein[1::3, 1] += ts.dimensions[1]
            shifted_center = periodic_center(shifted_protein, ts.dimensions)
            center_error = float(np.linalg.norm(minimize_vectors(
                (center-shifted_center)[None, :], ts.dimensions)[0]))
            max_wrap_center_diff = max(max_wrap_center_diff, center_error)
            shifted_leaf, _ = aligned_leaflets(
                shifted_xyz[:, 2], float(ts.dimensions[2]), initial)
            assert np.array_equal(shifted_leaf, current)
            shifted_enrichment = enrichment_frame(
                species, shifted_leaf, shifted_xyz, shifted_center,
                ts.dimensions, 0., 25., 2)
            count_difference = max(abs(a["count_local"]-b["count_local"])
                                   for a, b in zip(enrichment, shifted_enrichment))
            max_wrap_enrichment_diff = max(max_wrap_enrichment_diff, count_difference)
            shifted_angles, _ = fit_local_normals(
                shifted_xyz, shifted_leaf, ts.dimensions, radius=15.,
                min_neighbors=4, normal=2)
            angle_difference = float(np.nanmax(np.abs(angles-shifted_angles)))
            max_wrap_normal_diff = max(max_wrap_normal_diff, angle_difference)
            naive_enrichment = enrichment_frame(
                species, current, xyz, naive, ts.dimensions, 0., 25., 2)
            if any(a["count_local"] != b["count_local"]
                   for a, b in zip(enrichment, naive_enrichment)):
                naive_center_disagreement += 1
        frame_records.append({
            "frame": i, "time_ps": float(ts.time),
            "leaflet_changes": changes,
            "n_upper": int(current.sum()),
            "n_lower": int((~current).sum()),
            "valid_normal_fraction": float(valid.mean()),
            "normal_mean_deg": float(np.nanmean(angles)),
            "naive_center_error_A": center_delta,
            "local_lipid_count": int(sum(r["count_local"] for r in enrichment)),
        })
    assert len(frame_records) == len(u.trajectory)
    assert max_wrap_enrichment_diff == 0
    assert max_wrap_center_diff < 1e-3, max_wrap_center_diff
    assert max_wrap_normal_diff < 1e-3, max_wrap_normal_diff
    print(f"{system}: PBC wrap invariance center={max_wrap_center_diff:.5g} A, "
          f"enrichment count diff={max_wrap_enrichment_diff}, "
          f"local normal diff={max_wrap_normal_diff:.5g} deg", flush=True)
    frames = pd.DataFrame(frame_records)
    out = ROOT / system / "physics"
    out.mkdir(parents=True, exist_ok=True)
    frames.to_parquet(out / "geometry_per_frame.parquet", index=False)
    transition_details = transition_geometry(u, tracker, species, initial, changed_lipids, out)
    timer["geometry_total_s"] = time.perf_counter()-begin
    report = {
        "n_frames": len(frame_records),
        "n_anchors": len(species),
        "species_counts": dict(Counter(species)),
        "max_leaflet_changes": int(max(changed_counts)),
        "frames_with_leaflet_changes": int(sum(c > 0 for c in changed_counts)),
        "transition_geometry": transition_details,
        "changed_lipids": [dict(v, transitions=int(transition_counts[k]))
                           for k, v in changed_lipids.items()],
        "initial_upper_species": dict(Counter(species[initial])),
        "initial_lower_species": dict(Counter(species[~initial])),
        "final_upper_species": dict(Counter(species[previous_leaf])),
        "final_lower_species": dict(Counter(species[~previous_leaf])),
        "max_wrap_center_diff_A": max_wrap_center_diff,
        "max_wrap_enrichment_count_diff": max_wrap_enrichment_diff,
        "max_wrap_normal_diff_deg": max_wrap_normal_diff,
        "naive_center_wrong_on_sampled_frames": naive_center_disagreement,
        "max_naive_center_error_A": float(frames.naive_center_error_A.max()),
        "mean_valid_normal_fraction": float(frames.valid_normal_fraction.mean()),
        "mean_normal_deg": float(frames.normal_mean_deg.mean()),
        "timings_s": dict(timer),
        "rss_peak_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.,
    }
    print(f"{system}: geometry audit = {json.dumps(report, sort_keys=True)}", flush=True)
    del u
    gc.collect()
    return report


def transition_geometry(u, tracker, species, reference, changed_lipids, out):
    """Check sterol O3 crossings against phosphate-defined bilayer midplane.

    A change in the nearest leaflet alone does not establish sterol flip-flop.
    Record both O3 position relative to the membrane midplane and the O3-C17
    molecular orientation for manual inspection.
    """
    if not changed_lipids:
        return {}
    records = []
    phospho = np.isin(species, ["POPC", "POPE", "POPI", "POPS", "DOPC", "DPPC", "PSM"])
    def circ_mean(z, L):
        angle = 2 * np.pi * z / L
        return float(np.arctan2(np.sin(angle).mean(), np.cos(angle).mean())
                     * L / (2*np.pi) % L)
    def signed_delta(z, center, L):
        return float((z-center+L/2) % L - L/2)
    for lipid in changed_lipids.values():
        res = u.residues[int(lipid["resindex"])]
        head = res.atoms.select_atoms("name O3")
        tail = res.atoms.select_atoms("name C17")
        assert len(head) == 1, "Expected one cholesterol O3 anchor"
        assert len(tail) == 1, "Expected one cholesterol C17 orientation atom"
        for frame, ts in enumerate(u.trajectory):
            L = float(ts.dimensions[2])
            xyz = centers(tracker)
            current, _ = aligned_leaflets(xyz[:,2], L, reference)
            top = circ_mean(xyz[phospho & current,2], L)
            bottom = circ_mean(xyz[phospho & ~current,2], L)
            midpoint = (top + signed_delta(bottom, top, L)/2) % L
            head_z = float(head.positions[0,2])
            tail_z = float(tail.positions[0,2])
            records.append({
                "resindex": int(lipid["resindex"]), "species": str(lipid["species"]),
                "frame": frame, "time_ps": float(ts.time),
                "head_z_A": head_z, "tail_z_A": tail_z,
                "head_minus_tail_z_A": signed_delta(head_z, tail_z, L),
                "signed_head_midplane_A": signed_delta(head_z, midpoint, L),
                "upper_phosphate_center_z_A": top,
                "lower_phosphate_center_z_A": bottom,
                "midplane_z_A": midpoint,
            })
    data = pd.DataFrame(records)
    data.to_parquet(out / "changed_lipid_orientation.parquet", index=False)
    summary = {}
    for key, sub in data.groupby("resindex"):
        signed = sub.signed_head_midplane_A.to_numpy(float)
        orient = sub.head_minus_tail_z_A.to_numpy(float)
        summary[str(int(key))] = {
            "head_midplane_crossings": int(np.count_nonzero(np.sign(signed[1:]) != np.sign(signed[:-1]))),
            "orientation_sign_changes": int(np.count_nonzero(np.sign(orient[1:]) != np.sign(orient[:-1]))),
            "head_midplane_min_A": float(signed.min()),
            "head_midplane_max_A": float(signed.max()),
            "head_minus_tail_min_A": float(orient.min()),
            "head_minus_tail_max_A": float(orient.max()),
        }
    return summary


def hbond_audit(system, topology, trajectory):
    u = mda.Universe(str(topology), str(trajectory))
    lipid_sel = "chainID GM1a BGAL" if system == "carb_sdc4" else "chainID CER160"
    selected = u.select_atoms(lipid_sel)
    assert len(selected) > 0, f"{system}: no selected glycolipid/ceramide atoms"
    start = time.perf_counter()
    molecule_ids, species = lipid_molecule_map(u, selected)
    molecule_count = int(len(np.unique(molecule_ids[selected.indices])))
    assert molecule_count == 82, f"{system}: expected 82 glycolipid/ceramide molecules; found {molecule_count}"
    molecule_species = Counter()
    for mol in np.unique(molecule_ids[selected.indices]):
        idx = selected.indices[molecule_ids[selected.indices] == mol]
        labels = set(species[idx])
        assert len(labels) == 1, f"Mixed species in bonded molecule: {labels}"
        molecule_species[next(iter(labels))] += 1
    assert (molecule_species == Counter({"GM1a":24, "BGAL":58}) if system == "carb_sdc4"
            else molecule_species == Counter({"CER160":82}))
    map_s = time.perf_counter()-start
    selected_count = len(selected)
    del selected, u
    gc.collect()
    cfg = {
        "trajectory": {"start": 0, "stop": 20, "stride": 1},
        "time": {"dt_ps": None}, "progress": {"enabled": False},
        "analyses": {"lipid_hydrogen_bonds": {
            "lipid_selection": lipid_sel,
            "exclude_same_residue": True, "exclude_same_molecule": True,
            "update_selections": False,
            "d_h_cutoff_A": 1.2, "d_a_cutoff_A": 3.0,
            "angle_cutoff_deg": 150.,
        }},
    }
    rep = Replicate(system, 1, topology, trajectory, "rep1")
    out = ROOT / system / "physics"
    start = time.perf_counter()
    result = run_lipid_hydrogen_bonds(rep, cfg, out)
    elapsed = time.perf_counter()-start
    root = out / "interactions" / "lipid_hydrogen_bonds"
    audit = pd.read_parquet(root / "molecule_audit.parquet").iloc[0].to_dict()
    events = pd.read_parquet(root / "events.parquet")
    assert len(events) == result["events"]
    assert not (events.donor_molecule_id == events.acceptor_molecule_id).any()
    if len(events):
        assert set(events.donor_species).union(events.acceptor_species).issubset(
            set(molecule_species))
        assert events.distance_A.max() <= 3.0001
        assert events.angle_deg.min() >= 149.999
    report = {
        "n_frames": 20, "selected_atoms": selected_count,
        "n_molecules": molecule_count, "species": dict(molecule_species),
        "molecule_mapping_s": map_s, "hbond_20_frames_s": elapsed,
        "events": result["events"], "audit": audit,
        "rss_peak_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.,
    }
    print(f"{system}: glycolipid H-bond audit = {json.dumps(report, sort_keys=True)}",
          flush=True)
    return report


def main():
    summary = {}
    failed = False
    for system, (tpr, xtc) in SYSTEMS.items():
        summary[system] = {}
        try:
            summary[system]["geometry"] = geometry_audit(
                system, DATA / tpr, DATA / xtc)
            summary[system]["glycolipid_hbonds"] = hbond_audit(
                system, DATA / tpr, DATA / xtc)
            summary[system]["status"] = "passed"
        except Exception as exc:
            traceback.print_exc()
            summary[system]["status"] = "failed"
            summary[system]["error"] = repr(exc)
            failed = True
        finally:
            gc.collect()
    (ROOT / "physical_audit.json").write_text(json.dumps(summary, indent=2))
    print("PHYSICAL_AUDIT_JSON="+json.dumps(summary, sort_keys=True), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
