"""Audit SDC4 potassium coordination from one CHARMM TPR and optional XTC.

Outputs are per input trajectory. Keep conditions and replicate IDs separate.
Distances and topology masses are in Å and Da, respectively.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import re

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.spatial import cKDTree
import pyarrow.parquet as pq
import pyarrow as pa

from ..config import atomic_json, atomic_parquet, time_ps
from .hs_topology import discover_covalent_hs
from .mechanistic_common import sequence_matches
from ..observables.solvation import SPATIAL_OFFSETS, rdf_from_counts, water_donor_hbonds
from ..pipeline.progress import track

CONTACT_CUTOFFS = np.array([3., 3.5, 4., 5.])
STERIC_CUTOFFS = np.array([0., 1.8, 2.1, 2.4])
EDGES = np.linspace(0, 8, 81)


def topology_sites(u, config):
    """Identify terminal anionic oxygen sites by bonded graph and TPR charges."""
    masses, charges = u.atoms.masses, u.atoms.charges
    oxygen = np.abs(masses - 16) < .6
    carbon = np.abs(masses - 12) < .6
    nitrogen = np.abs(masses - 14) < .6
    sulfur = np.abs(masses - 32) < .6
    neighbors = [[] for _ in u.atoms]
    for a, b in u.bonds.indices:
        neighbors[a].append(b)
        neighbors[b].append(a)

    ions = u.select_atoms("resname POT and name POT")
    if (not len(ions) or not np.allclose(charges[ions.indices], 1., atol=.05)
            or not np.allclose(masses[ions.indices], 39.1, atol=.6)
            or any(neighbors[i] for i in ions.indices)
            or len(ions.residues) != len(ions)):
        raise ValueError("Expected unbonded +1 CHARMM POT ions with potassium mass")
    solvation = config["analyses"]["solvation"]
    water = u.select_atoms(solvation["water_selection"])
    if not len(water) or not oxygen[water.indices].all():
        raise ValueError("Water selector must resolve to oxygen atoms")
    hydrogen_ids = [[i for i in neighbors[o] if .5 < masses[i] < 2.]
                    for o in water.indices]
    if any(len(ids) != 2 for ids in hydrogen_ids):
        raise ValueError("Expected two bonded hydrogens per selected water oxygen")
    water_h = u.atoms[np.asarray(hydrogen_ids).ravel()]

    motif = solvation["targets"]["motif"]["acidic_patch"]
    matches = sequence_matches(u.select_atoms(motif["selection"]),
                               {"acid": re.escape(motif["sequence"])})
    if len(matches) != 2:
        raise ValueError(f"Expected two exact acidic-patch sequences, found {list(matches)}")
    chains, _, attachments = discover_covalent_hs(
        u, config["analyses"]["sdc4_mechanisms"]["hs_auto"])
    records = []
    for owner, group in ([(name, u.residues[ids].atoms) for name, ids in matches.items()]
                         + list(chains.items())):
        acid = owner in matches
        for residue in group.residues:
            if acid and residue.resname not in ("ASP", "GLU"):
                continue
            n_carboxylates = 0
            for center in residue.atoms:
                bonded = np.asarray(neighbors[center.index], dtype=int)
                oxygens = bonded[oxygen[bonded]]
                terminals = [i for i in oxygens if len(neighbors[i]) == 1]
                if (carbon[center.index] and len(bonded) == 3 and len(terminals) == 2
                        and carbon[bonded[~oxygen[bonded]]].sum() == 1):
                    kind = "protein_carboxylate" if acid else "hs_carboxylate"
                    n_carboxylates += 1
                elif (not acid and sulfur[center.index] and len(bonded) == 4
                      and len(terminals) == 3 and
                      (len(oxygens) == 4 or
                       (len(oxygens) == 3 and nitrogen[bonded].sum() == 1))):
                    kind = "hs_sulfate" if len(oxygens) == 4 else "hs_sulfamate"
                else:
                    continue
                for i in terminals:
                    atom = u.atoms[i]
                    if charges[i] >= 0 or atom.resindex != residue.ix:
                        raise ValueError(f"Unexpected oxygen charge or residue at atom {i}")
                    records.append((owner, kind, residue.ix, residue.resid,
                                    residue.resname, center.name, i, atom.name, float(charges[i])))
            if acid and n_carboxylates != 1:
                raise ValueError(f"Expected one deprotonated side-chain carboxylate in {owner} residue {residue}")
    sites = pd.DataFrame(records, columns=["owner", "kind", "resindex", "resid",
                                          "resname", "center_name", "atom_index", "atom_name", "charge_e"])
    if (sites.empty or set(sites.owner) != set(matches) | set(chains)
            or "hs_carboxylate" not in set(sites.kind)
            or not {"hs_sulfate", "hs_sulfamate"} & set(sites.kind)):
        raise ValueError("Incomplete protein/HS anionic site inventory")
    occluders = u.atoms[(masses > 2) & ~np.isin(u.atoms.resindices, water.resindices)
                         & (np.bincount(u.atoms.resindices)[u.atoms.resindices] > 1)]
    return sites, ions, water, water_h, occluders, attachments


def bridge_rows(positions, targets, occluder_tree, box):
    """Paired observed/shifted counts at identical cutoffs and topology geometry."""
    separations = {(a, b): targets[a].query(targets[b].data)[0].min()
                   for a in range(2) for b in range(2, len(targets))}
    for du, dv in ((0., 0.), *SPATIAL_OFFSETS):
        shifted = (positions + box[:3] * [du, dv, 0.]) % box[:3]
        nearest, _ = occluder_tree.query(shifted)
        proximity = np.array([tree.query(shifted)[0] for tree in targets])
        eligible = nearest[:, None] >= STERIC_CUTOFFS
        eligible_count = eligible.sum(axis=0)
        for (a, b), separation in separations.items():
            bridges = np.maximum(proximity[a], proximity[b])[:, None] <= CONTACT_CUTOFFS
            counts = bridges.astype(int).T @ eligible
            for i, cutoff in enumerate(CONTACT_CUTOFFS):
                for j, steric in enumerate(STERIC_CUTOFFS):
                    yield (a, b, du, dv, cutoff, steric, separation, eligible_count[j], counts[i, j])


def peak_rows(rdf, frames):
    """Locate descriptive peaks and intervening minimum with two smoothing widths."""
    r = (EDGES[:-1] + EDGES[1:]) / 2
    for kind, group in rdf.groupby("kind", sort=False):
        for sigma in (1., 2.):
            g = gaussian_filter1d(group.g_box.to_numpy(), sigma)
            if len(g) != len(r) or not np.isfinite(g).all() or not np.any(g > 0):
                continue
            peaks, _ = find_peaks(g, prominence=.05 * g.max())
            peaks = peaks[(r[peaks] >= 2.) & (r[peaks] <= 7.)]
            first, second = list(peaks[:2]) + [None] * (2 - len(peaks[:2]))
            minimum = (first + np.argmin(g[first:second+1])) if second is not None else None
            yield (kind, frames, sigma * .1, *[float(r[i]) if i is not None else np.nan
                                                for i in (first, minimum, second)])


def run_audit(config, tpr, xtc, output, *, start=None, stop=None, stride=1,
              initial_frame=False):
    """Write one replicate's topology and optional sampled-frame audit."""
    u = mda.Universe(str(tpr), *([str(xtc)] if xtc else []))
    sites, ions, water, water_h, occluders, attachments = topology_sites(u, config)
    output.mkdir(parents=True, exist_ok=False)
    atomic_parquet(output / "chemical_sites.parquet", sites)
    atomic_json(output / "topology_audit.json", {
        "tpr": str(tpr), "xtc": str(xtc) if xtc else None,
        "initial_frame": initial_frame, "potassium_ions": len(ions), "waters": len(water), "hs_attachments": attachments,
        "start": start, "stop": stop, "stride": stride, "time": config.get("time", {}),
        "site_counts": sites.groupby("kind").size().to_dict(),
    })
    if not (xtc or initial_frame):
        return

    groups = [((owner, kind), u.atoms[group.atom_index.to_numpy()])
              for (kind, owner), group in sites.groupby(["kind", "owner"], sort=False)]
    groups.sort(key=lambda item: item[0][1] != "protein_carboxylate")
    kinds = sorted(set(sites.kind)) + ["water_oxygen"]
    counts = {kind: np.zeros(len(EDGES)-1, int) for kind in kinds}
    shared = {kind: np.zeros(len(EDGES)-1, int) for kind in kinds if kind != "water_oxygen"}
    exposure = dict.fromkeys(kinds, 0.)
    keys = ["acid_patch", "hs_chain", "hs_kind", "contact_cutoff_A", "steric_cutoff_A", "reference"]
    columns = ["frame", "time_ps", "acid_patch", "hs_chain", "hs_kind", "shift_x", "shift_y",
               "contact_cutoff_A", "steric_cutoff_A", "min_target_separation_A",
               "eligible_ions", "bridging_ions"]
    totals = None
    frames = 0
    temporary = output / ".bridge_frames.tmp"
    writer = None
    try:
        window = u.trajectory[start:stop:stride]
        for ts in track(window, total=len(window), desc="Potassium bridge audit",
                        enabled=config.get("progress", {}).get("enabled", True)):
            box = ts.dimensions
            if (box is None or not np.isfinite(box).all()
                    or not np.allclose(box[3:], 90., atol=1e-5, rtol=0)
                    or np.min(box[:3]) < 2 * EDGES[-1]):
                raise ValueError(f"Frame {ts.frame} needs an orthorhombic box and an inscribed RDF sphere")
            k, w = ions.positions, water.positions
            density = len(k) / np.prod(box[:3])
            pairs, distances = capped_distance(k, w, EDGES[-1], box=box)
            counts["water_oxygen"] += np.histogram(distances, EDGES)[0]
            exposure["water_oxygen"] += len(w) * density
            hydrated = {}
            for ion, oxygen in pairs[distances <= 3.5]:
                hydrated.setdefault(ion, set()).add(oxygen)
            trees = []
            for (_, kind), atoms in groups:
                xyz = atoms.positions
                pair, dist = capped_distance(xyz, k, EDGES[-1], box=box)
                counts[kind] += np.histogram(dist, EDGES)[0]
                exposure[kind] += len(xyz) * density
                candidates = capped_distance(xyz, w, 3.5, box=box, return_distances=False)
                bonds = water_donor_hbonds(xyz, w, water_h.positions.reshape(-1, 2, 3),
                                          candidates, box, max_da_A=3.5, min_angle_deg=150.)
                donor = {}
                for water_id, site, _ in bonds:
                    donor.setdefault(site, set()).add(water_id)
                mediated = np.fromiter((bool(donor.get(site, set()) & hydrated.get(ion, set()))
                                        for site, ion in pair), bool, count=len(pair))
                shared[kind] += np.histogram(dist[mediated], EDGES)[0]
                trees.append(cKDTree(xyz % box[:3], boxsize=box[:3]))
            blockers = cKDTree(occluders.positions % box[:3], boxsize=box[:3])
            timestamp = np.nan if initial_frame else time_ps(ts, frames, config)
            records = ((ts.frame, timestamp, groups[a][0][0], *groups[b][0], *values)
                       for a, b, *values in bridge_rows(k, trees, blockers, box))
            bridge = pd.DataFrame.from_records(records, columns=columns)
            table = pa.Table.from_pandas(bridge, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(temporary, table.schema)
            writer.write_table(table)
            by_reference = bridge.assign(
                reference=np.where(bridge.shift_x.eq(0) & bridge.shift_y.eq(0), "observed", "shifted"),
                bridge_present=bridge.bridging_ions.gt(0), samples=1)
            frame_totals = by_reference.groupby(keys)[
                ["bridge_present", "eligible_ions", "bridging_ions", "samples"]].sum()
            totals = frame_totals if totals is None else totals + frame_totals
            frames += 1
        if not frames:
            raise ValueError("The trajectory frame window is empty")
        writer.close()
        writer = None
        temporary.replace(output / "bridge_frames.parquet")
    finally:
        if writer is not None:
            writer.close()
        temporary.unlink(missing_ok=True)
    summary = (totals.drop(columns="samples").div(totals.samples, axis=0)
               .rename(columns={"bridge_present": "occupancy"}).reset_index())
    summary["bridge_fraction_of_eligible_ion_samples"] = (
        summary.bridging_ions.div(summary.eligible_ions.where(summary.eligible_ions > 0)))
    atomic_parquet(output / "bridge_summary.parquet", summary)
    rdf = pd.concat([pd.DataFrame({"kind": kind, "r_inner_A": EDGES[:-1], "r_outer_A": EDGES[1:],
                                  "pair_count": counts[kind], "g_box": rdf_from_counts(counts[kind], EDGES,
                                                                                      exposure[kind]),
                                  "shared_water_pair_count": shared.get(kind, np.nan)})
                     for kind in kinds], ignore_index=True)
    atomic_parquet(output / "rdf.parquet", rdf)
    atomic_parquet(output / "rdf_peak_candidates.parquet", pd.DataFrame(
        peak_rows(rdf, frames), columns=["kind", "frames", "smoothing_A", "first_peak_A",
                                        "intervening_minimum_A", "second_peak_A"]))


def run_potassium_bridge_audit(rep, cfg, repout):
    output = repout / "extensions" / "potassium_bridge_audit"
    output.parent.mkdir(parents=True, exist_ok=True)
    window = cfg["trajectory"]
    with TemporaryDirectory(prefix=".potassium-audit-", dir=output.parent) as tmp:
        staging = Path(tmp) / "audit"
        run_audit(cfg, rep.topology, rep.trajectory, staging,
                  start=window["start"], stop=window["stop"], stride=window["stride"])
        staging.rename(output)
