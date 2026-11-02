"""SDC4 molecular-mechanism observers sharing the raw streaming trajectory pass.

All selections are explicit topology selections. No residue numbering, HS chain
identity, or GM1/BGAL molecular connectivity is guessed. Geometry is PBC-aware
and raw (never RMSD-aligned); distances are Å and time is ps.

Density maps are counts/area (not solvent-accessible-area-corrected), NxIP
steric occlusion is a proxy rather than binding affinity, and local curvature
is a small-slope quadratic height-field estimate.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import json
import re
import time

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance, minimize_vectors, calc_dihedrals

from ..config import DEFAULTS, atomic_json, atomic_parquet, ensure_dir, time_ps
from ..core.lipid_geometry import anchor_groups, aligned_leaflets, periodic_center
from .mechanistic_common import AA3, heavy, nearest, selected, sequence_matches, unit
from ..observables.membrane_order import _candidate_sites, _fixed_ch_pairs, _atomic_npz
from ..pipeline.progress import track
from .hs_topology import discover_covalent_hs


def enabled(cfg):
    return bool(cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("enabled", False))


def _id(label):
    value = re.sub(r"[^A-Za-z0-9_]+", "_", str(label)).strip("_")
    if not value:
        raise ValueError("Empty/invalid mechanism group label")
    return value


def _box(u):
    box = np.asarray(u.dimensions, float)
    if box.shape != (6,) or not np.all(np.isfinite(box)) or np.any(box[:3] <= 0):
        raise ValueError("SDC4 mechanisms require a finite periodic box")
    if not np.allclose(box[3:], 90., atol=.05):
        raise ValueError("SDC4 xy maps and local height fields require an orthorhombic box")
    return box


def _vector(a, b, box):
    return minimize_vectors(np.asarray(b, float) - np.asarray(a, float), box)


def _center(atoms, box):
    if not len(atoms):
        raise ValueError("Cannot calculate a center from an empty selection")
    return periodic_center(atoms.positions, box)


def _angle(a, b):
    a, b = unit(a), unit(b)
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b), -1., 1.)))) if (
        np.all(np.isfinite(a)) and np.all(np.isfinite(b))
    ) else np.nan


def _contacts(a, b, cutoff, box):
    if not len(a) or not len(b):
        return np.empty((0, 2), dtype=int)
    return capped_distance(
        a.positions, b.positions, max_cutoff=float(cutoff),
        box=box, return_distances=False,
    )


def _zmean(z, length):
    theta = np.asarray(z, float) * (2. * np.pi / length)
    return (np.arctan2(np.sin(theta).mean(), np.cos(theta).mean())
            * length / (2. * np.pi)) % length


def _local_surface(xy, z, minimum=8):
    """Quadratic z=a*x²+b*y²+c*xy+d*x+e*y+f, local small-slope H≈a+b."""
    if len(xy) < minimum:
        return np.nan
    x, y = xy[:, 0], xy[:, 1]
    scale = max(np.sqrt(np.mean(x*x+y*y)), 1.)
    X, Y = x / scale, y / scale
    design = np.column_stack((X*X, Y*Y, X*Y, X, Y, np.ones(len(X))))
    if np.linalg.matrix_rank(design) < 6 or np.linalg.cond(design) > 1e7:
        return np.nan
    coefficients = np.linalg.lstsq(design, z, rcond=None)[0]
    return float((coefficients[0] + coefficients[1]) / (scale * scale))


def _backbone_axis(atoms, box):
    """PBC-unwrapped PCA axis, oriented along the selected residue sequence."""
    xyz = atoms.positions.astype(float)
    if len(xyz) < 3:
        return np.full(3, np.nan)
    xyz = xyz[0] + minimize_vectors(xyz - xyz[0], box)
    _, _, vt = np.linalg.svd(xyz - xyz.mean(axis=0), full_matrices=False)
    direction = vt[0]
    if np.dot(direction, xyz[-1] - xyz[0]) < 0:
        direction = -direction
    return unit(direction)


def _normal_order_setup(u, spec):
    if not spec.get("enabled", False):
        return None
    options = dict(DEFAULTS["analyses"]["membrane_order"])
    options.update(spec)
    if not options.get("lipid_selection"):
        raise ValueError("sdc4_mechanisms.local_order.lipid_selection is required")
    sites = _candidate_sites(u, options)
    if not sites:
        raise ValueError("No acyl carbon sites matched local_order selection")
    sites, carbons, hydrogens, pairs = _fixed_ch_pairs(
        u, sites, float(options.get("hydrogen_cutoff_A", 1.25))
    )
    return {"sites": sites, "carbons": carbons, "hydrogens": hydrogens,
            "pairs": pairs, "resindices": np.array([s.residue_index for s in sites])}


def _order_by_residue(o, box, axis):
    pairs = o["pairs"]
    vec = minimize_vectors(
        o["hydrogens"].positions[pairs[:, 1]]
        - o["carbons"].positions[pairs[:, 0]], box,
    )
    norms = np.linalg.norm(vec, axis=1)
    valid = norms > 1e-8
    value = np.full(len(pairs), np.nan)
    value[valid] = -0.5 * (3 * (vec[valid, axis] / norms[valid])**2 - 1)
    sums = np.bincount(pairs[:, 0][valid], weights=value[valid],
                       minlength=len(o["sites"]))
    counts = np.bincount(pairs[:, 0][valid], minlength=len(o["sites"]))
    site = np.divide(sums, counts, out=np.full(len(o["sites"]), np.nan),
                     where=counts > 0)
    # Lipids with more tail sites are *not* allowed to dominate the mean.
    indices, inverse = np.unique(o["resindices"], return_inverse=True)
    good = np.isfinite(site)
    residue_sum = np.bincount(inverse[good], weights=site[good], minlength=len(indices))
    residue_count = np.bincount(inverse[good], minlength=len(indices))
    return indices, np.divide(
        residue_sum, residue_count, out=np.full(len(indices), np.nan),
        where=residue_count > 0,
    )


class SDC4Mechanisms:
    def __init__(self, universe, cfg, repout):
        self.u, self.cfg, self.repout = universe, cfg, Path(repout)
        self.spec = cfg["analyses"]["sdc4_mechanisms"]
        s = self.spec
        self.axis = {"x": 0, "y": 1, "z": 2}[str(s.get("axis", "z")).lower()]
        if self.axis != 2:
            raise ValueError("sdc4_mechanisms currently supports a z-normal membrane only")
        self.box = _box(universe)
        self.cache = ensure_dir(self.repout / "cache" / "sdc4_mechanisms")
        self.rows = defaultdict(list)
        self.frame_count = 0
        self.timings = defaultdict(float)
        self.tm = { _id(k): selected(universe, v, required=True)
                    for k, v in (s.get("tm_chains") or {}).items() }
        self.ecto = { _id(k): selected(universe, v, required=True)
                      for k, v in (s.get("ectodomains") or {}).items() }
        if len(self.tm) != 2 or set(self.tm) != set(self.ecto):
            raise ValueError("Supply exactly two matching tm_chains and ectodomains (A/B)")
        self.chain_names = tuple(self.tm)
        self.tm_heavy = {
            k: heavy(universe.atoms[np.isin(universe.atoms.resindices,
                                            v.resindices)])
            for k, v in self.tm.items()
        }
        self.hs_auto_audit = {}
        if (s.get("hs_auto") or {}).get("enabled", False):
            self.hs, self.attach, self.hs_auto_audit = discover_covalent_hs(
                universe, s["hs_auto"]
            )
            self.hs = {k: heavy(v) for k, v in self.hs.items()}
        else:
            self.hs = { _id(k): heavy(selected(universe, v, required=True))
                        for k, v in (s.get("hs_chains") or {}).items() }
            self.attach = { _id(k): selected(universe, v, required=True)
                           for k, v in (s.get("hs_attachment_sites") or {}).items() }
        if set(self.hs) != set(self.attach):
            raise ValueError("Each HS chain requires a matching hs_attachment_sites selection")
        self.lipid_spec = {
            "axis": "z", "species_selections": s.get("lipid_anchors") or {},
            "allow_absent_species": bool(s.get("allow_absent_species", False)),
        }
        self.species, self.resindices, self.tracker, self.initial_leaf = anchor_groups(
            universe, self.lipid_spec
        )
        # Thickness and curvature must use chemically comparable surface
        # anchors. Mixing phospholipid P, cholesterol O3 and ceramide NF
        # creates apparent height changes when species composition shifts.
        self.surface_species = tuple(s.get("surface_anchor_species") or (
            "POPC", "POPE", "POPI", "POPS", "DOPC", "DPPC", "PSM"
        ))
        self.surface_mask = np.isin(self.species, self.surface_species)
        if np.count_nonzero(self.surface_mask & self.initial_leaf) < 4 or (
            np.count_nonzero(self.surface_mask & ~self.initial_leaf) < 4
        ):
            raise ValueError("Need at least four comparable surface anchors per leaflet; check surface_anchor_species")
        self.lipid_contact = heavy(selected(
            universe, s.get("lipid_contact_selection"), required=True
        ))
        self.motif_lipid_residence_enabled = bool(s.get("motif_lipid_residence", False))
        self.lipid_molecule_by_atom = {}
        self.lipid_molecule_species = {}
        self.unmapped_lipid_fragment_count = 0
        if self.motif_lipid_residence_enabled:
            # Bonded fragments are stable physical molecules in a TPR. Unlike
            # residue numbers they correctly group GM1/BGAL ceramide + sugars.
            try:
                fragments = self.lipid_contact.fragments
            except Exception as exc:
                raise ValueError(
                    "Motif lipid residence requires lipid bond connectivity "
                    "to identify complete molecules"
                ) from exc
            anchor_by_residue = dict(zip(
                map(int, self.resindices), map(str, self.species)
            ))
            for fragment in fragments:
                # For glycolipids, chainID is the molecule species label;
                # the CER160 residue is only their ceramide anchor. Without
                # glycans, chainID CER160 denotes a distinct comparison group.
                chain_labels = set(map(str, fragment.chainIDs))
                glyco_labels = chain_labels & {"GM1a", "BGAL", "CER160"}
                if len(glyco_labels) > 1:
                    raise ValueError("A lipid molecule spans conflicting glyco chainIDs")
                anchors = {anchor_by_residue[int(ix)]
                           for ix in fragment.resindices
                           if int(ix) in anchor_by_residue}
                species = set(glyco_labels) if glyco_labels else anchors
                if glyco_labels and anchors and anchors != species:
                    raise ValueError("Glycolipid chainID conflicts with lipid anchor species")
                if not species:
                    self.unmapped_lipid_fragment_count += 1
                    continue
                if len(species) > 1:
                    raise ValueError("Lipid fragment has conflicting species anchors")
                identifier = int(np.min(fragment.indices))
                self.lipid_molecule_species[identifier] = species.pop()
                for index in fragment.indices:
                    self.lipid_molecule_by_atom[int(index)] = identifier
            if not self.lipid_molecule_by_atom:
                raise ValueError("No lipid molecules identified for motif residence")
        self.lipid_anchor_index_by_molecule = {}
        if self.motif_lipid_residence_enabled:
            anchor_index_by_residue = {
                int(res): i for i, res in enumerate(self.resindices)
            }
            for atom_index, molecule in self.lipid_molecule_by_atom.items():
                residue = int(universe.atoms[atom_index].resindex)
                if residue in anchor_index_by_residue:
                    anchor_index = anchor_index_by_residue[residue]
                    previous = self.lipid_anchor_index_by_molecule.get(molecule)
                    if previous is not None and previous != anchor_index:
                        raise ValueError("Lipid molecule has multiple distinct leaflet anchors")
                    self.lipid_anchor_index_by_molecule[molecule] = anchor_index
        self.protein = heavy(selected(
            universe, s.get("protein_selection") or cfg["selections"]["protein"],
            required=True,
        ))
        self.water = selected(universe, s.get("water_oxygen_selection"), required=bool(self.hs))
        self.ions = { _id(k): selected(universe, v, required=True)
                      for k, v in (s.get("ion_selections") or {}).items() }
        self.motifs = {}
        ca = selected(universe, s.get("motif_selection") or
                      (cfg["selections"]["protein"] + " and name CA"), required=True)
        for key, resindices in sequence_matches(
            ca, s.get("motifs", {"NxIP": "N.IP"})
        ).items():
            atoms = heavy(universe.atoms[np.isin(universe.atoms.resindices, resindices)])
            if len(atoms):
                self.motifs[_id(key)] = {"atoms": atoms, "resindices": resindices}
        # Prepare backbone torsions for the two prolines in each mapped
        # VVHPLVPL motif. Actual residue indices come from the topology, and
        # missing backbone atoms remain NaN rather than being invented.
        self.proline_torsions = {}
        for name, info in self.motifs.items():
            if not name.startswith("proline_hydrophobic_"):
                continue
            ids = [int(x) for x in info["resindices"]]
            for j in (3, 6):  # positions 4 and 7 of VVHPLVPL
                if len(ids) != 8 or self.u.residues[ids[j]].resname != "PRO":
                    raise ValueError(f"Unexpected proline mapping for {name}")
                before, current, after = (
                    self.u.residues[ids[j-1]], self.u.residues[ids[j]],
                    self.u.residues[ids[j+1]],
                )
                def atom(residue, label):
                    matches = residue.atoms.select_atoms(f"name {label}")
                    return matches[0] if len(matches) == 1 else None
                self.proline_torsions[f"{name}_P{j+1}"] = (
                    (atom(before, "C"), atom(current, "N"),
                     atom(current, "CA"), atom(current, "C")),
                    (atom(current, "N"), atom(current, "CA"),
                     atom(current, "C"), atom(after, "N")),
                )
        # Fail explicitly if the functional endpoint was requested but the
        # topology's actual sequence does not contain the motif.
        if s.get("require_motifs", True) and not self.motifs:
            raise ValueError("No NxIP motif matched; verify the mature SDC4 sequence/chain selections")
        self.gxxxg = {}
        self.gxxxg_sequences = {}
        for chain, tm in self.tm.items():
            ca_tm = tm.select_atoms("name CA")
            if not len(ca_tm):
                raise ValueError(f"TM selection for {chain} must contain CA atoms")
            # Match the known SDC4 GGIVG sequence; generic G...G matches can
            # select an unrelated glycine pair and mask a wrong TM selection.
            matches = sequence_matches(ca_tm, {"GGIVG": "GGIVG"})
            indices = np.unique(np.concatenate(list(matches.values()))) if matches else np.empty(0, int)
            sequence = "".join(
                AA3.get(atom.resname.upper(), "X") for atom in ca_tm
            )
            self.gxxxg_sequences[chain] = sequence
            if not len(indices) and s.get("require_gxxxg", False):
                raise ValueError(
                    f"No GGIVG found in TM chain {chain} (sequence={sequence!r}); "
                    "verify topology chain IDs and TM residue range"
                )
            self.gxxxg[chain] = heavy(universe.atoms[np.isin(universe.atoms.resindices, indices)])
        self.local_order = _normal_order_setup(universe, s.get("local_order") or {})
        self.extent = float(s.get("map_extent_A", 40.))
        self.bin_A = float(s.get("map_bin_A", 2.))
        n = int(np.ceil(2 * self.extent / self.bin_A))
        if n < 2 or n > 400:
            raise ValueError("Invalid lipid map resolution (2 <= bins <= 400)")
        self.edges = np.linspace(-self.extent, self.extent, n + 1)
        self.map_names = [(species, leaf) for species in sorted(set(self.species))
                          for leaf in ("upper", "lower")]
        self.map_counts = np.zeros((len(self.map_names), n, n), np.int64)
        self.map_baselines = np.zeros(len(self.map_names), float)
        self.map_frames = 0
        self.contact_A = float(s.get("contact_cutoff_A", 5.))
        self.hydration_A = float(s.get("hydration_cutoff_A", 3.5))
        self.ion_A = float(s.get("ion_cutoff_A", 4.))
        self.local_A = float(s.get("local_radius_A", 20.))
        self.min_lipids = int(s.get("local_min_lipids", 8))
        if min(self.contact_A, self.hydration_A, self.ion_A, self.local_A) <= 0:
            raise ValueError("All SDC4 mechanism distance cutoffs must be positive")
        self.hs_distal = {}
        for key, ag in self.hs.items():
            # Determine distal terminus by proximity of end residues to the
            # attachment selection on the initial frame, never by atom index.
            residues = list(ag.residues)
            if len(residues) < 2:
                raise ValueError(f"HS {key} must contain at least two glycan residues")
            first = heavy(residues[0].atoms)
            last = heavy(residues[-1].atoms)
            dfirst = nearest(first, self.attach[key], self.box)
            dlast = nearest(last, self.attach[key], self.box)
            self.hs_distal[key] = first if dfirst > dlast else last
        self.audit = {
            "tm_atom_counts": {k: len(v) for k, v in self.tm.items()},
            "ectodomain_atom_counts": {k: len(v) for k, v in self.ecto.items()},
            "hs_atom_counts": {k: len(v) for k, v in self.hs.items()},
            "hs_covalent_attachment_audit": self.hs_auto_audit,
            "hs_distal_resindices": {k: v.resindices.tolist()
                                     for k, v in self.hs_distal.items()},
            "motif_resindices": {k: v["resindices"].tolist()
                                 for k, v in self.motifs.items()},
            "gxxxg_resindices": {k: v.residues.ix.tolist()
                                 for k, v in self.gxxxg.items()},
            "gxxxg_tm_sequences": self.gxxxg_sequences,
            "motif_lipid_residence_unmapped_fragments": self.unmapped_lipid_fragment_count,
            "motif_lipid_residence_molecules_by_species": {
                species: sum(v == species for v in self.lipid_molecule_species.values())
                for species in sorted(set(self.lipid_molecule_species.values()))
            },
            "surface_anchor_species": list(self.surface_species),
            "lipid_counts": {str(k): int(np.count_nonzero(self.species == k))
                             for k in sorted(set(self.species))},
            "map_rotation": "x axis TM chain 1 -> TM chain 2; y is +90 degree rotation",
            "map_density_units": "lipids/nm^2; no protein excluded-area correction",
            "local_curvature_units": "A^-1, small-slope quadratic approximation",
        }

    def update(self, ts, frame_index):
        start = time.perf_counter()
        box = _box(self.u)
        t = time_ps(ts, frame_index, self.cfg)
        anchors = self.tracker.positions()
        leaf, changed = aligned_leaflets(anchors[:, 2], box[2], self.initial_leaf)
        tm_centers = {k: _center(g, box) for k, g in self.tm.items()}
        first, second = self.chain_names
        midpoint = tm_centers[first] + .5 * _vector(
            tm_centers[first], tm_centers[second], box
        )
        axis_vec = _vector(tm_centers[first], tm_centers[second], box)[:2]
        axis_norm = np.linalg.norm(axis_vec)
        if axis_norm < 1e-6:
            raise ValueError("TM centers coincide laterally; cannot orient lipid map")
        xaxis = axis_vec / axis_norm
        yaxis = np.array([-xaxis[1], xaxis[0]])
        dxy = minimize_vectors(anchors - midpoint, box)[:, :2]
        rotated = np.column_stack((dxy @ xaxis, dxy @ yaxis))
        area_A2 = box[0] * box[1]
        for j, (species, leaflet) in enumerate(self.map_names):
            mask = (self.species == species) & (leaf if leaflet == "upper" else ~leaf)
            h, _, _ = np.histogram2d(
                rotated[mask, 0], rotated[mask, 1],
                bins=(self.edges, self.edges),
            )
            self.map_counts[j] += h.astype(np.int64)
            self.map_baselines[j] += int(mask.sum()) / area_A2
        self.map_frames += 1
        row = {"time_ps": t, "sdc4_lipid_leaflet_changes": changed}
        upper_z = _zmean(anchors[leaf & self.surface_mask, 2], box[2])
        lower_z = _zmean(anchors[~leaf & self.surface_mask, 2], box[2])
        mid_z = lower_z + .5 * _vector(
            [0., 0., lower_z], [0., 0., upper_z], box
        )[2]
        tm_z = float(midpoint[2])
        local = np.linalg.norm(dxy, axis=1) <= self.local_A
        upper_local = local & leaf & self.surface_mask
        lower_local = local & ~leaf & self.surface_mask
        for label, mask in (("upper", upper_local), ("lower", lower_local)):
            row[f"sdc4_local_{label}_n_lipids"] = int(mask.sum())
            if np.count_nonzero(mask) >= self.min_lipids:
                zz = minimize_vectors(
                    anchors[mask] - np.array([midpoint[0], midpoint[1], mid_z]), box
                )[:, 2]
                row[f"sdc4_local_{label}_height_A"] = float(np.mean(zz))
                row[f"sdc4_local_{label}_curvature_Ainv"] = _local_surface(
                    dxy[mask], zz, minimum=max(8, self.min_lipids)
                )
            else:
                row[f"sdc4_local_{label}_height_A"] = np.nan
                row[f"sdc4_local_{label}_curvature_Ainv"] = np.nan
        a, b = row["sdc4_local_upper_height_A"], row["sdc4_local_lower_height_A"]
        row["sdc4_local_thickness_A"] = float(a-b) if np.isfinite(a) and np.isfinite(b) else np.nan
        if self.local_order is not None:
            indices, values = _order_by_residue(self.local_order, box, self.axis)
            local_res = set(self.resindices[local].tolist())
            keep = np.isin(indices, list(local_res)) & np.isfinite(values)
            row["sdc4_local_minus_SCH"] = float(values[keep].mean()) if np.any(keep) else np.nan
            row["sdc4_bulk_minus_SCH"] = float(np.nanmean(values)) if np.any(np.isfinite(values)) else np.nan
        for chain, ag in self.ecto.items():
            ecto_center = _center(ag, box)
            # Determine the extracellular side from the distal ectodomain's
            # signed distance to the bilayer midplane, not from z alone.
            sign = np.sign(_vector(
                [midpoint[0], midpoint[1], mid_z], ecto_center, box
            )[2]) or 1.
            leaflet_z = upper_z if sign > 0 else lower_z
            height = _vector(
                [midpoint[0], midpoint[1], leaflet_z], ecto_center, box
            )[2] * sign
            ca = ag.select_atoms("name CA")
            if len(ca) < 4:
                raise ValueError(f"Ectodomain {chain} needs at least 4 CA atoms")
            # N-terminal end is distal in mature SDC4; use end-to-end axis
            # with first/last 10% CA to reduce single-atom noise.
            nend = max(2, len(ca) // 10)
            distal = _center(ca[:nend], box)
            proximal = _center(ca[-nend:], box)
            direction = _vector(proximal, distal, box)
            tilt = _angle(direction, [0., 0., sign])
            hits = _contacts(heavy(ag), self.lipid_contact, self.contact_A, box)
            prefix = f"sdc4_ectodomain_{chain}_"
            row[prefix + "height_A"] = float(height)
            row[prefix + "tilt_deg"] = tilt
            row[prefix + "membrane_contact_state"] = int(len(hits) > 0)
            row[prefix + "lipid_contact_atom_fraction"] = (
                len(np.unique(hits[:, 0])) / len(heavy(ag)) if len(hits) else 0.
            )
            row[prefix + "lipid_contact_residue_count"] = (
                len(np.unique(self.lipid_contact.resindices[hits[:, 1]])) if len(hits) else 0
            )
        axisA = _backbone_axis(self.tm[first].select_atoms("name CA"), box)
        axisB = _backbone_axis(self.tm[second].select_atoms("name CA"), box)
        row["sdc4_tm_crossing_angle_deg"] = _angle(axisA, axisB)
        cross = np.cross(axisA, axisB)
        separation = unit(_vector(tm_centers[first], tm_centers[second], box))
        row["sdc4_tm_crossing_handedness"] = (
            float(np.sign(np.dot(cross, separation)))
            if np.all(np.isfinite(cross)) else np.nan
        )
        tm1, tm2 = self.tm_heavy[first], self.tm_heavy[second]
        interface_hits = _contacts(tm1, tm2, self.contact_A, box)
        row["sdc4_tm_interface_residue_pairs"] = (
            len(np.unique(np.column_stack((
                tm1.resindices[interface_hits[:, 0]],
                tm2.resindices[interface_hits[:, 1]],
            )), axis=0)) if len(interface_hits) else 0
        )
        g1, g2 = self.gxxxg[first], self.gxxxg[second]
        for name, motif, opposite in (
            (first, g1, tm2), (second, g2, tm1)
        ):
            motif_hits = _contacts(motif, opposite, self.contact_A, box)
            row[f"sdc4_gxxxg_{name}_to_other_tm_residue_pairs"] = (
                len(np.unique(np.column_stack((
                    motif.resindices[motif_hits[:, 0]],
                    opposite.resindices[motif_hits[:, 1]],
                )), axis=0)) if len(motif_hits) else 0
            )
        ghits = _contacts(g1, g2, self.contact_A, box)
        row["sdc4_gxxxg_interface_atom_contacts"] = len(ghits)
        row["sdc4_gxxxg_interface_residue_pairs"] = (
            len(np.unique(np.column_stack((
                g1.resindices[ghits[:, 0]], g2.resindices[ghits[:, 1]]
            )), axis=0)) if len(ghits) else 0
        )
        for label, (phi_atoms, psi_atoms) in self.proline_torsions.items():
            for suffix, atom_group in (("phi_deg", phi_atoms), ("psi_deg", psi_atoms)):
                value = np.nan
                if all(atom is not None for atom in atom_group):
                    points = [atom.position for atom in atom_group]
                    value = float(np.degrees(calc_dihedrals(
                        points[0], points[1], points[2], points[3], box=box
                    )))
                row[f"sdc4_proline_{label}_{suffix}"] = value
        for name, info in self.motifs.items():
            ag = info["atoms"]
            own = info["resindices"]
            # Avoid calling adjacent backbone atoms "occlusion".
            # Adjacent backbone residues are covalent/near-covalent neighbors,
            # not evidence that an NxIP motif is functionally occluded.
            near_res = set(
                int(j) for ix in own for j in range(max(0, int(ix)-2), int(ix)+3)
            )
            protein_other = self.protein[
                ~np.isin(self.protein.resindices, list(near_res))
            ]
            contacts = _contacts(ag, protein_other, self.contact_A, box)
            glycan_all = self.u.atoms[np.unique(np.concatenate([
                g.indices for g in self.hs.values()
            ]))] if self.hs else self.u.atoms[[]]
            gly_hits = _contacts(ag, glycan_all, self.contact_A, box)
            lip_hits = _contacts(ag, self.lipid_contact, self.contact_A, box)
            occluded = np.unique(np.concatenate([
                contacts[:, 0], gly_hits[:, 0], lip_hits[:, 0]
            ])) if (len(contacts) or len(gly_hits) or len(lip_hits)) else np.empty(0, int)
            prefix = f"sdc4_motif_{name}_"
            row[prefix + "steric_occlusion_fraction"] = float(len(occluded) / len(ag))
            row[prefix + "protein_neighbor_residues"] = (
                len(np.unique(protein_other.resindices[contacts[:, 1]])) if len(contacts) else 0
            )
            row[prefix + "glycan_contacts"] = len(gly_hits)
            row[prefix + "lipid_contacts"] = len(lip_hits)
            if self.motif_lipid_residence_enabled and len(lip_hits):
                molecules = {
                    self.lipid_molecule_by_atom[int(self.lipid_contact.indices[j])]
                    for j in lip_hits[:, 1]
                    if int(self.lipid_contact.indices[j]) in self.lipid_molecule_by_atom
                }
                for molecule in sorted(molecules):
                    self.rows["lipid_residence"].append({
                        "time_ps": t, "motif": name,
                        "molecule_id": molecule,
                        "species": self.lipid_molecule_species[molecule],
                    })
        if self.motif_lipid_residence_enabled:
            # Species-specific molecule occupancy, with leaflet denominator
            # from the *same frame*. Absent species use NaN, never zero.
            regions = {
                **{f"TM_{chain}": atoms for chain, atoms in self.tm_heavy.items()},
                **{f"motif_{name}": item["atoms"] for name, item in self.motifs.items()},
            }
            availability = {}
            for species in sorted(set(self.lipid_molecule_species.values())):
                for leaflet, mask in (("upper", leaf), ("lower", ~leaf)):
                    availability[(species, leaflet)] = int(np.count_nonzero(
                        (self.species == species) & mask
                    ))
            for region, atoms in regions.items():
                hits = _contacts(atoms, self.lipid_contact, self.contact_A, box)
                observed = {}
                for j in hits[:, 1]:
                    lipid_atom = self.lipid_contact[int(j)]
                    molecule = self.lipid_molecule_by_atom.get(int(lipid_atom.index))
                    if molecule is None:
                        continue
                    anchor = self.lipid_anchor_index_by_molecule.get(molecule)
                    if anchor is None:
                        continue
                    species = self.lipid_molecule_species[molecule]
                    leaflet = "upper" if leaf[anchor] else "lower"
                    observed.setdefault((species, leaflet, "any"), set()).add(molecule)
                    if species in ("GM1a", "BGAL"):
                        # Glycan versus whole ceramide residue. 'Ceramide'
                        # includes its polar atoms: it is NOT a tail-only call.
                        moiety = ("ceramide" if lipid_atom.resname == "CER160"
                                  else "glycan")
                        observed.setdefault((species, leaflet, moiety), set()).add(molecule)
                for (species, leaflet), denominator in availability.items():
                    moieties = ("any", "glycan", "ceramide") if species in (
                        "GM1a", "BGAL"
                    ) else ("any",)
                    for moiety in moieties:
                        count = len(observed.get((species, leaflet, moiety), ()))
                        self.rows["lipid_specificity"].append({
                            "time_ps": t, "region": region, "species": species,
                            "leaflet": leaflet, "moiety": moiety,
                            "contacted_molecules": count,
                            "available_molecules": denominator,
                            "contact_fraction": (
                                count / denominator if denominator else np.nan
                            ),
                        })
        for name, ag in self.hs.items():
            attachment = _center(self.attach[name], box)
            distal = _center(self.hs_distal[name], box)
            vec = _vector(attachment, distal, box)
            sign = np.sign(_vector(
                [midpoint[0], midpoint[1], mid_z], attachment, box
            )[2]) or 1.
            lip_hits = _contacts(ag, self.lipid_contact, self.contact_A, box)
            site_hits = _contacts(self.attach[name], self.lipid_contact, self.contact_A, box)
            wat_hits = _contacts(ag, self.water, self.hydration_A, box)
            prefix = f"sdc4_hs_{name}_"
            row[prefix + "orientation_cos"] = float(np.dot(unit(vec), [0., 0., sign]))
            row[prefix + "extension_A"] = float(np.linalg.norm(vec))
            row[prefix + "lateral_extension_A"] = float(np.linalg.norm(vec[:2]))
            row[prefix + "lipid_contact_atom_pairs"] = len(lip_hits)
            row[prefix + "lipid_contact_lipids"] = (
                len(np.unique(self.lipid_contact.resindices[lip_hits[:, 1]]))
                if len(lip_hits) else 0
            )
            row[prefix + "attachment_lipid_contact_pairs"] = len(site_hits)
            row[prefix + "hydration_water_count"] = (
                len(np.unique(self.water.indices[wat_hits[:, 1]])) if len(wat_hits) else 0
            )
            row[prefix + "hydration_water_per_heavy_atom"] = (
                row[prefix + "hydration_water_count"] / len(ag)
            )
            row[prefix + "membrane_contact_state"] = int(len(lip_hits) > 0)
            for ion, atoms in self.ions.items():
                hits = _contacts(ag, atoms, self.ion_A, box)
                row[prefix + f"ion_{ion}_count"] = (
                    len(np.unique(atoms.indices[hits[:, 1]])) if len(hits) else 0
                )
        self.rows["timeseries"].append(row)
        self.frame_count += 1
        self.timings["compute"] += time.perf_counter() - start

    def flush(self, part):
        start = time.perf_counter()
        if self.rows["timeseries"]:
            atomic_parquet(
                self.cache / f"part-{part:06d}.parquet",
                pd.DataFrame(self.rows["timeseries"]),
            )
            self.rows["timeseries"].clear()
        if self.motif_lipid_residence_enabled:
            atomic_parquet(
                self.cache / f"specificity-{part:06d}.parquet",
                pd.DataFrame(self.rows["lipid_specificity"], columns=[
                    "time_ps", "region", "species", "leaflet", "moiety",
                    "contacted_molecules", "available_molecules",
                    "contact_fraction",
                ]),
            )
            self.rows["lipid_specificity"].clear()
        if self.motif_lipid_residence_enabled:
            atomic_parquet(
                self.cache / f"residence-{part:06d}.parquet",
                pd.DataFrame(self.rows["lipid_residence"], columns=[
                    "time_ps", "motif", "molecule_id", "species",
                ]),
            )
            self.rows["lipid_residence"].clear()
        _atomic_npz(
            self.cache / f"map-{part:06d}.npz",
            counts=self.map_counts, baselines=self.map_baselines,
            nframes=np.asarray(self.map_frames, dtype=int),
        )
        self.map_counts.fill(0)
        self.map_baselines.fill(0)
        self.map_frames = 0
        self.timings["chunk_write"] += time.perf_counter() - start

    def prune_from(self, part):
        for path in self.cache.glob("*"):
            match = re.match(r"(?:part|map|residence|specificity)-(\d+)\.(?:parquet|npz)$", path.name)
            if match and int(match.group(1)) >= part:
                path.unlink()

    def finalize(self):
        start = time.perf_counter()
        parts = sorted(self.cache.glob("part-*.parquet"))
        if not parts:
            raise ValueError("No SDC4 mechanism frames were processed")
        df = pd.concat((pd.read_parquet(p) for p in parts), ignore_index=True)
        df = df.sort_values("time_ps").reset_index(drop=True)
        # Optional direct motif SASA endpoint, using the *existing*
        # FreeSASA calculation instead of recomputing SASA or reading XTC.
        # Interpretation follows the atom selection used for global SASA.
        sasa_file = self.repout / "profiles/sasa/residue_sasa.parquet"
        if sasa_file.exists() and self.motifs:
            sasa = pd.read_parquet(sasa_file)
            if not sasa.empty:
                for name, info in self.motifs.items():
                    sub = sasa[sasa["topology_residue_index"].isin(
                        info["resindices"]
                    )]
                    if not sub.empty:
                        label = f"sdc4_motif_{name}_sasa_nm2"
                        values = sub.groupby("time_ps", as_index=False)[
                            "sasa_nm2"
                        ].sum().rename(columns={"sasa_nm2": label})
                        df = df.merge(values, on="time_ps", how="left",
                                      validate="one_to_one")
        out = ensure_dir(self.repout / "extensions" / "sdc4_mechanisms")
        atomic_parquet(out / "timeseries.parquet", df)
        counts = np.zeros_like(self.map_counts)
        baselines = np.zeros_like(self.map_baselines)
        nframes = 0
        for path in sorted(self.cache.glob("map-*.npz")):
            with np.load(path, allow_pickle=False) as data:
                counts += data["counts"]
                baselines += data["baselines"]
                nframes += int(data["nframes"])
        if nframes != len(df):
            raise ValueError("Lipid density map frame count disagrees with scalar timeseries")
        area = (self.edges[1] - self.edges[0])**2
        density = counts / nframes / area * 100.  # lipids per nm^2
        bulk = baselines / nframes * 100.
        fold = np.divide(density, bulk[:, None, None],
                         out=np.full_like(density, np.nan),
                         where=bulk[:, None, None] > 0)
        _atomic_npz(
            out / "lipid_density_maps.npz",
            counts=counts, density_per_nm2=density, enrichment_fold=fold,
            bulk_density_per_nm2=bulk,
            species=np.asarray([x[0] for x in self.map_names], str),
            leaflet=np.asarray([x[1] for x in self.map_names], str),
            x_edges_A=self.edges, y_edges_A=self.edges,
            nframes=np.asarray(nframes, int),
        )
        atomic_json(out / "selection_audit.json", self.audit)
        if self.motif_lipid_residence_enabled:
            from .motif_residence import summarize_lipid_residence
            parts = sorted(self.cache.glob("residence-*.parquet"))
            events = pd.concat((pd.read_parquet(p) for p in parts),
                               ignore_index=True) if parts else pd.DataFrame()
            atomic_parquet(out / "motif_lipid_contact_events.parquet", events)
            residence = summarize_lipid_residence(
                events, df["time_ps"].to_numpy(float)
            )
            atomic_parquet(out / "motif_lipid_residence.parquet", residence)
            pieces = sorted(self.cache.glob("specificity-*.parquet"))
            specificity = pd.concat(
                (pd.read_parquet(p) for p in pieces), ignore_index=True
            ) if pieces else pd.DataFrame()
            atomic_parquet(out / "lipid_region_specificity.parquet", specificity)

        # Contact-conditioned hydration is a *within-replicate* descriptive
        # comparison, not an independent-frame significance test.
        conditioned = []
        for name in self.hs:
            prefix = f"sdc4_hs_{name}_"
            for contact_state, sub in df.groupby(prefix + "membrane_contact_state"):
                item = {
                    "hs_chain": name, "membrane_contact_state": int(contact_state),
                    "frames": len(sub),
                    "hydration_water_per_heavy_atom_mean": float(
                        sub[prefix + "hydration_water_per_heavy_atom"].mean()
                    ),
                }
                for ion in self.ions:
                    item[f"ion_{ion}_count_mean"] = float(
                        sub[prefix + f"ion_{ion}_count"].mean()
                    )
                conditioned.append(item)
        atomic_parquet(out / "hs_conditioned_solvation.parquet",
                       pd.DataFrame(conditioned, columns=(
                           ["hs_chain", "membrane_contact_state", "frames",
                            "hydration_water_per_heavy_atom_mean"]
                           + [f"ion_{ion}_count_mean" for ion in self.ions]
                       )))
        self.timings["finalize"] += time.perf_counter() - start
        atomic_json(out / "timings.json", {
            "frames_this_invocation": self.frame_count,
            "compute_seconds": dict(self.timings),
        })
        return {"frames": nframes, "n_motifs": len(self.motifs),
                "n_hs_chains": len(self.hs), "n_map_species": len(set(self.species))}


def run_sdc4_mechanisms(rep, cfg, repout):
    """Standalone fallback when regular protein streaming is disabled."""
    if not enabled(cfg):
        return {}
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    tr = cfg["trajectory"]
    start = 0 if tr.get("start") is None else int(tr["start"])
    u.trajectory[start]
    observer = SDC4Mechanisms(u, cfg, repout)
    observer.prune_from(0)
    frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
    chunk = max(1, int(tr.get("chunk_frames", 1000)))
    count = 0
    for count, ts in enumerate(track(
        frames, total=len(frames), desc=f"{rep.system}/{rep.label}: SDC4 mechanisms",
        enabled=cfg.get("progress", {}).get("enabled", True),
    ), start=1):
        observer.update(ts, count - 1)
        if count % chunk == 0:
            observer.flush((count - 1) // chunk)
    if count % chunk:
        observer.flush(count // chunk)
    return observer.finalize()
