"""Shared periodic solvation neighbor search, compact RDFs and optional identity events.

Box-density-normalized radial distributions are NOT corrected for inaccessible volume.
All distances are Å and times ps; the actual water/ion atom selections are
explicit to avoid guessing CHARMM lipid and glycan atom chemistry.
"""
from collections import defaultdict
from itertools import combinations, product
from pathlib import Path
import os
import tempfile

import pyarrow.parquet as pq

import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.lib.mdamath import triclinic_vectors
from ..core.geometry import minimum_image_displacements

from ..config import atomic_json, atomic_parquet, ensure_dir

PRESETS = {"standard": (False, False), "detailed": (True, True),
           "kinetics": (True, True)}
COUNT_COLUMNS = ["frame", "time_ps", "solvent", "target_id", "target_class",
                 "site_class", "n_solvent"]
EVENT_COLUMNS = ["frame", "time_ps", "solvent", "target_id", "solvent_id"]
BRIDGE_COLUMNS = ["frame", "time_ps", "solvent", "solvent_id", "target_a", "target_b"]

# A fixed, symmetric in-plane grid makes the reference reproducible on restart.
SPATIAL_OFFSETS = tuple(pair for pair in product((0., 1/3, 2/3), repeat=2) if any(pair))
SPATIAL_COLUMNS = ["frame", "time_ps", "solvent", "target_a", "target_b",
                   "offset_u_fraction", "offset_v_fraction", "bridge_entity_count"]

EPISODE_COLUMNS = ["solvent", "target_id", "solvent_id", "start_frame", "end_frame",
                   "duration_ps", "left_censored", "right_censored",
                   "sampled_duration_ps", "duration_convention"]


def rdf_from_counts(counts, edges, density_exposure):
    """Observed pairs / expected pairs at the supplied reference density exposure."""
    shell_volume = (4 * np.pi / 3) * np.diff(np.asarray(edges, float) ** 3)
    if density_exposure <= 0:
        return np.full(len(counts), np.nan)
    return np.asarray(counts, float) / (shell_volume * density_exposure)


def aqueous_slab_density(solvent_positions, membrane_positions, box, half_width_A, axis=2):
    """Selected entities / geometric volume outside an explicit membrane slab.

    This is an aqueous-region reference, not a shell accessible-volume or
    protein/HS excluded-volume correction. Only axis-aligned orthorhombic
    bilayers are supported. An empty reference has unavailable density.
    """
    box = np.asarray(box, float)
    solvent = np.asarray(solvent_positions, float)
    membrane = np.asarray(membrane_positions, float)
    if (box.shape != (6,) or not np.isfinite(box).all() or np.any(box[:3] <= 0)
            or not np.allclose(box[3:], 90., rtol=0, atol=1e-5)):
        raise ValueError("aqueous_reference requires an orthorhombic periodic box")
    if axis not in (0, 1, 2) or not np.isfinite(half_width_A) or not 0 < half_width_A < box[axis]/2:
        raise ValueError("aqueous_reference half_width_A must leave a nonempty aqueous region")
    if (membrane.ndim != 2 or membrane.shape[1] != 3 or not len(membrane)
            or solvent.ndim != 2 or solvent.shape[1] != 3
            or not np.isfinite(membrane).all() or not np.isfinite(solvent).all()):
        raise ValueError("aqueous_reference needs finite membrane and solvent coordinates")
    length = box[axis]
    angles = 2*np.pi*membrane[:, axis]/length
    center_vector = np.mean(np.exp(1j*angles))
    if abs(center_vector) < 1e-5:
        raise ValueError("aqueous_reference membrane center is ambiguous; narrow membrane_selection")
    center = (np.angle(center_vector) % (2*np.pi))*length/(2*np.pi)
    distance = (solvent[:, axis]-center+length/2) % length-length/2
    count = int(np.count_nonzero(np.abs(distance) > half_width_A))
    volume = float(np.prod(box[:3])*(1-2*half_width_A/length))
    return (count/volume if count else np.nan), count, volume, center


def contact_episodes(events, observed_frames):
    """Continuous sampled-frame runs with both observed and imputed durations.

    Legacy duration_ps remains the observed first-to-last span. The added
    sampled_duration_ps follows the lipid-residence sampled interval rule.
    """
    from ..core.residence import sampled_episode_durations
    columns = EPISODE_COLUMNS
    if len(observed_frames) == 0:
        return pd.DataFrame(columns=columns)
    frames = [int(f) for f, _ in observed_frames]
    frame_times = np.asarray([float(t) for _, t in observed_frames])
    if len(set(frames)) != len(frames):
        raise ValueError("Duplicate observation frames")
    if not np.isfinite(frame_times).all() or np.any(np.diff(frame_times) <= 0):
        raise ValueError("Episode observation times must increase")
    order = {f: i for i, f in enumerate(frames)}
    pairs = defaultdict(set)
    for f, solvent, target, molecule in events:
        if int(f) not in order:
            raise ValueError("Contact event outside observation frames")
        pairs[(solvent, target, molecule)].add(int(f))
    rows, bounds = [], []
    for (solvent, target, molecule), active in pairs.items():
        previous = start = None
        for frame in sorted(active, key=order.__getitem__):
            if start is None:
                start = previous = frame
            elif order[frame] == order[previous] + 1:
                previous = frame
            else:
                bounds.append((order[start], order[previous]))
                rows.append((solvent, target, molecule, start, previous,
                             np.nan, start == frames[0], previous == frames[-1],
                             np.nan, "next_sample_interval_or_terminal_median"))
                start = previous = frame
        bounds.append((order[start], order[previous]))
        rows.append((solvent, target, molecule, start, previous,
                     np.nan, start == frames[0], previous == frames[-1],
                     np.nan, "next_sample_interval_or_terminal_median"))
    result = pd.DataFrame(rows, columns=columns)
    if bounds:
        first, last = np.asarray(bounds, dtype=int).T
        result["duration_ps"], result["sampled_duration_ps"] = sampled_episode_durations(
            frame_times, first, last
        )
    return result


def _targets(u, spec):
    result = []
    seen = set()
    for kind, groups in (spec.get("targets") or {}).items():
        for label, definition in groups.items():
            motif = definition.get("sequence") if isinstance(definition, dict) else None
            if motif is not None:
                from ..mechanisms.mechanistic_common import sequence_matches
                import re
                # Exact protein sequence, not regex: prevent accidental
                # partial/wildcard matching and verify chain-specific identity.
                pattern = str(motif).upper()
                if not re.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+", pattern):
                    raise ValueError(f"solvation motif {kind}/{label} must be an exact amino-acid sequence")
                selection = definition.get("selection", "protein and name CA")
                ca = u.select_atoms(selection)
                ca = ca.select_atoms("name CA")
                matches = sequence_matches(ca, {label: re.escape(pattern)})
                if not matches:
                    raise ValueError(f"solvation motif {kind}/{label} sequence {pattern} not found in topology")
                for match_id, residues in matches.items():
                    # Whole heavy-atom sidechains and backbone; not CA-only.
                    atoms = u.atoms[np.isin(u.atoms.resindices, residues)]
                    atoms = atoms.select_atoms(definition.get("atom_selection", "not name H*"))
                    if not len(atoms):
                        raise ValueError(f"solvation motif {kind}/{label} has no heavy atoms")
                    overlap = seen.intersection(map(int, atoms.indices))
                    if overlap:
                        raise ValueError(f"solvation motif {kind}/{label} overlaps another target")
                    seen.update(map(int, atoms.indices))
                    result.append((f"{kind}/{label}/{match_id}", kind, label, atoms))
                continue
            selection = definition if isinstance(definition, str) else definition["selection"]
            grouping = "residue" if isinstance(definition, str) else definition.get("group_by", "residue")
            atoms = u.select_atoms(selection)
            if not len(atoms):
                raise ValueError(f"solvation target {kind}/{label} has no atoms")
            if grouping == "residue":
                chunks = [(str(res.ix), atoms[atoms.resindices == res.ix])
                          for res in atoms.residues]
            elif grouping == "hs_chain":
                # Topology-derived HS chain identities (never chainID alone).
                # HS is attached to the protein, so generic bonded fragments
                # cannot distinguish HS chains connected to the same protomer.
                auto = spec.get("hs_auto") or {}
                if not auto.get("glycan_selection") or not auto.get("attachment_selection"):
                    raise ValueError("solvation hs_chain grouping requires hs_auto glycan and attachment selections")
                from ..mechanisms.hs_topology import discover_covalent_hs
                chains, _, _ = discover_covalent_hs(u, auto)
                chunks = []
                covered = set()
                selected = set(map(int, atoms.indices))
                for identity, glycan in chains.items():
                    atom_ids = sorted(selected.intersection(map(int, glycan.indices)))
                    if atom_ids:
                        chunks.append((identity, u.atoms[atom_ids]))
                        covered.update(atom_ids)
                if covered != selected:
                    raise ValueError("solvation hs_chain selected atoms not covered by discovered HS chains")
            elif grouping == "selection":
                chunks = [("all", atoms)]
            elif grouping == "fragment":
                try:
                    chunks = [(str(min(f.indices)), f) for f in atoms.fragments]
                except Exception as exc:
                    raise ValueError("solvation fragment grouping requires topology bonds") from exc
            else:
                raise ValueError("solvation group_by: use residue, fragment or selection")
            for ident, group in chunks:
                overlap = seen.intersection(map(int, group.indices))
                if overlap:
                    raise ValueError("solvation target selections overlap")
                seen.update(map(int, group.indices))
                result.append((f"{kind}/{label}/{ident}", kind, label, group))
    if not result:
        raise ValueError("solvation.targets must contain nonempty target groups")
    return result




def unique_target_solvent_pairs(target_indices, solvent_indices, n_targets):
    """Count unique solvent entities per target without Python per-pair sets.

    Return (counts, sorted unique local index pairs). No atom identities are
    discarded: the caller can map local solvent indices to topology indices
    before persistence. Duplicated coordination atoms are counted once.
    """
    target_indices = np.asarray(target_indices, dtype=np.intp)
    solvent_indices = np.asarray(solvent_indices, dtype=np.intp)
    if target_indices.shape != solvent_indices.shape:
        raise ValueError("target/solvent arrays must have matching shape")
    if len(target_indices) == 0:
        return np.zeros(n_targets, dtype=np.int64), np.empty((0, 2), dtype=np.intp)
    order = np.lexsort((solvent_indices, target_indices))
    targets = target_indices[order]
    solvent = solvent_indices[order]
    first = np.r_[True, (targets[1:] != targets[:-1]) |
                  (solvent[1:] != solvent[:-1])]
    unique_targets = targets[first]
    pairs = np.column_stack((unique_targets, solvent[first]))
    return np.bincount(unique_targets, minlength=n_targets), pairs


def water_donor_hbonds(site_xyz, water_xyz, water_hydrogens, pairs, box,
                       max_da_A=3.5, min_angle_deg=150.):
    """Return water-index/site-index/H-index triples satisfying D-H...A.

    Explicit water O-H bonds and target-acceptor atom selections are required.
    This deliberately does not classify target-donor/water-acceptor contacts.
    """
    pairs = np.asarray(pairs, dtype=int).reshape(-1, 2)
    if not len(pairs):
        return np.empty((0, 3), dtype=int)
    donor = np.asarray(water_xyz, float)[pairs[:, 1]]
    acceptor = np.asarray(site_xyz, float)[pairs[:, 0]]
    da = np.linalg.norm(minimum_image_displacements(acceptor-donor, box), axis=1)
    output = []
    for hi in (0, 1):
        hyd = np.asarray(water_hydrogens, float)[pairs[:, 1], hi]
        dh = minimum_image_displacements(donor-hyd, box)
        ah = minimum_image_displacements(acceptor-hyd, box)
        norms = np.linalg.norm(dh, axis=1)*np.linalg.norm(ah, axis=1)
        cosine = np.divide(np.einsum("ij,ij->i", dh, ah), norms,
                           out=np.zeros(len(da)), where=norms > 1e-12)
        angle = np.degrees(np.arccos(np.clip(cosine, -1., 1.)))
        output.extend((int(pairs[i,1]), int(pairs[i,0]), hi)
                      for i in np.flatnonzero((da <= max_da_A) &
                                               (angle >= min_angle_deg)))
    return np.asarray(output, dtype=int).reshape(-1, 3)


class SolvationObserver:
    """Fused observer. Histograms and compact tables are checkpointed per chunk."""

    def __init__(self, universe, cfg, repout):
        self.u = universe
        self.spec = cfg["analyses"]["solvation"]
        self.repout = Path(repout)
        preset = self.spec.get("preset", "standard")
        if preset not in PRESETS:
            raise ValueError("solvation.preset: expected standard, detailed or kinetics")
        events, bridges = PRESETS[preset]
        opts = self.spec.get("options") or {}
        self.events_on = bool(opts.get("events", events))
        self.bridges_on = bool(opts.get("bridges", bridges))
        self.water_hbond_bridges = bool(opts.get("water_hbond_bridges", False))
        if self.water_hbond_bridges:
            self.bridges_on = True
        if self.bridges_on:
            self.events_on = True
        hs_auto = self.spec.get("hs_auto") or (
            cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("hs_auto") or {})
        self.targets = _targets(universe, {**self.spec, "hs_auto": hs_auto})
        self.track_classes = opts.get("track_target_classes")
        if self.track_classes is not None:
            if not isinstance(self.track_classes, (tuple, list)) or not self.track_classes:
                raise ValueError("solvation.options.track_target_classes must be a nonempty list")
            known = {kind for _, kind, _, _ in self.targets}
            unknown = set(self.track_classes) - known
            if unknown:
                raise ValueError(f"unknown solvation tracking target classes: {sorted(unknown)}")
        self.tracked_target_indices = {
            i for i, (_, kind, _, _) in enumerate(self.targets)
            if self.track_classes is None or kind in self.track_classes}

        self.eligible_bridge_pairs = opts.get("bridge_target_pairs") or []
        tracked = {self.targets[i][0] for i in self.tracked_target_indices}
        normalized_pairs = []
        for pair in self.eligible_bridge_pairs:
            if (not isinstance(pair, (list, tuple)) or len(pair) != 2
                    or pair[0] == pair[1] or not set(pair) <= tracked):
                raise ValueError("bridge_target_pairs must name two distinct tracked target IDs")
            normalized_pairs.append(tuple(sorted(pair)))
        if len(set(normalized_pairs)) != len(normalized_pairs):
            raise ValueError("bridge_target_pairs contains duplicate pairs")
        if normalized_pairs and not self.bridges_on:
            raise ValueError("bridge_target_pairs requires bridges enabled")

        indices = np.concatenate([g.indices for _, _, _, g in self.targets])
        self.site_atoms = universe.atoms[indices]
        self.target_of_site = np.concatenate(
            [np.full(len(g), i, dtype=np.int32)
             for i, (_, _, _, g) in enumerate(self.targets)])
        self.solvents = {}
        selectors = {}
        if self.spec.get("water_selection"):
            selectors["water"] = self.spec["water_selection"]
        ion_specs = self.spec.get("ion_selections") or {}
        if ion_specs and self.spec.get("ion_selection"):
            raise ValueError("solvation: specify ion_selection OR ion_selections")
        if ion_specs:
            for label, selection in ion_specs.items():
                if not label or "/" in str(label):
                    raise ValueError("solvation ion species labels must be nonempty and exclude '/'")
                selectors["ions/" + str(label)] = selection
        elif self.spec.get("ion_selection"):
            selectors["ions"] = self.spec["ion_selection"]
        self.solvent_selections = selectors
        assigned = set()
        for name, query in selectors.items():
            group = universe.select_atoms(query)
            if not len(group) or len(np.unique(group.resindices)) != len(group):
                raise ValueError(f"solvation.{name}: select exactly one atom per solvent entity")
            overlap = set(map(int, group.indices)) & assigned
            if overlap:
                raise ValueError(f"solvation solvent species selections overlap: {name}")
            assigned.update(map(int, group.indices))
            if set(map(int, group.indices)).intersection(map(int, indices)):
                raise ValueError(f"solvation target sites overlap solvent selection {name}")
            self.solvents[name] = group
        if not self.solvents:
            raise ValueError("solvation requires water_selection and/or ion_selection")
        self.spatial_ion = opts.get("bridge_spatial_ion")
        self.spatial_rows = []
        if self.spatial_ion is not None:
            if (not isinstance(self.spatial_ion, str)
                    or not self.spatial_ion.startswith("ions/")
                    or self.spatial_ion not in self.solvents or not normalized_pairs):
                raise ValueError("bridge_spatial_ion requires a named ion species and bridge_target_pairs")
            axis = self.spec.get("aqueous_reference") or {}
            self.spatial_plane = [i for i, label in enumerate("xyz") if label != axis.get("axis", "z")]
            targets_by_id = {name: i for i, (name, _, _, _) in enumerate(self.targets)}
            selected = [targets_by_id[name] for pair in normalized_pairs for name in pair]
            mask = np.isin(self.target_of_site, selected)
            self.spatial_sites = self.site_atoms[mask]
            self.spatial_site_targets = self.target_of_site[mask]
            self.spatial_pairs = [(a, b, targets_by_id[a], targets_by_id[b])
                                  for a, b in normalized_pairs]
        self.aqueous_reference = self.spec.get("aqueous_reference")
        self.reference_membrane = None
        if self.aqueous_reference:
            if (not isinstance(self.aqueous_reference, dict)
                    or not self.aqueous_reference.get("membrane_selection")
                    or "half_width_A" not in self.aqueous_reference
                    or self.aqueous_reference.get("axis", "z") not in ("x", "y", "z")):
                raise ValueError("aqueous_reference requires membrane_selection, half_width_A and axis x/y/z")
            self.reference_membrane = universe.select_atoms(self.aqueous_reference["membrane_selection"])
            if not len(self.reference_membrane):
                raise ValueError("aqueous_reference membrane_selection matched zero atoms")
            self.reference_axis = {"x": 0, "y": 1, "z": 2}[self.aqueous_reference.get("axis", "z")]
            self.reference_half_width = float(self.aqueous_reference["half_width_A"])
            aqueous_slab_density(next(iter(self.solvents.values())).positions,
                self.reference_membrane.positions, universe.dimensions,
                self.reference_half_width, self.reference_axis)
        self.reference_rows = []
        self.aqueous_exposure = defaultdict(float)
        self.reference_invalid = defaultdict(int)
        self.hbond_rows = []
        self.water_hydrogens = None
        self.hbond_da_A = float(opts.get("hbond_da_cutoff_A", 3.5))
        self.hbond_angle = float(opts.get("hbond_angle_deg", 150.))
        if self.water_hbond_bridges:
            if "water" not in self.solvents:
                raise ValueError("water H-bond bridges require water_selection")
            if not (0 < self.hbond_da_A <= float(self.spec.get("search_radius_A", 8.))
                    and 0 < self.hbond_angle <= 180):
                raise ValueError("invalid water H-bond bridge thresholds")
            try:
                bonds = self.u.bonds.indices
            except Exception as exc:
                raise ValueError("water H-bond bridging requires O-H topology bonds") from exc
            waters = self.solvents["water"]
            water_owner = np.full(len(self.u.atoms), -1, dtype=int)
            water_owner[waters.indices] = np.arange(len(waters))
            directed = np.concatenate((bonds, bonds[:, ::-1]))
            owners = water_owner[directed[:, 0]]
            neighbors = directed[owners >= 0, 1]
            owners = owners[owners >= 0]
            hydrogen = np.fromiter(
                (str(name).upper().lstrip("0123456789").startswith("H")
                 for name in self.u.atoms[neighbors].names),
                dtype=bool, count=len(neighbors),
            )
            owners = owners[hydrogen]
            if not np.all(np.bincount(owners, minlength=len(waters)) == 2):
                raise ValueError("water oxygen must have precisely two bonded H atoms")
            self.water_hydrogens = self.u.atoms[neighbors[hydrogen][np.argsort(owners, kind="stable")]]


        radius = float(self.spec.get("search_radius_A", 8.0))
        width = float(self.spec.get("rdf_bin_width_A", 0.2))
        if not np.isfinite(radius) or not np.isfinite(width) or min(radius, width) <= 0:
            raise ValueError("solvation RDF radius and bin width must be positive")
        self.neighbor_method = self.spec.get("neighbor_method", "auto")
        if self.neighbor_method not in ("auto", "nsgrid", "pkdtree", "bruteforce"):
            raise ValueError("solvation.neighbor_method must be auto, nsgrid, pkdtree or bruteforce")
        self.edges = np.linspace(0, radius, int(np.ceil(radius / width)) + 1)
        self.shells = {"water": 3.5, "ions": 4.0, **(self.spec.get("shell_cutoffs_A") or {})}
        self.shells.update({key: self.shells.get(key, self.shells["ions"])
                            for key in self.solvents if key.startswith("ions/")})
        if any(not np.isfinite(float(x)) or not 0 < float(x) <= radius for x in self.shells.values()):
            raise ValueError("solvation shell cutoffs must be positive and <= search radius")
        # Optional local 3D WATER density maps. The 3 anchor atoms define a
        # co-moving orthonormal frame for a named target, not a global grid.
        self.grid_frames = {}
        self.grid_width = float(self.spec.get("grid_width_A", 12.0))
        self.grid_spacing = float(self.spec.get("grid_spacing_A", 0.75))
        if self.spec.get("grid_frames"):
            if "water" not in self.solvents:
                raise ValueError("solvation 3D hydration maps require water_selection")
            if not 0 < self.grid_spacing <= self.grid_width <= 24.0:
                raise ValueError("solvation grid width must be <=24 Å and spacing positive")
            n = int(np.ceil(self.grid_width / self.grid_spacing))
            if n > 64:
                raise ValueError("solvation grid exceeds 64 bins per axis")
            self.grid_edges = np.linspace(-self.grid_width/2, self.grid_width/2, n+1)
            known = {target for target, _, _, _ in self.targets}
            for target, selection in self.spec["grid_frames"].items():
                if target not in known:
                    raise ValueError(f"unknown 3D hydration target {target}")
                if isinstance(selection, list):
                    anchors = [universe.select_atoms(query) for query in selection]
                    if len(anchors) != 3 or any(len(group) != 1 for group in anchors):
                        raise ValueError(f"3D hydration frame {target} needs three single-atom selectors")
                    frame = universe.atoms[[int(group.indices[0]) for group in anchors]]
                else:
                    frame = universe.select_atoms(selection)
                if len(frame) != 3 or len(set(frame.indices)) != 3:
                    raise ValueError(f"3D hydration frame {target} needs exactly 3 ordered atoms")
                target_atoms = next(g for name, _, _, g in self.targets if name == target)
                if int(frame.indices[0]) not in set(map(int, target_atoms.indices)):
                    raise ValueError("3D grid origin must belong to the target site selection")
                self.grid_frames[target] = frame
        self.grid_counts = {key: np.zeros((len(self.grid_edges)-1,)*3, dtype=np.int64)
                            for key in self.grid_frames} if self.grid_frames else {}
        self.grid_exposure = 0.0
        self.hist = defaultdict(lambda: np.zeros(len(self.edges) - 1, np.int64))
        self.exposure = defaultdict(float)
        self.rows, self.events, self.bridges, self.frame_rows = [], [], [], []
        self.frames = 0
        self.class_names = list(dict.fromkeys(
            (kind, label) for _, kind, label, _ in self.targets))
        class_idx = {name: i for i, name in enumerate(self.class_names)}
        self.class_of_target = np.array(
            [class_idx[(kind, label)] for _, kind, label, _ in self.targets],
            dtype=np.int32)
        self.class_atoms = np.bincount(
            self.class_of_target,
            weights=[len(group) for _, _, _, group in self.targets],
            minlength=len(self.class_names))

    def update(self, ts, frame_index, time_ps):
        box = np.asarray(self.u.dimensions, float)
        if box.shape != (6,) or not np.all(np.isfinite(box)) or np.any(box[:3] <= 0):
            raise ValueError("solvation needs a valid periodic box")
        volume = abs(float(np.linalg.det(triclinic_vectors(box))))
        if volume <= 0:
            raise ValueError("solvation box volume must be positive")
        self.frame_rows.append((int(frame_index), float(time_ps)))
        for solvent, atoms in self.solvents.items():
            density_reference = np.nan
            if self.reference_membrane is not None:
                if self.edges[-1] > np.min(box[:3])/2:
                    raise ValueError("aqueous_reference RDF radius exceeds the inscribed periodic sphere")
                density_reference, count_reference, volume_reference, center_reference = aqueous_slab_density(
                    atoms.positions, self.reference_membrane.positions, box,
                    self.reference_half_width, self.reference_axis)
                self.reference_rows.append((int(frame_index), float(time_ps), solvent,
                    count_reference, volume_reference, density_reference, center_reference,
                    "valid" if np.isfinite(density_reference) else "empty_reference"))
            pair, distances = capped_distance(
                self.site_atoms.positions, atoms.positions,
                max_cutoff=float(self.edges[-1]), box=box, return_distances=True,
                method=None if self.neighbor_method == "auto" else self.neighbor_method)
            indices = self.target_of_site[pair[:, 0]] if len(pair) else np.empty(0, int)
            shell = distances <= float(self.shells[solvent])
            counts, unique_pairs = unique_target_solvent_pairs(
                indices[shell], pair[shell, 1], len(self.targets))
            # Only build per-solvent identity sets when events or bridges
            # actually require them. The standard preset stays NumPy-only.
            contacts = None
            if self.events_on or self.bridges_on:
                contacts = defaultdict(set)
                for target_idx, solvent_local in unique_pairs:
                    if int(target_idx) in self.tracked_target_indices:
                        contacts[int(target_idx)].add(
                            int(atoms.indices[solvent_local]))
            # Bin all site–solvent pairs in a single vectorized pass.
            # Per-target boolean masking here would scale as O(targets*pairs)
            # and is particularly expensive for hundreds of lipid molecules.
            radial_bin = np.searchsorted(self.edges, distances, side="right")-1
            valid = (radial_bin >= 0) & (radial_bin < len(self.edges)-1)
            n_bins = len(self.edges)-1
            flat = self.class_of_target[indices[valid]] * n_bins + radial_bin[valid]
            histogram = np.bincount(flat, minlength=len(self.class_names)*n_bins).reshape(
                len(self.class_names), n_bins)
            for class_index, (kind, label) in enumerate(self.class_names):
                key = (solvent, kind, label)
                self.hist[key] += histogram[class_index]
                self.exposure[key] += self.class_atoms[class_index] * len(atoms)/volume
                if self.reference_membrane is not None:
                    if np.isfinite(density_reference):
                        self.aqueous_exposure[key] += self.class_atoms[class_index]*density_reference
                    else:
                        self.reference_invalid[key] += 1
            for i, (target, kind, label, group) in enumerate(self.targets):
                self.rows.append((int(frame_index), float(time_ps), solvent,
                                  target, kind, label, int(counts[i])))
                if self.events_on and i in self.tracked_target_indices:
                    self.events.extend((int(frame_index), float(time_ps), solvent,
                                        target, solvent_id)
                                       for solvent_id in sorted(contacts[i]))
            if solvent == "water" and self.grid_frames:
                # Reuse candidate pairs from the same periodic neighbor search.
                # For each target, map water positions to its local anchor frame.
                water_xyz = atoms.positions
                for target_i, (target_id, _, _, _) in enumerate(self.targets):
                    anchors = self.grid_frames.get(target_id)
                    if anchors is None:
                        continue
                    origin, second, third = anchors.positions
                    v1, v2 = minimum_image_displacements(
                        np.asarray([second-origin, third-origin], dtype=float), box)
                    norm1 = np.linalg.norm(v1)
                    if norm1 < 1e-7:
                        raise ValueError(f"degenerate 3D grid frame for {target_id}")
                    e1 = v1/norm1
                    v2 = v2-np.dot(v2,e1)*e1
                    norm2 = np.linalg.norm(v2)
                    if norm2 < 1e-7:
                        raise ValueError(f"collinear 3D grid frame for {target_id}")
                    e2 = v2/norm2
                    rotation = np.column_stack((e1, e2, np.cross(e1,e2)))
                    js = np.unique(pair[indices == target_i, 1])
                    if len(js):
                        delta = minimum_image_displacements(
                            np.asarray(water_xyz[js]-origin,dtype=float),box)
                        local = delta @ rotation
                        # Inscribed sphere is fully covered by the original
                        # site-to-water search because origin belongs to target.
                        # Cube corners outside that sphere have no sampling
                        # guarantee and must never be interpreted as zero density.
                        local = local[np.linalg.norm(local, axis=1) <= self.grid_width/2]
                        self.grid_counts[target_id] += np.histogramdd(
                            local, bins=(self.grid_edges,)*3)[0].astype(np.int64)
                self.grid_exposure += len(atoms)/volume
            if solvent == "water" and self.water_hbond_bridges:
                candidate = pair[distances <= self.hbond_da_A]
                qualified = water_donor_hbonds(
                    self.site_atoms.positions, atoms.positions,
                    self.water_hydrogens.positions.reshape(-1, 2, 3),
                    candidate, box, max_da_A=self.hbond_da_A,
                    min_angle_deg=self.hbond_angle)
                by_water = defaultdict(lambda: defaultdict(set))
                for water_idx, site_idx, hydrogen_idx in qualified:
                    target_idx = int(self.target_of_site[site_idx])
                    if target_idx in self.tracked_target_indices:
                        by_water[int(water_idx)][target_idx].add(int(hydrogen_idx))
                for water_idx, target_hydrogens in by_water.items():
                    for a, b in combinations(sorted(target_hydrogens), 2):
                        # Two distinct donor H atoms must participate.
                        if any(h1 != h2 for h1 in target_hydrogens[a]
                               for h2 in target_hydrogens[b]):
                            self.hbond_rows.append(
                                (int(frame_index), float(time_ps),
                                 int(atoms.indices[water_idx]),
                                 self.targets[a][0], self.targets[b][0]))
            if solvent == self.spatial_ion:
                if not np.allclose(box[3:], 90., rtol=0, atol=1e-5):
                    raise ValueError("bridge_spatial_ion requires an orthorhombic periodic box")
                positions = atoms.positions
                sites = self.spatial_sites.positions
                cutoff = float(self.shells[solvent])
                method = None if self.neighbor_method == "auto" else self.neighbor_method
                for du, dv in SPATIAL_OFFSETS:
                    shift = np.zeros(3)
                    shift[self.spatial_plane] = np.array([du, dv]) * box[self.spatial_plane]
                    neighbors = capped_distance(sites, (positions + shift) % box[:3],
                        max_cutoff=cutoff, box=box, return_distances=False, method=method)
                    target_ids = self.spatial_site_targets[neighbors[:, 0]]
                    ion_ids = neighbors[:, 1]
                    for left, right, a, b in self.spatial_pairs:
                        count = len(np.intersect1d(ion_ids[target_ids == a], ion_ids[target_ids == b]))
                        self.spatial_rows.append((int(frame_index), float(time_ps), solvent,
                                                  left, right, du, dv, count))
            if self.bridges_on:
                inverse = defaultdict(set)
                for target_index, solvent_ids in contacts.items():
                    if target_index not in self.tracked_target_indices:
                        continue
                    for sid in solvent_ids:
                        inverse[sid].add(target_index)
                for sid, target_indices in inverse.items():
                    for left, right in combinations(sorted(target_indices), 2):
                        self.bridges.append((int(frame_index), float(time_ps), solvent, sid,
                                             self.targets[left][0], self.targets[right][0]))
        self.frames += 1

    def flush(self, part):
        root = ensure_dir(self.repout / "cache" / "solvation")
        for name, values, columns in (
            ("counts", self.rows, COUNT_COLUMNS),
            ("events", self.events, EVENT_COLUMNS),
            ("bridges", self.bridges, BRIDGE_COLUMNS),
            ("water_hbond_bridges", self.hbond_rows,
             ["frame", "time_ps", "water_id", "target_a", "target_b"]),
            ("bridge_spatial_reference", self.spatial_rows, SPATIAL_COLUMNS),
            ("frames", self.frame_rows, ["frame", "time_ps"]),
        ):
            if values or name in {"counts", "frames"}:
                atomic_parquet(root / f"{name}-{part:06d}.parquet",
                               pd.DataFrame(values, columns=columns))
            values.clear()
        history = [(solvent, kind, label, *map(int, bins), float(self.exposure[(solvent, kind, label)]))
                   for (solvent, kind, label), bins in sorted(self.hist.items())]
        atomic_parquet(root / f"hist-{part:06d}.parquet",
            pd.DataFrame(history, columns=["solvent", "target_class", "site_class"]
                         + [f"bin_{i}" for i in range(len(self.edges)-1)] + ["exposure"]))
        if self.reference_membrane is not None:
            atomic_parquet(root / f"aqueous_reference-{part:06d}.parquet", pd.DataFrame(
                self.reference_rows, columns=["frame", "time_ps", "solvent", "n_reference",
                    "reference_volume_A3", "density_per_A3", "membrane_center_A", "status"]))
            self.reference_rows.clear()
            atomic_parquet(root / f"aqueous_exposure-{part:06d}.parquet", pd.DataFrame([
                (*key, self.aqueous_exposure[key], self.reference_invalid[key]) for key in sorted(self.hist)
            ], columns=["solvent", "target_class", "site_class", "exposure", "invalid_frames"]))
            self.aqueous_exposure.clear()
            self.reference_invalid.clear()
        if self.grid_frames:
            records = []
            for target, grid in self.grid_counts.items():
                i,j,k = np.nonzero(grid)
                records.extend((target,int(x),int(y),int(z),int(grid[x,y,z]))
                               for x,y,z in zip(i,j,k))
                grid.fill(0)
            atomic_parquet(root / f"grid-{part:06d}.parquet", pd.DataFrame(
                records, columns=["target_id","ix","iy","iz","count"]))
            atomic_parquet(root / f"grid_exposure-{part:06d}.parquet",
                           pd.DataFrame({"exposure":[self.grid_exposure]}))
            self.grid_exposure = 0.
        self.hist.clear()
        self.exposure.clear()

    def prune_from(self, part):
        for path in (self.repout / "cache" / "solvation").glob("*.parquet"):
            try:
                number = int(path.stem.rsplit("-", 1)[1])
            except (IndexError, ValueError):
                continue
            if number >= part:
                path.unlink()

    def finalize(self):
        cache = self.repout / "cache" / "solvation"
        dest = ensure_dir(self.repout / "extensions" / "solvation")
        def read(name):
            files = sorted(cache.glob(f"{name}-*.parquet"))
            return pd.concat((pd.read_parquet(p) for p in files), ignore_index=True) if files else pd.DataFrame()
        for name in ("counts", "events", "bridges", "water_hbond_bridges"):
            files = sorted(cache.glob(f"{name}-*.parquet"))
            if files:
                # Keep the single-Parquet public output but stream chunk
                # row groups without concatenating every frame in RAM.
                target = dest / f"{name}.parquet"
                with tempfile.NamedTemporaryFile(dir=dest, suffix=".parquet",
                                                 delete=False) as handle:
                    temporary = Path(handle.name)
                try:
                    with pq.ParquetWriter(str(temporary),
                                          pq.read_schema(files[0])) as writer:
                        for chunk in files:
                            writer.write_table(pq.read_table(chunk))
                    os.replace(temporary, target)
                finally:
                    temporary.unlink(missing_ok=True)
            elif name == "counts":
                atomic_parquet(dest / "counts.parquet",
                               pd.DataFrame(columns=COUNT_COLUMNS))
        hist = read("hist")
        rows = []
        if not hist.empty:
            for (solvent, kind, label), data in hist.groupby(["solvent", "target_class", "site_class"]):
                counts = data[[f"bin_{i}" for i in range(len(self.edges)-1)]].sum(axis=0).to_numpy(dtype=float)
                g = rdf_from_counts(counts, self.edges, float(data.exposure.sum()))
                rows.extend(zip([solvent]*len(g), [kind]*len(g), [label]*len(g),
                                self.edges[:-1], self.edges[1:], counts, g))
        atomic_parquet(dest / "rdf.parquet", pd.DataFrame(rows, columns=[
            "solvent", "target_class", "site_class", "r_inner_A", "r_outer_A",
            "pair_count", "g_bulk"]))
        if self.reference_membrane is not None:
            reference = read("aqueous_exposure")
            rdf_path = dest / "rdf.parquet"
            rdf = pd.read_parquet(rdf_path)
            rdf["g_aqueous_reference"] = np.nan
            for key, group in reference.groupby(["solvent", "target_class", "site_class"]):
                mask = np.logical_and.reduce([rdf[column].eq(value) for column, value in
                    zip(("solvent", "target_class", "site_class"), key)])
                if not group.invalid_frames.sum():
                    rdf.loc[mask, "g_aqueous_reference"] = rdf_from_counts(
                        rdf.loc[mask, "pair_count"], self.edges, float(group.exposure.sum()))
            atomic_parquet(rdf_path, rdf)
            atomic_parquet(dest / "aqueous_reference.parquet", read("aqueous_reference"))
        if self.grid_frames:
            grids = read("grid")
            exposure = read("grid_exposure").exposure.sum()
            if not grids.empty:
                grids = grids.groupby(["target_id","ix","iy","iz"],as_index=False).agg(
                    count=("count","sum"))
                width = float(self.grid_edges[1]-self.grid_edges[0])
                for dim in "xyz":
                    grids[dim+"_A"] = self.grid_edges[grids["i"+dim].to_numpy(int)] + width/2
                grids["g_bulk"] = grids["count"]/(width**3*exposure)
                grids["valid_spherical_domain"] = True
            atomic_parquet(dest / "water_density_3d.parquet", grids)
        # Persist target-to-topology-residue membership for SASA reuse.
        membership = [(target, int(res.ix)) for target, _, _, group in self.targets
                      for res in group.residues]
        atomic_parquet(dest / "target_residues.parquet", pd.DataFrame(
            membership, columns=["target_id", "topology_residue_index"]))
        atomic_parquet(dest / "target_atoms.parquet", pd.DataFrame([
            (target, kind, label, int(atom.index), int(atom.resindex),
             str(atom.resname), str(atom.name), str(getattr(atom, "element", "")),
             float(getattr(atom, "charge", np.nan)))
            for target, kind, label, atoms in self.targets for atom in atoms
        ], columns=["target_id", "target_class", "site_class", "atom_index", "topology_residue_index",
                    "resname", "atom_name", "element", "charge_e"]))
        if self.events_on:
            frames = read("frames").drop_duplicates("frame").sort_values("frame")
            from .solvation_postprocess import write_contact_episodes
            write_contact_episodes(sorted(cache.glob("events-*.parquet")),
                list(frames[["frame", "time_ps"]].itertuples(index=False, name=None)),
                dest / "episodes.parquet")
        from .solvation_postprocess import write_solvation_derivatives
        write_solvation_derivatives(dest, self.repout,
            sasa_run=self.spec.get("sasa_run", "default"))
        if self.eligible_bridge_pairs:
            from .solvation_postprocess import bridge_pair_occupancy, bridge_pair_evidence
            frame_times = read("frames").drop_duplicates("frame").sort_values("frame").time_ps.to_numpy(float)
            pair_summary = bridge_pair_occupancy(
                read("bridges"), frame_times, self.eligible_bridge_pairs, self.solvents)
            observations = read("frames").drop_duplicates("frame").sort_values("frame")
            evidence = bridge_pair_evidence(sorted(cache.glob("events-*.parquet")),
                list(observations[["frame", "time_ps"]].itertuples(index=False, name=None)),
                self.eligible_bridge_pairs, self.solvents)
            keys = ["solvent", "target_a", "target_b"]
            reconstructed = evidence[0].groupby(keys, as_index=False).agg(
                frames=("bridge_present_fraction", "sum"), count=("bridge_entity_count", "sum"))
            matched = pair_summary.merge(reconstructed, on=keys, validate="one_to_one")
            if (len(matched) != len(pair_summary)
                    or not np.array_equal(matched.bridge_frames, matched.frames)
                    or not np.array_equal(np.rint(matched.mean_simultaneous_bridging_ions*
                        matched.observed_frames).astype(np.int64), matched['count'])):
                raise ValueError("Bridge records disagree with persisted target-contact identities")
            if self.spatial_ion is not None:
                spatial = read("bridge_spatial_reference")
                trace = evidence[0]
                frame_keys = keys + ["frame", "time_ps"]
                expected = trace.loc[trace.solvent == self.spatial_ion].merge(
                    observations, on="time_ps", validate="many_to_one")[frame_keys]
                if spatial.empty or spatial.duplicated(frame_keys + SPATIAL_COLUMNS[5:7]).any():
                    raise ValueError("Missing or duplicate spatial bridge reference records")
                fractions = (spatial.assign(present=spatial.bridge_entity_count.gt(0))
                             .groupby(frame_keys, as_index=False)
                             .agg(shifted_bridge_fraction=("present", "mean"),
                                  offsets=("present", "size")))
                coverage = expected.merge(fractions, on=frame_keys, how="outer", validate="one_to_one",
                                          indicator=True)
                if not coverage._merge.eq("both").all() or not coverage.offsets.eq(len(SPATIAL_OFFSETS)).all():
                    raise ValueError("Incomplete spatial bridge reference coverage")
                trace = trace.merge(fractions.drop(columns=["offsets", "frame"]),
                                    on=keys + ["time_ps"], how="left",
                                    validate="one_to_one")
                trace["observed_minus_spatial_bridge_fraction"] = (
                    trace.bridge_present_fraction - trace.shifted_bridge_fraction)
                atomic_parquet(dest / "bridge_spatial_reference.parquet", spatial)
                evidence = (trace, *evidence[1:])
            atomic_parquet(dest / "bridge_pair_summary.parquet", pair_summary)
            for name, data in zip(("bridge_pair_timeseries", "bridge_pair_shift_reference",
                                   "bridge_pair_control_summary"), evidence):
                atomic_parquet(dest / f"{name}.parquet", data)
        atomic_json(dest / "metadata.json", {
            "preset": self.spec.get("preset", "standard"), "edges_A": self.edges.tolist(),
            "neighbor_method": self.neighbor_method,
            "shell_cutoffs_A": self.shells, "events": self.events_on,
            "3d_grid_frames": list(self.grid_frames),
            "3d_grid_domain": "Only observed voxels within an inscribed sphere around a target origin; other voxels are not sampled.",
            "3d_grid_edges_A": self.grid_edges.tolist() if self.grid_frames else None,
            "bridges": self.bridges_on,
            "bridge_target_pairs": self.eligible_bridge_pairs,
            **({"bridge_spatial_ion": self.spatial_ion,
                "bridge_spatial_axes": self.spatial_plane,
                "bridge_spatial_offsets": SPATIAL_OFFSETS,
                "bridge_spatial_limitation": "Rigid periodic in-plane ion shifts preserve height and ion–ion distances but may put ions in excluded volume. Descriptive reference, not a chemical cooperativity or affinity test."}
               if self.spatial_ion is not None else {}),
            "bridge_control": "Deterministic nonzero circular offsets of target-B contact traces, at most 100; regular windows with >=4 observations only. Reference span is not a CI or a spatial/chemical cooperativity test.",
            "track_target_classes": self.track_classes,
            "solvent_selection_audit": {
                name: {"selection": self.solvent_selections[name], "n_entities": len(atoms),
                       "atom_names": sorted(set(map(str, atoms.names))),
                       "resnames": sorted(set(map(str, atoms.resnames)))}
                for name, atoms in self.solvents.items()},
            "water_hbond_bridges": self.water_hbond_bridges,
            "density_reference": "selected solvent entities / instantaneous full periodic cell volume",
            "aqueous_reference": self.aqueous_reference,
            "aqueous_reference_limitation": "Geometric volume outside a configured axis-aligned membrane slab; protein/HS excluded volume and shell accessibility are not corrected. Empty reference in any frame makes the entire aqueous-reference curve unavailable. Not evidence of a bulk plateau or affinity.",
            "rdf_value_column": "g_bulk (legacy name; full-cell average density reference)",
            "rdf_limitation": "Not independently measured aqueous bulk density; no excluded-volume correction. Evaluate alongside SASA and membrane geometry.",
            "residence_limitation": "Episodes represent consecutive sampled contacts; temporal resolution is trajectory frame interval.",
        })
        return {"frames": len(read("frames")), "fused": True}

