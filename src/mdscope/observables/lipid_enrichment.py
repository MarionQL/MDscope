"""Lipid composition and proximity relative to the TM footprint in each leaflet."""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from MDAnalysis.exceptions import NoDataError
from scipy.spatial import cKDTree

from ..config import atomic_json, atomic_parquet
from ..core.atom_selection import select_heavy_atoms
from ..core.lipid_geometry import circular_axis_mean

FOOTPRINT_OUTPUT = "extensions/lipid_footprint"


def footprint_bands(anchors, protein, lengths, axis, edges, radius, half_width, grid):
    """Assign anchors and quadrature samples using the same periodic disk union."""
    lateral = np.arange(3) != axis
    cell = lengths[lateral]
    height = circular_axis_mean(anchors[:, axis], lengths[axis])
    offset = (protein[:, axis] - height + lengths[axis] / 2) % lengths[axis] - lengths[axis] / 2
    selected = protein[np.abs(offset) <= half_width][:, lateral]
    if not len(selected):
        raise ValueError("No footprint atoms in the leaflet slab. Check footprint_selection and slab_half_width_A.")
    origin = selected[0]
    tree = cKDTree((selected - origin) % cell, boxsize=cell)
    points = np.concatenate(((anchors[:, lateral] - origin) % cell, grid * cell))
    distance = tree.query(points, k=1)[0] - radius
    band = np.searchsorted(edges, distance, side="right") - 1
    sample_bands = band[len(anchors):]
    accessible = np.bincount(sample_bands[sample_bands >= 0], minlength=len(edges))
    pixel_area = cell.prod() / len(grid)
    return band[:len(anchors)], accessible[:-1] * pixel_area, accessible.sum() * pixel_area


def shell_fold(observed, local_total, bulk_count, bulk_total):
    """Local mole fraction divided by the same leaflet's mole fraction."""
    return (observed * bulk_total / (local_total * bulk_count)
            if local_total and bulk_count else np.nan)


def anchored_lipid_owners(universe, residues):
    """Return per-atom anchor IDs using native molecule IDs when available."""
    try:
        molecule_ids = universe.atoms.molnums
        anchors = universe.residues[residues].molnums
    except NoDataError:
        selected_residues = universe.residues[residues]
        residue_atoms = selected_residues.atoms
        owner = np.full(len(universe.atoms), -1, dtype=int)
        owner[residue_atoms.indices] = np.repeat(
            np.arange(len(residues)), [len(res.atoms) for res in selected_residues],
        )
        try:
            fragments = residue_atoms.fragments
        except NoDataError:
            return owner
        for fragment in fragments:
            members = owner[fragment.indices]
            selected = members[members >= 0]
            if len(selected) and selected.min() != selected.max():
                raise ValueError("Lipid footprint requires one anchor per bonded molecule")
            if len(selected):
                owner[fragment.indices] = selected[0]
        return owner
    if len(np.unique(anchors)) != len(residues):
        raise ValueError("Lipid footprint requires one anchor per topology molecule")
    lookup = np.full(molecule_ids.max() + 1, -1, dtype=int)
    lookup[anchors] = np.arange(len(residues))
    return lookup[molecule_ids]


def footprint_reference(spec):
    return {
        "footprint_selection": spec["footprint_selection"], "axis": str(spec.get("axis", "z")).lower(),
        "species_selections": spec["species_selections"],
        "allow_absent_species": bool(spec.get("allow_absent_species", False)),
        "shell_edges_A": list(map(float, spec.get("shell_edges_A", [0, 5, 10, 20]))),
        "exclusion_radius_A": float(spec.get("exclusion_radius_A", 3.0)),
        "slab_half_width_A": float(spec.get("slab_half_width_A", 8.0)),
        "area_grid_size": int(spec.get("area_grid_size", 64)),
        "contact_cutoff_A": float(spec.get("contact_cutoff_A", 5.0)),
    }


def footprint_pending(path, spec):
    reference = path.parent / "reference.json"
    if not path.is_file() or not reference.is_file():
        return True
    return json.loads(reference.read_text()) != footprint_reference(spec)


class LipidFootprint:
    def __init__(self, universe, spec, species, residues, axis):
        self.reference = footprint_reference(spec)
        self.axis = axis
        self.protein = select_heavy_atoms(
            universe.select_atoms(spec["footprint_selection"]), policy="mechanistic",
        )
        if not len(self.protein):
            raise ValueError("lipid_enrichment.footprint_selection has no heavy atoms")
        self.edges = np.asarray(self.reference["shell_edges_A"])
        self.radius = self.reference["exclusion_radius_A"]
        self.half_width = self.reference["slab_half_width_A"]
        self.cutoff = self.reference["contact_cutoff_A"]
        size = self.reference["area_grid_size"]
        fraction = (np.arange(size) + .5) / size
        self.grid = np.stack(np.meshgrid(fraction, fraction), axis=-1).reshape(-1, 2)
        self.categories = sorted(spec["species_selections"])
        self.species_codes = np.searchsorted(self.categories, species)
        self.n_lipids = len(residues)

        owner = anchored_lipid_owners(universe, residues)
        if np.any(owner[self.protein.indices] >= 0):
            raise ValueError("Lipid footprint and anchored lipid atoms must be disjoint")
        self.lipids = select_heavy_atoms(universe.atoms[owner >= 0], policy="mechanistic")
        self.owners = owner[self.lipids.indices]
        if len(np.unique(self.owners)) != self.n_lipids:
            raise ValueError("Each lipid requires a heavy atom for proximity measurement")

    def frame(self, time, anchors, upper, changes, box):
        if not np.allclose(box[3:], 90., rtol=0, atol=1e-5):
            raise ValueError("Lipid footprint requires an orthorhombic box")
        lengths = box[:3]
        protein = self.protein.positions
        tree = cKDTree(protein % lengths, boxsize=lengths)
        distance = np.full(self.n_lipids, np.inf)
        np.minimum.at(distance, self.owners, tree.query(self.lipids.positions % lengths)[0])
        in_contact = distance <= self.cutoff
        series = {"time_ps": time, "lipid_footprint_leaflet_changes": changes}
        rows = []
        n_species = len(self.categories)
        n_bands = len(self.edges) - 1
        for leaf, mask in (("upper", upper), ("lower", ~upper)):
            band, areas, total_area = footprint_bands(
                anchors[mask], protein, lengths, self.axis, self.edges,
                self.radius, self.half_width, self.grid,
            )
            codes = self.species_codes[mask]
            bulk = np.bincount(codes, minlength=n_species)
            inside = (band >= 0) & (band < n_bands)
            local = np.bincount(band[inside] * n_species + codes[inside],
                                minlength=n_bands * n_species).reshape(n_bands, n_species)
            totals = local.sum(axis=1)
            contacted = np.bincount(self.species_codes[mask & in_contact], minlength=n_species)
            leaf_distances = distance[mask]
            for j, name in enumerate(self.categories):
                prefix = f"lipid_footprint_{name}_{leaf}"
                series[prefix + "_contact_fraction"] = contacted[j] / bulk[j] if bulk[j] else np.nan
                species_distances = leaf_distances[codes == j]
                series[prefix + "_mean_min_distance_A"] = species_distances.mean() if bulk[j] else np.nan
                series[prefix + "_closest_distance_A"] = species_distances.min() if bulk[j] else np.nan
                for i, (low, high) in enumerate(zip(self.edges[:-1], self.edges[1:])):
                    count = local[i, j]
                    fold = shell_fold(count, totals[i], bulk[j], len(codes))
                    density = (count * total_area / (bulk[j] * areas[i])
                               if bulk[j] and areas[i] > 0 and total_area > 0 else np.nan)
                    label = f"{prefix}_{low:g}_{high:g}A"
                    series[label + "_fold"] = fold
                    series[label + "_area_density_ratio"] = density
                    rows.append(dict(time_ps=time, leaflet=leaf, species=name,
                                     lower_A=low, upper_A=high, count_local=count,
                                     n_local=totals[i], count_bulk=bulk[j], n_bulk=len(codes),
                                     count_contact=contacted[j],
                                     mean_min_distance_A=series[prefix + "_mean_min_distance_A"],
                                     bulk_fraction=bulk[j] / len(codes), enrichment_fold=fold,
                                     area_density_ratio=density, shell_accessible_area_A2=areas[i],
                                     leaflet_accessible_area_A2=total_area))
        return series, rows

    def write_reference(self, out):
        atomic_json(out / "reference.json", self.reference)


def run_lipid_enrichment(rep, cfg, repout):
    from ..pipeline.fused_membrane import run_fused_membrane

    return run_fused_membrane(rep, cfg, repout, names=("lipid_enrichment",))["lipid_enrichment"]


def summarize_lipid_composition(repout):
    """Offline population inventory from counted molecular anchors."""
    directory = Path(repout) / FOOTPRINT_OUTPUT
    source = directory / "counts.parquet"
    if not source.exists():
        return
    counts = pd.read_parquet(source)
    if counts.empty:
        return
    keys = ["leaflet", "species", "lower_A", "upper_A"]
    if counts.duplicated(["time_ps", *keys]).any():
        raise ValueError(f"Duplicate lipid composition observations: {source}")
    summary = counts.assign(available=counts.count_bulk > 0).groupby(keys, as_index=False).agg(
        sampled_frames=("time_ps", "size"),
        mean_molecules=("count_bulk", "mean"),
        minimum_molecules=("count_bulk", "min"), maximum_molecules=("count_bulk", "max"),
        mean_leaflet_fraction=("bulk_fraction", "mean"),
        fraction_frames_available=("available", "mean"),
        mean_local_molecules=("count_local", "mean"),
    )
    atomic_parquet(directory / "composition_summary.parquet", summary)
