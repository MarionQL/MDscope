"""Leaflet-normalized lipid enrichment/depletion around a protein selection.

Raw counts and replicate-wise scalar series feed existing MDscope statistics.
A lipid is counted once using a single representative anchor per residue.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from ..config import atomic_json,atomic_parquet,ensure_dir,time_ps
from ..core.lipid_geometry import anchor_groups,centers,lateral_distances,periodic_center,aligned_leaflets,accessible_area_grid
from ..pipeline.progress import progress
from ..core.atom_selection import select_heavy_atoms

def shell_fold(observed,local_total,bulk_count,bulk_total):
    """Local mole fraction divided by bulk mole fraction in SAME leaflet."""
    if min(local_total,bulk_total)<=0 or bulk_count<=0:
        return np.nan
    return float((observed/local_total)/(bulk_count/bulk_total))

def accessible_plane_areas(obstacles, center, plane_height, box, inner, outer,
                           normal, exclusion_radius_A, slab_half_width_A, quadrature_order=64):
    """Geometric areas outside the existing protein-disk/leaflet-slab model.

    Deterministic midpoint quadrature integrates the periodic cell and annulus
    separately. This explicit probe model is not a lipid binding free energy,
    atom-specific van der Waals surface, or curved-membrane area correction.
    """
    box = np.asarray(box, float)
    obstacles = np.asarray(obstacles, float).reshape(-1, 3)
    center = np.asarray(center, float)
    n = quadrature_order
    if (box.shape != (6,) or not np.isfinite(box).all() or np.any(box[:3] <= 0)
            or not np.allclose(box[3:], 90., rtol=0, atol=1e-5)):
        raise ValueError("lipid area reference requires an orthorhombic planar box")
    if (normal not in (0, 1, 2) or center.shape != (3,) or not np.isfinite(center).all()
            or not np.isfinite(plane_height) or not np.isfinite(obstacles).all()):
        raise ValueError("lipid area reference requires finite coordinates and a valid normal")
    lateral = [i for i in range(3) if i != normal]
    if (not np.isfinite([inner, outer, exclusion_radius_A]).all()
            or not 0 <= inner < outer <= min(box[lateral])/2
            or not 0 < exclusion_radius_A < min(box[lateral])/2
            or not np.isfinite(slab_half_width_A) or not 0 < slab_half_width_A < box[normal]/2):
        raise ValueError("lipid area reference radii must fit the periodic cell")
    if type(n) is not int or not 16 <= n <= 256:
        raise ValueError("lipid area quadrature_order must be an integer from 16 to 256")
    fraction = (np.arange(n)+.5)/n
    x, y = np.meshgrid((fraction-.5)*box[lateral[0]],
                       (fraction-.5)*box[lateral[1]], indexing='ij')
    radius, angle = np.meshgrid(np.sqrt(inner**2+fraction*(outer**2-inner**2)),
                               fraction*2*np.pi, indexing='ij')
    areas = []
    for dx, dy, total in ((x, y, float(np.prod(box[lateral]))),
                         (radius*np.cos(angle), radius*np.sin(angle), np.pi*(outer**2-inner**2))):
        samples = np.column_stack((dx.ravel(), dy.ravel()))
        areas.append(accessible_area_grid(center[lateral], plane_height, obstacles,
            box, samples, np.sqrt(total/(n*n)), exclusion_radius_A, slab_half_width_A,
            normal=normal))
    return areas[0], areas[1]  # periodic-cell area, annular area (A^2)


def enrichment_frame(species,upper,positions,protein_center,box,inner,outer,normal,
                     *, categories=None, obstacles=None, area_reference=None):
    """Return long-form per-leaflet species metrics and raw denominators."""
    distance=lateral_distances(positions,protein_center,box,normal)
    local=(distance>=inner)&(distance<outer)
    categories=np.unique(species) if categories is None else sorted(categories)
    rows=[]
    for leaf,mask in (("upper",upper),("lower",~upper)):
        bulk_total=int(mask.sum())
        near=mask&local
        n_local=int(near.sum())
        area = {}
        if area_reference:
            # Circular mean preserves the leaflet plane across periodic wraps.
            phase = np.mean(np.exp(2j*np.pi*positions[mask, normal]/box[normal])) if bulk_total else 0j
            if abs(phase) < 1e-5:
                raise ValueError("lipid area reference needs an unambiguous populated leaflet plane")
            height = (np.angle(phase) % (2*np.pi))*box[normal]/(2*np.pi)
            full, shell = accessible_plane_areas(obstacles, protein_center, height, box,
                inner, outer, normal, float(area_reference["exclusion_radius_A"]),
                float(area_reference["slab_half_width_A"]), area_reference.get("quadrature_order", 64))
            area = dict(leaflet_accessible_area_A2=full, shell_accessible_area_A2=shell,
                        leaflet_accessible_fraction=full/float(np.prod(np.delete(np.asarray(box[:3], float), normal))),
                        shell_accessible_fraction=shell/(np.pi*(outer**2-inner**2)))
        for name in categories:
            smask=(species==name)
            count=int((near&smask).sum())
            bulk=int((mask&smask).sum())
            fold=shell_fold(count,n_local,bulk,bulk_total)
            rows.append({"leaflet":leaf,"species":name,
                         "count_local":count,"count_bulk":bulk,
                         "n_local":n_local,"n_bulk":bulk_total,
                         "local_fraction":count/n_local if n_local else np.nan,
                         "bulk_fraction":bulk/bulk_total if bulk_total else np.nan,
                         "enrichment_fold":fold,
                         "depleted":bool(fold<1.) if np.isfinite(fold) else False,
                         **area})
            if area_reference:
                rows[-1]["area_density_ratio"] = (count*full/(bulk*shell)
                    if bulk > 0 and full > 0 and shell > 0 else np.nan)
                rows[-1]["area_reference_status"] = (
                    "unavailable_absent_species" if not bulk else
                    "unavailable_zero_accessible_area" if min(full, shell) <= 0 else "available")
    return rows

def run_lipid_enrichment(rep,cfg,repout):
    spec=cfg["analyses"]["lipid_enrichment"]
    u=mda.Universe(str(rep.topology),str(rep.trajectory))
    tr=cfg["trajectory"]
    start=tr.get("start");stop=tr.get("stop");stride=int(tr.get("stride",1))
    u.trajectory[0 if start is None else int(start)]
    species,indices,groups,leaf=anchor_groups(u,spec)
    protein=u.select_atoms(str(spec.get("protein_selection") or cfg["selections"]["protein"]))
    if not len(protein):raise ValueError("lipid_enrichment.protein_selection is empty")
    obstacles = select_heavy_atoms(protein, policy="mechanistic") if spec.get("area_reference") else None
    outer=float(spec.get("outer_radius_A",15.0))
    inner=float(spec.get("inner_radius_A",0.0))
    if not 0<=inner<outer:raise ValueError("Require 0 <= inner_radius_A < outer_radius_A")
    normal={"x":0,"y":1,"z":2}[str(spec.get("axis","z")).lower()]
    rows=[];tables=[]
    frames=u.trajectory[slice(start,stop,stride)]
    with progress(frames,total=len(frames),desc=f"{rep.system}/{rep.label}: lipid enrichment",
                  enabled=cfg.get("progress",{}).get("enabled",True)) as iterator:
        for k,ts in enumerate(iterator):
            t=time_ps(ts,k,cfg)
            box=np.asarray(u.dimensions,float)
            anchors = centers(groups)
            current_leaf, changes = aligned_leaflets(anchors[:, normal], float(box[normal]), leaf)
            # A Cartesian center_of_geometry is incorrect for a protein dimer
            # straddling a periodic boundary. Use a circular lattice center.
            center = periodic_center(protein.positions, box)
            vals=enrichment_frame(species,current_leaf,anchors,center,box,inner,outer,normal,
                categories=spec["species_selections"],
                obstacles=obstacles.positions if obstacles is not None else None,
                area_reference=spec.get("area_reference"))
            row={"time_ps":t, "lipid_enrichment_leaflet_changes": changes}
            for value in vals:
                name=f"lipid_enrichment_{value['species']}_{value['leaflet']}"
                row[name+"_fold"]=value["enrichment_fold"]
                row[name+"_local_fraction"]=value["local_fraction"]
                row[name+"_contact_count"]=value["count_local"]
                if spec.get("area_reference"):
                    row[name+"_area_density_ratio"] = value["area_density_ratio"]
                tables.append({"time_ps":t,**value})
            rows.append(row)
    out=ensure_dir(repout/"extensions"/"lipid_enrichment")
    atomic_parquet(out/"timeseries.parquet",pd.DataFrame(rows))
    atomic_parquet(out/"counts.parquet",pd.DataFrame(tables))
    atomic_parquet(out/"selection_audit.parquet",pd.DataFrame({
        "species":species,"resindex":indices,"leaflet":np.where(leaf,"upper","lower")}))
    write_area_reference_metadata(out, spec, spec.get("protein_selection") or cfg["selections"]["protein"])
    return {"frames":len(rows),"n_lipids":len(groups),"n_species":len(np.unique(species))}


def write_area_reference_metadata(out, spec, protein_selection):
    """Persist the optional geometric definition shared by both schedulers."""
    if spec.get("area_reference"):
        atomic_json(Path(out)/"area_reference.json", {
            **spec["area_reference"], "quadrature_order": spec["area_reference"].get("quadrature_order", 64),
            "protein_selection": protein_selection,
            "hydrogen_policy": "mechanistic heavy atoms",
            "model": "equal-radius projected heavy-atom disks within a slab around each circular-mean anchor plane",
            "quadrature": "deterministic midpoint grids in cell lateral coordinates and annular r^2/theta",
            "limitation": "Geometric probe model; validate radius, resolution and planar approximation; not affinity or atom-specific VDW area"})


def summarize_lipid_composition(repout):
    """Offline population inventory from already counted molecular anchors."""
    directory = Path(repout)/"extensions/lipid_enrichment"
    source = directory/"counts.parquet"
    if not source.exists():
        return
    counts = pd.read_parquet(source)
    if counts.empty:
        return
    keys = ["leaflet", "species"]
    if counts.duplicated(["time_ps", *keys]).any():
        raise ValueError(f"Duplicate lipid composition observations: {source}")
    summary = counts.assign(available=counts.count_bulk > 0).groupby(keys, as_index=False).agg(
        sampled_frames=("time_ps", "size"),
        mean_molecules=("count_bulk", "mean"),
        minimum_molecules=("count_bulk", "min"), maximum_molecules=("count_bulk", "max"),
        mean_leaflet_fraction=("bulk_fraction", "mean"),
        fraction_frames_available=("available", "mean"),
        mean_local_molecules=("count_local", "mean"))
    atomic_parquet(directory/"composition_summary.parquet", summary)
