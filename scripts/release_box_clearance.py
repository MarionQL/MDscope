#!/usr/bin/env python3
"""Evaluate real MDscope release XTCs for protein–HS image-clearance evidence.

Read-only trajectory analysis; released systems may differ from AIM2 preparation.
No arbitrary conformer generation or approval of production box size.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from itertools import product
from scipy.spatial import cKDTree

def image_minimum(a, b, box, *, nonzero=True, z_image_only=False):
    """Exact nearest sampled points across relevant translated image cells.

    Coordinates may lie outside the primary box. The translation search is
    centered on their mean offset; includes adjacent translations, plus the
    nonzero immediate neighbors so nonzero image is never inadvertently skipped.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    box = np.asarray(box, dtype=float)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != 3 or b.shape[1] != 3:
        raise ValueError('Coordinates must be nonempty N x 3 arrays')
    if not len(a) or not len(b) or box.shape != (3,) or np.any(box <= 0):
        raise ValueError('Missing coordinates or invalid box')
    if not (np.isfinite(a).all() and np.isfinite(b).all() and np.isfinite(box).all()):
        raise ValueError('Nonfinite coordinates')
    # Build the KD tree over the larger set and query the smaller one.
    # This avoids repeatedly querying ~100k membrane atoms against a small
    # protein tree for every periodic translation.
    swapped = len(a) < len(b)
    if swapped:
        a, b = b, a
    center = np.rint((a.mean(axis=0)-b.mean(axis=0))/box).astype(int)
    # Include enough translations to cover the full spread, not just centers.
    extent = np.ceil((np.ptp(a, axis=0)+np.ptp(b, axis=0))/box).astype(int)+1
    candidates = set(product(*(range(int(c-e), int(c+e+1))
                              for c,e in zip(center,extent))))
    candidates.update(product((-1, 0, 1), repeat=3))
    tree = cKDTree(a)
    mean_a, mean_b = a.mean(axis=0), b.mean(axis=0)
    radii = (np.max(np.linalg.norm(a-mean_a, axis=1)) +
             np.max(np.linalg.norm(b-mean_b, axis=1)))
    best = float('inf')
    chosen = None
    best_pair = None
    for shift in sorted(candidates):
        if nonzero and shift == (0, 0, 0):
            continue
        if z_image_only and shift[2] == 0:
            continue
        # Lower bound on distance between bounding spheres permits skipping
        # many distant images, without missing a closer pair.
        center_distance = np.linalg.norm(mean_a-(mean_b+np.asarray(shift)*box))
        if center_distance-radii >= best:
            continue
        ds, neighbors = tree.query(b+np.asarray(shift)*box, k=1)
        local = int(np.argmin(ds))
        dist = float(ds[local])
        if dist < best:
            best, chosen = dist, list(shift)
            best_pair = [local, int(neighbors[local])] if swapped else [int(neighbors[local]),local]
    if chosen is None:
        raise ValueError('No periodic images selected')
    if swapped:
        chosen = [-int(x) for x in chosen]
    return {'distance_nm':best, 'shift':chosen, 'atom_indices':best_pair}



GLYCANS = "AGLCNA AIDOA BGLCA BGAL BXYL"
MEMBRANE = ("resname CHL1 POPC DPPC POPE POPI POPS PSM DOPC CER1 CER160 "
            "or chainID GM1a BGAL CER160")
AIM2_BOX_NM = np.array([20.55790, 20.55790, 24.71667])


def report_decision(rows, cutoff=1.2, margin=0.3):
    if not rows:
        raise ValueError("No frames analyzed")
    minima = np.asarray([min(r["assembly_image_nm"],r["membrane_z_image_nm"])
                         for r in rows])
    result = {"sampled_frames":len(rows),
        "first_time_ps":rows[0]["time_ps"],"last_time_ps":rows[-1]["time_ps"],
        "minimum_assembly_image_nm":min(r["assembly_image_nm"] for r in rows),
        "minimum_membrane_z_image_nm":min(r["membrane_z_image_nm"] for r in rows),
        "minimum_box_nm":[min(r[k] for r in rows) for k in
                          ("box_x_nm","box_y_nm","box_z_nm")],
        "contact_frames":int(np.count_nonzero(minima<=cutoff)),
        "near_contact_frames":int(np.count_nonzero(minima<=cutoff+margin)),
        "min_distance_quantiles_nm":dict(zip(("min","p05","median","p95","max"),
           map(float,np.quantile(minima,[0,.05,.5,.95,1]))))}
    result["conclusion"] = ("contact_observed_in_released_trajectory" if result["contact_frames"]
                           else "near_contact_observed_in_released_trajectory" if result["near_contact_frames"]
                           else "no_sampled_contacts_not_box_approval")
    return result


def assess(tpr, xtc, out, max_frames=50, cutoff=1.2, margin=.3):
    import MDAnalysis as mda
    from MDAnalysis.transformations import unwrap
    if max_frames<1: raise ValueError("max_frames must be >=1")
    u=mda.Universe(str(tpr),str(xtc))
    atoms=u.atoms
    if not hasattr(atoms,"chainIDs"):
        raise ValueError("Released TPR lacks chainIDs needed to map PROA/PROB")
    groups=[]
    for chain in ("PROA","PROB"):
        sel=u.select_atoms("chainID "+chain)
        heavy=sel[sel.masses>2.]
        proteins=sel.select_atoms("protein")
        glycans=sel.select_atoms("resname "+GLYCANS)
        if len(heavy)<100 or len(proteins)<100 or len(glycans)<100:
            raise ValueError(f"{chain}: expected protein and covalently associated HS; "
                             f"found {len(proteins)} protein and {len(glycans)} glycan atoms")
        groups.append((chain,heavy,glycans[glycans.masses>2.]))
    membrane=u.select_atoms("("+MEMBRANE+") and not name H*")
    membrane=membrane[membrane.masses>2.]
    if len(membrane)<1000:
        raise ValueError("Could not resolve intact membrane lipid selection")
    # The released XTC can store covalently bonded atoms across PBC.
    # Make each whole protein-HS molecule contiguous before asking whether
    # its distinct periodic COPY is within direct-space interaction range.
    proteins_hs=u.select_atoms("chainID PROA PROB")
    fragment_count=len(proteins_hs.fragments)
    if fragment_count!=2:
        raise ValueError(f"Expected exactly two covalent protein-HS fragments; found {fragment_count}")
    u.trajectory.add_transformations(unwrap(proteins_hs))
    # These native residue ranges come from MDscope's two-chain TPR selections.
    tm_ca={}; ecto_ca={}
    for chain,tm,ecto in (("PROA","150:170","1:149"),("PROB","420:440","270:419")):
        tm_ca[chain]=u.select_atoms(f"protein and chainID {chain} and name CA and resid {tm}")
        ecto_ca[chain]=u.select_atoms(f"protein and chainID {chain} and name CA and resid {ecto}")
        if len(tm_ca[chain])<5 or len(ecto_ca[chain])<10:
            raise ValueError(f"Missing verified TM/ectodomain CA ranges for {chain}")
    glycan_set=set(GLYCANS.split())
    n=len(u.trajectory)
    indices=np.unique(np.linspace(0,n-1,min(n,max_frames),dtype=int))
    rows=[]
    for i in indices:
        ts=u.trajectory[int(i)]
        dims=np.asarray(ts.dimensions,dtype=float)
        if dims.shape!=(6,) or np.any(np.abs(dims[3:]-90)>1e-2):
            raise ValueError("Orthorhombic trajectory box required")
        box=dims[:3]/10.
        if np.any(box<=0):
            raise ValueError("Nonpositive trajectory box")
        m=membrane.positions.astype(float)/10.
        positions=[heavy.positions.astype(float)/10. for _,heavy,_ in groups]
        pair_results=[]
        for a in range(len(positions)):
            for b in range(a,len(positions)):
                pair_results.append((groups[a][0],groups[b][0],
                    image_minimum(positions[a],positions[b],box)))
        mem_results=[(name,image_minimum(pos,m,box,z_image_only=True))
                     for (name,_,_),pos in zip(groups,positions)]
        pa=min(pair_results,key=lambda x:x[2]["distance_nm"])
        pm=min(mem_results,key=lambda x:x[1]["distance_nm"])
        # Record the limiting atom identities to distinguish HS contacts
        # from protein contacts at each periodic image crossing.
        a_atoms=next(group for name,group,_ in groups if name==pa[0])
        b_atoms=next(group for name,group,_ in groups if name==pa[1])
        aa=a_atoms[pa[2]["atom_indices"][0]]
        bb=b_atoms[pa[2]["atom_indices"][1]]
        def region(atom):
            return "HS" if atom.resname in glycan_set else "protein"
        categories=region(aa)+"-"+region(bb)
        ecto_reach=max(float(np.max(np.linalg.norm(ecto_ca[name].positions-
            tm_ca[name].positions.mean(axis=0),axis=1))/10.) for name,_,_ in groups)
        glycan_reach=max(float(np.max(np.linalg.norm(g.positions-
            tm_ca[name].positions.mean(axis=0),axis=1))/10.) for name,_,g in groups)
        rows.append({"frame":int(i),"time_ps":float(ts.time),
            "box_x_nm":float(box[0]),"box_y_nm":float(box[1]),"box_z_nm":float(box[2]),
            "assembly_image_nm":pa[2]["distance_nm"],
            "assembly_pair":pa[0]+"/"+pa[1],
            "limiting_atom_a":f"{aa.resname}:{aa.resid}:{aa.name}",
            "limiting_atom_b":f"{bb.resname}:{bb.resid}:{bb.name}",
            "limiting_regions":categories,
            "max_ectodomain_TM_distance_nm":ecto_reach,
            "max_HS_TM_distance_nm":glycan_reach,
            "assembly_shift":str(pa[2]["shift"]),
            "membrane_z_image_nm":pm[1]["distance_nm"],
            "membrane_pair":pm[0]+"/membrane","membrane_shift":str(pm[1]["shift"]),
            "max_glycan_z_span_nm":max(float(np.ptp(g.positions[:,2])/10.) for _,_,g in groups)})
    summary=report_decision(rows,cutoff,margin)
    summary.update({"source_tpr":Path(tpr).name,"source_xtc":Path(xtc).name,
        "trajectory_frames":n,"covalent_protein_HS_fragments":fragment_count,"periodic_unwrap":"MDAnalysis bonded-fragment whole-molecule", "time_window_ps":[rows[0]["time_ps"],rows[-1]["time_ps"]],
        "glycan_residue_counts":{name:len(g.residues) for name,_,g in groups},
        "max_glycan_z_span_nm":max(r["max_glycan_z_span_nm"] for r in rows),
        "max_ectodomain_TM_distance_nm":max(r["max_ectodomain_TM_distance_nm"] for r in rows),
        "max_HS_TM_distance_nm":max(r["max_HS_TM_distance_nm"] for r in rows),
        "contact_region_counts":{region:sum(r["limiting_regions"]==region and
            min(r["assembly_image_nm"],r["membrane_z_image_nm"])<=cutoff for r in rows)
            for region in ("HS-HS","HS-protein","protein-HS","protein-protein")},
        "AIM2_target_box_nm":AIM2_BOX_NM.tolist(),
        "sampled_box_relative_difference_from_AIM2":(np.array(summary["minimum_box_nm"])/AIM2_BOX_NM-1).tolist(),
        "scope":"Actual release box only. Released system can differ in sequence, signal "
                "peptide, HS glycosylation, membrane packing, amyloid absence, force field and "
                "thermal history. Full PROA/PROB include TM and HS; image minima do not "
                "quantify independently converged ectodomain states. Representative rep1 "
                "test trajectories do not bound rare configurations."})
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    with (out/"frames.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tpr",required=True);p.add_argument("--xtc",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--max-frames",type=int,default=50)
    p.add_argument("--cutoff",type=float,default=1.2)
    p.add_argument("--margin",type=float,default=.3)
    a=p.parse_args()
    assess(a.tpr,a.xtc,a.out,a.max_frames,a.cutoff,a.margin)
if __name__=="__main__":main()
