#!/usr/bin/env python3
"""Hypothetical lattice sizing from actual MDscope protein-HS conformers.

NOT a physical membrane resize. Shift each intact PROA/PROB into consistent
native-cell images using TM reference atoms, then assess periodic copies at
candidate box lengths. Coordinate perturbations, NPT re-equilibration and
changed membrane packing are NOT simulated. Never claim these are new MD runs.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from release_box_clearance import image_minimum

GLYCANS="AGLCNA AIDOA BGLCA BGAL BXYL"
MEMBRANE=("resname CHL1 POPC DPPC POPE POPI POPS PSM DOPC CER1 CER160 "
          "or chainID GM1a BGAL CER160")
# Candidate EQUILIBRATED lateral dimensions; later estimate BUILD sizes
# from the measured ~10.5% NPT lateral contraction in related reference runs.
LATERAL=(20.5579,22.,24.,26.,28.,30.,32.,36.)
VERTICAL=(24.71667,28.,32.,36.,40.,44.)
TARGET_Z=24.71667
CLEARANCE=1.5

def placement_shift(tm_xyz, native_box):
    """Choose a whole fragment image with TM CA centroid in the native cell."""
    center=np.mean(tm_xyz,axis=0)
    return -np.floor(center/native_box).astype(int)*native_box

def evaluate_candidates(positions, membrane, lateral, vertical, cutoff=CLEARANCE):
    # The same whole-assembly coordinates must be used for all trial lattices.
    # Keeping molecule structure unchanged is deliberate and is not new MD.
    out=[]
    for axis,values in (("lateral",lateral),("vertical",vertical)):
        for size in values:
            box=np.asarray([size,size,TARGET_Z] if axis=="lateral"
                           else [36.,36.,size],dtype=float)
            pairs=[]
            for a in range(len(positions)):
                for b in range(a,len(positions)):
                    d=image_minimum(positions[a],positions[b],box)
                    pairs.append((a,b,d))
            limiting=min(pairs,key=lambda x:x[2]["distance_nm"])
            pm=min((image_minimum(a,membrane,box,z_image_only=True)
                    for a in positions),key=lambda x:x["distance_nm"])
            cross=min((v for v in pairs if v[0]!=v[1]),key=lambda x:x[2]["distance_nm"])
            out.append({"axis":axis,"candidate_nm":float(size),
                "box_x_nm":float(box[0]),"box_y_nm":float(box[1]),"box_z_nm":float(box[2]),
                "min_assembly_image_nm":limiting[2]["distance_nm"],
                "min_cross_assembly_image_nm":cross[2]["distance_nm"],
                "min_same_assembly_image_nm":min(v[2]["distance_nm"] for v in pairs if v[0]==v[1]),
                "limiting_pair":f"PRO{chr(65+limiting[0])}-PRO{chr(65+limiting[1])}",
                "limiting_shift":str(limiting[2]["shift"]),
                "min_z_membrane_image_nm":pm["distance_nm"],
                "pass_1p5nm_snapshot":min(limiting[2]["distance_nm"],pm["distance_nm"])>cutoff})
    return out

def run(tpr,xtc,out,max_frames=40,lateral=None,vertical=None):
    import MDAnalysis as mda
    lateral=LATERAL if lateral is None else tuple(map(float,lateral))
    vertical=VERTICAL if vertical is None else tuple(map(float,vertical))
    if not lateral or not vertical or max_frames<1:
        raise ValueError("Nonempty candidate lists and positive frame count required")
    from MDAnalysis.transformations import unwrap
    u=mda.Universe(str(tpr),str(xtc))
    if not hasattr(u.atoms,"chainIDs"):raise ValueError("Missing release chainIDs")
    proteins_hs=u.select_atoms("chainID PROA PROB")
    if len(proteins_hs.fragments)!=2:
        raise ValueError("Expected exactly two covalent protein-HS fragments")
    u.trajectory.add_transformations(unwrap(proteins_hs))
    groups=[];tm=[]
    for chain,resids in (("PROA","150:170"),("PROB","420:440")):
        atoms=u.select_atoms("chainID "+chain);heavy=atoms[atoms.masses>2.]
        gly=atoms.select_atoms("resname "+GLYCANS)
        helix=u.select_atoms(f"protein and chainID {chain} and name CA and resid {resids}")
        if len(heavy)<100 or len(gly)<100 or len(helix)<5:
            raise ValueError("Unrecognized complete protein-HS reference: "+chain)
        groups.append(heavy);tm.append(helix)
    membrane=u.select_atoms("("+MEMBRANE+")")
    membrane=membrane[membrane.masses>2.]
    if len(membrane)<1000:raise ValueError("Could not identify membrane")
    indices=np.unique(np.linspace(0,len(u.trajectory)-1,min(len(u.trajectory),max_frames),dtype=int))
    records=[]
    for frame in indices:
        ts=u.trajectory[int(frame)]
        native=np.asarray(ts.dimensions[:3],dtype=float)/10.
        if np.any(np.abs(np.asarray(ts.dimensions[3:])-90)>1e-2):
            raise ValueError("Only orthorhombic unit cells supported")
        # Consistent anchor placement. Treat fragments as WHOLE, not atomwise wrapped.
        tm_centers=[np.mean(x.positions,axis=0)/10. for x in tm]
        coords=[g.positions.astype(float)/10. for g in groups]
        moves=[placement_shift(c,native) for c in tm_centers]
        # Both fragments are now in a consistent primary-native-cell image.
        positions=[p+d for p,d in zip(coords,moves)]
        membrane_pos=membrane.positions.astype(float)/10.
        for v in evaluate_candidates(positions,membrane_pos,lateral,vertical):
            records.append({"frame":int(frame),"time_ps":float(ts.time),
                            "native_lateral_x_nm":float(native[0]),
                            "native_lateral_y_nm":float(native[1]),
                            "native_z_nm":float(native[2]),**v})
        print(f"frame {int(frame)+1} sampled",flush=True)
    output=Path(out);output.mkdir(parents=True,exist_ok=True)
    with (output/"candidate-clearance.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]))
        writer.writeheader();writer.writerows(records)
    results=[]
    for axis,sizes in (("lateral",lateral),("vertical",vertical)):
        for size in sizes:
            vals=[r for r in records if r["axis"]==axis and r["candidate_nm"]==size]
            results.append({"axis":axis,"candidate_nm":size,"frames":len(vals),
                "worst_assembly_image_nm":min(r["min_assembly_image_nm"] for r in vals),
                "worst_cross_assembly_image_nm":min(r["min_cross_assembly_image_nm"] for r in vals),
                "worst_same_assembly_image_nm":min(r["min_same_assembly_image_nm"] for r in vals),
                "worst_z_membrane_image_nm":min(r["min_z_membrane_image_nm"] for r in vals),
                "frames_below_1p5nm":sum(not r["pass_1p5nm_snapshot"] for r in vals),
                "zero_observed_near_contacts":all(r["pass_1p5nm_snapshot"] for r in vals)})
    report={"source_tpr":Path(tpr).name,"source_xtc":Path(xtc).name,
       "sampled_frames":len(indices),"time_window_ps":[records[0]["time_ps"],records[-1]["time_ps"]],
       "native_lateral_min_nm":min(min(r["native_lateral_x_nm"],r["native_lateral_y_nm"]) for r in records),
       "candidate_lateral_at_fixed_z_nm":TARGET_Z,
       "candidate_vertical_at_fixed_xy_nm":36.,
       "comparison_cutoff_nm":CLEARANCE,"candidate_results":results,
       "limitations":"Hypothetical lattice at fixed historical reference conformations. "
         "These frames originate from a smaller equilibrated membrane with already "
         "present image interactions; replacing box vectors does not repack lipids, "
         "generate a real larger bilayer, remove conformational bias, or relax "
         "the membrane and protein. No isolated AIM2 box is certified. "
         "Source and AIM2 topology/amyloid differ; no time-independent confidence interval."}
    (output/"candidate-summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tpr",required=True);p.add_argument("--xtc",required=True)
    p.add_argument("--out",required=True);p.add_argument("--max-frames",type=int,default=40)
    p.add_argument("--lateral",nargs="+",type=float)
    p.add_argument("--vertical",nargs="+",type=float)
    a=p.parse_args()
    run(a.tpr,a.xtc,a.out,a.max_frames,a.lateral,a.vertical)
if __name__=="__main__":main()
