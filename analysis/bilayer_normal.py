"""Local bilayer-normal angular deviations from a laboratory membrane normal.

Fit z = a*x + b*y + c to lipid headgroup anchors within each leaflet's
periodic neighbor shell, then theta = atan(sqrt(a²+b²)). Optimized sparse
neighbor accumulation avoids dense N_lipid² matrices.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance,minimize_vectors
from .config import ensure_dir,atomic_parquet,time_ps
from .lipid_geometry import anchor_groups,centers,aligned_leaflets
from .progress import progress

def fit_local_normals(xyz,upper,box,radius=15.,min_neighbors=5,normal=2):
    xyz=np.asarray(xyz,float)
    upper=np.asarray(upper,bool)
    if len(xyz)!=len(upper):raise ValueError("Mismatch between anchor and leaflet arrays")
    p=capped_distance(xyz,xyz,max_cutoff=radius,box=box,return_distances=False)
    if len(p)==0:return np.full(len(xyz),np.nan),np.zeros(len(xyz),int)
    i=p[:,0];j=p[:,1]
    choose=(i!=j)&(upper[i]==upper[j])
    i=i[choose];j=j[choose]
    if not len(i):return np.full(len(xyz),np.nan),np.zeros(len(xyz),int)
    delta=minimize_vectors(xyz[j]-xyz[i],box)
    tangent=[k for k in range(3) if k!=normal]
    x=delta[:,tangent[0]].astype(float)
    y=delta[:,tangent[1]].astype(float)
    z=delta[:,normal].astype(float)
    N=len(xyz)
    sums=lambda v:np.bincount(i,weights=v,minlength=N)
    n=np.bincount(i,minlength=N)
    safe=np.maximum(n,1)
    sx,sy,sz=sums(x),sums(y),sums(z)
    a=sums(x*x)-sx*sx/safe
    b=sums(y*y)-sy*sy/safe
    c=sums(x*y)-sx*sy/safe
    d=sums(x*z)-sx*sz/safe
    e=sums(y*z)-sy*sz/safe
    det=a*b-c*c
    valid=(n>=min_neighbors)&(det>1e-8)&np.isfinite(det)
    tilt=np.full(N,np.nan)
    slope_x=np.zeros(N);slope_y=np.zeros(N)
    slope_x[valid]=(b[valid]*d[valid]-c[valid]*e[valid])/det[valid]
    slope_y[valid]=(a[valid]*e[valid]-c[valid]*d[valid])/det[valid]
    tilt[valid]=np.degrees(np.arctan(np.hypot(slope_x[valid],slope_y[valid])))
    return tilt,n

def run_bilayer_normal(rep,cfg,repout):
    spec=cfg["analyses"]["bilayer_normal"]
    u=mda.Universe(str(rep.topology),str(rep.trajectory))
    tr=cfg["trajectory"]
    u.trajectory[0 if tr.get("start") is None else int(tr["start"])]
    species,indices,groups,upper=anchor_groups(u,spec)
    axis={"x":0,"y":1,"z":2}[str(spec.get("axis","z")).lower()]
    radius=float(spec.get("neighbor_radius_A",15.0))
    minimum=int(spec.get("min_neighbors",5))
    if radius<=0 or minimum<3:raise ValueError("neighbor_radius_A > 0 and min_neighbors >= 3")
    rows=[]; long=[]
    out=ensure_dir(repout/"extensions"/"bilayer_normal")
    save_detail=bool(spec.get("save_lipid_angles",False))
    frame_chunk=max(1,int(spec.get("angle_chunk_frames",100)))
    species_ids=np.unique(species)
    frames=u.trajectory[slice(tr.get("start"),tr.get("stop"),int(tr["stride"]))]
    with progress(frames,total=len(frames),desc=f"{rep.system}/{rep.label}: bilayer normals",
                  enabled=cfg.get("progress",{}).get("enabled",True)) as iterator:
        for k,ts in enumerate(iterator):
            t=time_ps(ts,k,cfg)
            xyz = centers(groups)
            current_upper, changes = aligned_leaflets(
                xyz[:, axis], float(u.dimensions[axis]), upper
            )
            angles,neighbors=fit_local_normals(xyz,current_upper,u.dimensions,radius,minimum,axis)
            frame={"time_ps":t,"bilayer_normal_valid_fraction":float(np.mean(np.isfinite(angles))),
                   "bilayer_normal_leaflet_changes": changes}
            for leaf,mask in (("upper",current_upper),("lower",~current_upper)):
                values=angles[mask];valid=values[np.isfinite(values)]
                frame[f"bilayer_normal_{leaf}_mean_deg"]=float(np.mean(valid)) if len(valid) else np.nan
                frame[f"bilayer_normal_{leaf}_sd_deg"]=float(np.std(valid)) if len(valid)>1 else np.nan
                for lipid in species_ids:
                    idx=mask&(species==lipid)
                    v=angles[idx];v=v[np.isfinite(v)]
                    frame[f"bilayer_normal_{lipid}_{leaf}_mean_deg"]=float(np.mean(v)) if len(v) else np.nan
                if save_detail:
                    long.extend({"time_ps":t,"resindex":int(indices[i]),"species":species[i],
                        "leaflet":leaf,"normal_deviation_deg":float(angles[i]),
                        "n_neighbors":int(neighbors[i])} for i in np.flatnonzero(mask))
            rows.append(frame)
            if save_detail and (k+1)%frame_chunk==0:
                atomic_parquet(ensure_dir(out/"per_lipid")/f"angles_part_{k//frame_chunk:06d}.parquet",
                               pd.DataFrame(long))
                long.clear()
    if long:
        atomic_parquet(ensure_dir(out/"per_lipid")/f"angles_part_{len(rows)//frame_chunk:06d}.parquet",
                       pd.DataFrame(long))
    atomic_parquet(out/"timeseries.parquet",pd.DataFrame(rows))
    return {"frames":len(rows),"n_lipids":len(groups)}
