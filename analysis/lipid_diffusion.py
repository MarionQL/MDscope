"""Species-resolved, multi-time-origin lateral lipid diffusion.

Lipid COMs are unwrapped incrementally; translational membrane drift is
subtracted. MSDs use multiple origins and no framewise inferential statistics.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import minimize_vectors
from .biophysical_common import box_of
from .config import atomic_parquet,ensure_dir,time_ps

def multi_origin_msd(x, lag_frames):
    """Average over all starting times and all lipid residues, Å²."""
    x=np.asarray(x,float)
    return np.array([np.mean(np.sum((x[k:]-x[:-k])**2,axis=-1)) for k in lag_frames])

def diffusion_fit(lag_ps,msd,begin=.25,end=.75):
    """Return lateral D in nm²/ns, slope=4D for 2D Brownian motion."""
    t=np.asarray(lag_ps,float);m=np.asarray(msd,float)
    a=t.min()+begin*np.ptp(t); b=t.min()+end*np.ptp(t)
    keep=(t>=a)&(t<=b)&np.isfinite(m)
    if keep.sum()<3:return np.nan
    return float(np.polyfit(t[keep],m[keep],1)[0]*2.5)

def run_lipid_diffusion(rep,cfg,repout):
    import MDAnalysis as mda
    spec=cfg["analyses"].get("lipid_diffusion",{})
    u=mda.Universe(str(rep.topology),str(rep.trajectory))
    head=u.select_atoms(spec.get("selection","(resname POPC POPE POPI POPS DOPC DPPC PSM and name P) or (resname CHL1 CER160 CER1 BGAL GM1a and name O*)"))
    if not len(head):raise ValueError("No tracked lipid atoms for diffusion")
    groups=[head & r.atoms for r in head.residues]
    names=np.asarray([r.resname for r in head.residues],str)
    t=[];coords=[];prev=None;displacement=None
    tr=cfg["trajectory"]
    axis={"x":0,"y":1,"z":2}[str(spec.get("axis","z")).lower()]
    lateral=[k for k in (0,1,2) if k!=axis]
    for i,ts in enumerate(u.trajectory[slice(tr.get("start"),tr.get("stop"),int(tr.get("stride",1)))]):
        b=box_of(u)
        if b is None:raise ValueError("Diffusion requires valid periodic box")
        centers=np.stack([a.center_of_geometry() for a in groups])
        if prev is None:
            displacement=np.zeros_like(centers)
        else:
            delta=minimize_vectors(centers-prev,b)
            displacement+=delta
        prev=centers
        t.append(time_ps(ts,i,cfg))
        coords.append(displacement[:,lateral].copy())
    if len(coords)<8:raise ValueError("Diffusion needs >=8 analyzed frames")
    x=np.stack(coords)
    x-=x.mean(axis=1,keepdims=True) # remove global lipid drift
    times=np.asarray(t,float)
    step=np.median(np.diff(times))
    if step<=0 or not np.allclose(np.diff(times),step,rtol=1e-3,atol=1e-3):
        raise ValueError("Diffusion requires approximately regularly sampled trajectory")
    max_lag=int(min(len(x)//3,int(spec.get("max_lag_frames",len(x)//3))))
    lag=np.unique(np.linspace(1,max_lag,min(70,max_lag)).astype(int))
    rows=[];summ=[]
    for name in np.unique(names):
        arr=x[:,names==name,:]
        if not len(arr[0]):continue
        msd=multi_origin_msd(arr,lag)
        lagtime=lag*step
        D=diffusion_fit(lagtime,msd,float(spec.get("fit_start_fraction",.25)),float(spec.get("fit_end_fraction",.75)))
        summ.append({"species":str(name),"diffusion_nm2_ns":D,"n_lipids":arr.shape[1],
                     "fit_start_fraction":spec.get("fit_start_fraction",.25),"fit_end_fraction":spec.get("fit_end_fraction",.75)})
        rows.extend({"species":str(name),"lag_ps":float(tt),"msd_lateral_A2":float(mm)} for tt,mm in zip(lagtime,msd))
    out=ensure_dir(repout/"extensions"/"lipid_diffusion")
    atomic_parquet(out/"msd_by_species.parquet",pd.DataFrame(rows))
    atomic_parquet(out/"coefficients.parquet",pd.DataFrame(summ))
    if spec.get("plots",True):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(7.3,4.8))
        for name,part in pd.DataFrame(rows).groupby("species"):
            ax.plot(part["lag_ps"]/1000,part["msd_lateral_A2"]*.01,label=name)
        ax.set(xlabel="Lag time (ns)",ylabel="Lateral MSD (nm²)",title="Multi-origin lateral diffusion")
        ax.legend(frameon=False,ncol=2)
        fig.tight_layout()
        fig.savefig(out/"diffusion_msd.png",dpi=250)
        plt.close(fig)
    return {"species":len(summ)}
