"""Leaflet and periodic neighbor primitives shared by all membrane analyses.

Anchor selections must resolve one representative residue per physical lipid.
For full glycosphingolipids with several carbohydrate residues, select an
atom belonging to the ceramide anchor, not every carbohydrate residue.
"""
from __future__ import annotations
import numpy as np
from MDAnalysis.lib.distances import minimize_vectors

def classify_leaflets(z, length):
    """Two-cluster circular partition robust to a bilayer crossing periodic z."""
    z=np.asarray(z,float)
    if len(z)<4 or length<=0:raise ValueError("Need >=4 lipid anchors and positive box length")
    angle=(2*np.pi*z/length)%(2*np.pi)
    unit=np.stack([np.cos(angle),np.sin(angle)],axis=1)
    seed0=unit[0]
    seed1=unit[np.argmin(unit@seed0)]
    centers=np.stack([seed0,seed1])
    groups=np.zeros(len(z),bool)
    for _ in range(30):
        assignment=np.argmax(unit@centers.T,axis=1)
        if np.any(np.bincount(assignment,minlength=2)<2):
            raise ValueError("Cannot resolve two populated bilayer leaflets")
        fresh=np.stack([unit[assignment==i].mean(axis=0) for i in range(2)])
        fresh/=np.linalg.norm(fresh,axis=1)[:,None]
        if np.allclose(fresh,centers,rtol=0,atol=1e-7):break
        centers=fresh
    # Assign conventional upper/lower using shortest arc from cluster 0 -> 1.
    angle_centers=np.arctan2(centers[:,1],centers[:,0])
    delta=(angle_centers[1]-angle_centers[0]+np.pi)%(2*np.pi)-np.pi
    upper_label=1 if delta>0 else 0
    return assignment==upper_label

def anchor_groups(u,spec):
    """Return indexed per-species lipid anchors and first-frame leaflet labels."""
    selections=spec.get("species_selections") or {}
    if not selections:raise ValueError("Define species_selections: {POPC: 'resname POPC and name P', ...}")
    species=[];atom_groups=[];resindices=[]
    for name,selection in selections.items():
        atoms=u.select_atoms(str(selection))
        if not len(atoms):
            raise ValueError(f"Lipid species {name!r} selection matched zero atoms: {selection}")
        for residue in atoms.residues:
            part=atoms.intersection(residue.atoms)
            species.append(str(name));resindices.append(int(residue.ix));atom_groups.append(part)
    if len(resindices)!=len(set(resindices)):
        raise ValueError("species_selections overlap on lipid residues")
    coords=np.stack([g.center_of_geometry() for g in atom_groups])
    if np.any(~np.isfinite(coords)):raise ValueError("Non-finite anchor coordinates")
    normal={"x":0,"y":1,"z":2}[str(spec.get("axis","z")).lower()]
    box=np.asarray(u.dimensions,float)
    if box.shape!=(6,) or np.any(box[:3]<=0):
        raise ValueError("Lipid anchors require a periodic simulation box")
    leaf=classify_leaflets(coords[:,normal],float(box[normal]))
    return np.asarray(species,str),np.asarray(resindices,int),atom_groups,leaf

def centers(groups):
    return np.stack([g.center_of_geometry() for g in groups])

def planar_offsets(points,reference,box,normal=2):
    delta=minimize_vectors(np.asarray(points,float)-np.asarray(reference,float),np.asarray(box,float))
    return np.delete(delta,normal,axis=-1)

def lateral_distances(points,reference,box,normal=2):
    return np.linalg.norm(planar_offsets(points,reference,box,normal),axis=-1)
