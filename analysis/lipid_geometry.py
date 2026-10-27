"""Leaflet and periodic neighbor primitives shared by all membrane analyses.

Anchor selections must resolve one representative residue per physical lipid.
For full glycosphingolipids with several carbohydrate residues, select an
atom belonging to the ceramide anchor, not every carbohydrate residue.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from MDAnalysis.lib.distances import minimize_vectors

@dataclass
class AnchorTracker:
    """Fixed atom->lipid reduction, reuses NumPy bincount across frames."""
    atoms: object
    owner: np.ndarray
    counts: np.ndarray

    def __len__(self):
        return len(self.counts)

    def positions(self):
        xyz=self.atoms.positions.astype(float,copy=False)
        return np.column_stack([
            np.bincount(self.owner,weights=xyz[:,axis],minlength=len(self.counts))/self.counts
            for axis in range(3)
        ])

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
    indexes=np.concatenate([g.indices for g in atom_groups])
    owner=np.repeat(np.arange(len(atom_groups)),[len(g) for g in atom_groups])
    tracker=AnchorTracker(u.atoms[indexes],owner,
                          np.bincount(owner,minlength=len(atom_groups)))
    coords=tracker.positions()
    if np.any(~np.isfinite(coords)):raise ValueError("Non-finite anchor coordinates")
    normal={"x":0,"y":1,"z":2}[str(spec.get("axis","z")).lower()]
    box=np.asarray(u.dimensions,float)
    if box.shape!=(6,) or np.any(box[:3]<=0):
        raise ValueError("Lipid anchors require a periodic simulation box")
    leaf=classify_leaflets(coords[:,normal],float(box[normal]))
    return np.asarray(species,str),np.asarray(resindices,int),tracker,leaf

def centers(groups):
    return groups.positions() if isinstance(groups,AnchorTracker) else np.stack([g.center_of_geometry() for g in groups])

def planar_offsets(points,reference,box,normal=2):
    delta=minimize_vectors(np.asarray(points,float)-np.asarray(reference,float),np.asarray(box,float))
    return np.delete(delta,normal,axis=-1)

def lateral_distances(points,reference,box,normal=2):
    return np.linalg.norm(planar_offsets(points,reference,box,normal),axis=-1)
