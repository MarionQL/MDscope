"""Common tools for neuronal membrane observables. PBC distances in Angstrom."""
from __future__ import annotations
import re
import numpy as np
from MDAnalysis.lib.distances import capped_distance, distance_array
from .core.atom_selection import select_heavy_atoms

AA = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E",
"GLN":"Q","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M",
"PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def box_of(u):
    b = np.asarray(u.dimensions, dtype=float)
    return b if b.shape==(6,) and np.all(b[:3]>0) else None

def heavy(ag):
    return select_heavy_atoms(ag, policy="biophysical")

def atom_contacts(a,b,cutoff,box=None):
    if len(a)==0 or len(b)==0: return np.empty((0,2),int)
    return capped_distance(a.positions,b.positions,max_cutoff=float(cutoff),box=box,return_distances=False)

def residue_pairs(a,b,cutoff,box=None):
    pairs=atom_contacts(a,b,cutoff,box)
    if not len(pairs):return np.empty((0,2),int)
    return np.unique(np.stack((a.resindices[pairs[:,0]],b.resindices[pairs[:,1]]),axis=1),axis=0)

def lateral_distance(a,b,box,axis=2):
    """Planar distances for orthorhombic boxes; fails explicitly for triclinic."""
    if not np.allclose(box[3:],90,atol=0.05):
        raise ValueError("Lateral planar analysis requires orthorhombic box")
    dims=[i for i in range(3) if i!=axis]
    d=np.asarray(a)[:,None,dims]-np.asarray(b)[None,:,dims]
    lengths=np.asarray(box[:3])[dims]
    d-=np.rint(d/lengths)*lengths
    return np.linalg.norm(d,axis=-1)

def sequence_motifs(residues,patterns):
    """Return chain-local sequence matches; never match across chain boundaries."""
    chains={}
    for r in residues:
        key=str(getattr(r,"chainID","") or getattr(r,"segid","") or "unknown")
        chains.setdefault(key,[]).append(r)
    out=[]
    for chain,rs in chains.items():
        seq="".join(AA.get(str(r.resname).upper(),"X") for r in rs)
        for name,pat in patterns.items():
            for hit in re.finditer(pat,seq):
                out.append({"chain":chain,"motif":name,"sequence":hit.group(),
                            "start_resid":int(rs[hit.start()].resid),
                            "end_resid":int(rs[hit.end()-1].resid),
                            "resindices":[int(r.ix) for r in rs[hit.start():hit.end()]]})
    return out

def leaflet_groups(atoms,axis=2):
    """Assign lipid residues by circular z position on the first frame."""
    residues=list(atoms.residues)
    centers=np.array([r.atoms.center_of_geometry() for r in residues])
    z=centers[:,axis]; L=float(atoms.universe.dimensions[axis])
    phase=np.exp(2j*np.pi*z/L)
    center=(np.angle(phase.mean())/(2*np.pi)*L)%L
    wrapped=(z-center+L/2)%L-L/2
    # For a centered bilayer leaflets are opposing z bands; use largest sorted gap
    vals=np.sort(wrapped); gaps=np.diff(vals)
    if len(gaps)==0:raise ValueError("Need at least two lipid residues for leaflets")
    cut=(vals[np.argmax(gaps)]+vals[np.argmax(gaps)+1])/2
    return {int(r.ix):("upper" if w>cut else "lower") for r,w in zip(residues,wrapped)}
