"""Reproducible batched lipid-contact workload; no timing pass thresholds.

Run each batch in a fresh process for comparable peak RSS:
python scripts/benchmark_lipid_contacts.py --batch-atoms 256
python scripts/benchmark_lipid_contacts.py --batch-atoms 4096
"""
import argparse
import hashlib
import json
import resource
import tempfile
import time
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from mdscope.observables.lipid_contacts import LipidContacts

if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--batch-atoms', type=int,default=4096)
    args=parser.parse_args(); rng=np.random.default_rng(2026)
    molecules, per_molecule=1024,12
    u=mda.Universe.empty(molecules*per_molecule,n_residues=molecules,
        atom_resindex=np.repeat(np.arange(molecules),per_molecule),trajectory=True)
    u.add_TopologyAttr('names', ['P']+['C']*(per_molecule-1) if molecules==1 else (['P']+['C']*(per_molecule-1))*molecules)
    u.add_TopologyAttr('resnames', np.where(np.arange(molecules)%2,'BBB','AAA'))
    u.add_TopologyAttr('resids',np.arange(molecules)+1)
    u.add_TopologyAttr('bonds',[(i*per_molecule+j,i*per_molecule+j+1) for i in range(molecules) for j in range(per_molecule-1)])
    xyz=np.repeat(rng.uniform([0,0,0],[160,160,0],(molecules,3)),per_molecule,axis=0)
    xyz[:,2]=np.repeat(np.where(np.arange(molecules)<molecules//2,30.,90.),per_molecule)
    xyz+=rng.normal(0,.5,xyz.shape)
    u.atoms.positions=xyz;u.dimensions=[160,160,160,90,90,90]
    spec={'selection':'all','species_selections':{'A':'resname AAA and name P','B':'resname BBB and name P'},'batch_atoms':args.batch_atoms}
    with tempfile.TemporaryDirectory() as tmp:
        started=time.perf_counter(); observer=LipidContacts(u,spec,tmp)
        for frame in range(3):observer.update(u.trajectory.ts,frame,float(frame))
        observer.flush(0);observer.finalize();seconds=time.perf_counter()-started
        table=pd.read_parquet(Path(tmp)/'extensions/lipid_contacts/matrix.parquet')
        digest=hashlib.sha256(pd.util.hash_pandas_object(table,index=False).values.tobytes()).hexdigest()
        print(json.dumps({'molecules':molecules,'atoms':len(u.atoms),'frames':3,'batch_atoms':args.batch_atoms,
            'seconds':seconds,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            'matrix_digest':digest},indent=2))
