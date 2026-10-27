"""Read-only TPR glycan-chain and covalent-attachment audit (no XTC required)."""
from pathlib import Path
import os
import numpy as np
import MDAnalysis as mda

root=Path(os.environ.get("MDSCOPE_REAL_DATA","real_test_data"))
for system,filename in (
    ("carb_sdc4","carb-sdc4-step7_production.tpr"),
    ("nocarb_sdc4","nocarb-sdc4-step7_production.tpr"),
):
    u=mda.Universe(str(root/filename))
    print(f"HS_TOPOLOGY {system}: atoms={len(u.atoms)} bonds={len(u.bonds)}",flush=True)
    residues, counts = np.unique(u.atoms.resnames, return_counts=True)
    print(f"ION_AUDIT {system}: all resnames={list(zip(residues.tolist(), counts.tolist()))}", flush=True)
    for selection in (
        "resname SOD", "resname CLA", "resname NA", "resname CL",
        "resname POT", "resname CAL", "name SOD", "name CLA",
        "name NA", "name CL", "name POT", "name CAL",
    ):
        print(f"ION_AUDIT {system}: {selection}={len(u.select_atoms(selection))}",flush=True)
    prot=u.select_atoms("protein")
    prot_ids=set(prot.indices.tolist())
    candidate=u.select_atoms(
        "resname AGLCNA AIDOA BGLCA BGAL BXYL "
        "and not chainID GM1a BGAL CER160"
    )
    print(f"HS_TOPOLOGY {system}: candidate_atoms={len(candidate)} "
          f"candidate_residues={len(candidate.residues)}",flush=True)
    if not len(candidate):
        continue
    bonds=np.asarray(u.bonds.indices,int)
    for chain in sorted(set(candidate.chainIDs)):
        sub=candidate[candidate.chainIDs==chain]
        gids=set(sub.indices.tolist())
        boundary=[(int(a),int(b)) for a,b in bonds
                  if (a in gids and b in prot_ids) or (b in gids and a in prot_ids)]
        attaches=[]
        for a,b in boundary:
            i=a if a in prot_ids else b
            atom=u.atoms[i]
            attaches.append((str(atom.chainID),int(atom.resid),str(atom.resname),str(atom.name)))
        print(f"HS_TOPOLOGY {system} chain={chain} n_residues={len(sub.residues)} "
              f"resnames={sorted(set(sub.resnames))} protein_bonds={len(boundary)} "
              f"attachment={sorted(set(attaches))}",flush=True)
