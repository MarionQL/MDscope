"""Bond-aware intermolecular lipid–lipid hydrogen bonds.

A glycosphingolipid may span several CHARMM-GUI residues. Intermolecular
exclusion must therefore use bonded connected components (whole molecules),
not residue IDs. Chain IDs identify GM1a/BGAL/CER160 species in the
published Syndecan-4 TPRs, but are *not* unique molecule identifiers.
"""
from __future__ import annotations

from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from ..config import ensure_dir, atomic_parquet
from .interaction_tables import atom_label, count_timeseries


GLYCOLIPID_CHAIN_IDS = ("GM1a", "BGAL", "CER160")


def ordered_pair(a, b):
    return "__".join(sorted((str(a), str(b))))


def lipid_molecule_map(universe, selected_atoms):
    """Return per-atom bonded molecule IDs and lipid species.

    Graph nodes are the selected lipid atoms. Bonds come from the full TPR,
    so inter-residue covalent links within GM1a/BGAL are retained. An
    incomplete selection can split molecules: select *all* atoms of each
    lipid type, not just headgroup atoms, when running this analysis.
    """
    atoms = selected_atoms
    n = len(universe.atoms)
    local = np.full(n, -1, dtype=np.int32)
    local[atoms.indices] = np.arange(len(atoms), dtype=np.int32)
    try:
        bonds = np.asarray(universe.bonds.indices, dtype=np.int32)
    except (AttributeError, mda.exceptions.NoDataError) as exc:
        raise ValueError(
            "Lipid molecule identification requires topology bonds (e.g. TPR); "
            "cannot reliably distinguish multi-residue glycolipids without bonds"
        ) from exc
    a = local[bonds[:, 0]]
    b = local[bonds[:, 1]]
    valid = (a >= 0) & (b >= 0)
    a, b = a[valid], b[valid]
    graph = coo_matrix(
        (np.ones(len(a), dtype=np.int8), (a, b)),
        shape=(len(atoms), len(atoms)),
    ).tocsr()
    _, component = connected_components(graph, directed=False)
    molecule_ids = np.full(n, -1, dtype=np.int32)
    molecule_ids[atoms.indices] = component
    species = np.asarray(universe.atoms.resnames, dtype=str).copy()
    try:
        chainids = np.asarray(universe.atoms.chainIDs, dtype=str)
        glycolipid = np.isin(chainids, GLYCOLIPID_CHAIN_IDS)
        species[glycolipid] = chainids[glycolipid]
    except (AttributeError, mda.exceptions.NoDataError):
        pass
    return molecule_ids, species


def tabulate_events(
    event_array, atoms, frames, times, exclude_same_residue=True,
    molecule_ids=None, atom_species=None, exclude_same_molecule=True,
):
    cols = [
        "time_ps", "donor", "hydrogen", "acceptor",
        "donor_molecule_id", "acceptor_molecule_id",
        "donor_species", "acceptor_species", "species_pair", "scope",
        "interaction_class", "distance_A", "angle_deg",
    ]
    data = np.asarray(event_array)
    if data.size == 0:
        return pd.DataFrame(columns=cols)
    donor = data[:, 1].astype(int)
    hydrogen = data[:, 2].astype(int)
    acceptor = data[:, 3].astype(int)
    kept = np.ones(len(data), dtype=bool)
    if exclude_same_residue:
        kept &= atoms[donor].resindices != atoms[acceptor].resindices
    if exclude_same_molecule and molecule_ids is not None:
        # A single molecule can contain CER160 plus many glycan residues.
        kept &= molecule_ids[donor] != molecule_ids[acceptor]
    data = data[kept]
    donor, hydrogen, acceptor = donor[kept], hydrogen[kept], acceptor[kept]
    frame_map = dict(zip(np.asarray(frames, int), np.asarray(times, float)))
    if atom_species is None:
        sp1 = np.asarray(atoms[donor].resnames, str)
        sp2 = np.asarray(atoms[acceptor].resnames, str)
    else:
        sp1, sp2 = atom_species[donor], atom_species[acceptor]
    pairs = [ordered_pair(a, b) for a, b in zip(sp1, sp2)]
    return pd.DataFrame({
        "time_ps": [frame_map[int(f)] for f in data[:, 0]],
        "donor": [atom_label(a) for a in atoms[donor]],
        "hydrogen": [atom_label(a) for a in atoms[hydrogen]],
        "acceptor": [atom_label(a) for a in atoms[acceptor]],
        "donor_molecule_id": molecule_ids[donor] if molecule_ids is not None else atoms[donor].resindices,
        "acceptor_molecule_id": molecule_ids[acceptor] if molecule_ids is not None else atoms[acceptor].resindices,
        "donor_species": sp1, "acceptor_species": sp2, "species_pair": pairs,
        "scope": ["intraspecies" if a == b else "interspecies" for a, b in zip(sp1, sp2)],
        "interaction_class": pairs,
        "distance_A": data[:, 4].astype(float),
        "angle_deg": data[:, 5].astype(float),
    }, columns=cols)


def run_lipid_hydrogen_bonds(rep, cfg, repout):
    spec = cfg["analyses"]["lipid_hydrogen_bonds"]
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    selection = str(spec.get("lipid_selection") or (
        "(resname POPC POPE POPI POPS DOPC DPPC PSM CHL1) "
        "or (chainID GM1a BGAL CER160)"
    ))
    atoms = universe.select_atoms(selection)
    if not len(atoms):
        raise ValueError("lipid_hydrogen_bonds.lipid_selection matched no atoms")
    molecule_ids, atom_species = lipid_molecule_map(universe, atoms)
    hydrogen_sel = spec.get("hydrogens_sel")
    acceptor_sel = spec.get("acceptors_sel")
    if hydrogen_sel is None or acceptor_sel is None:
        guesser = HydrogenBondAnalysis(
            universe, donors_sel=None, hydrogens_sel="name __MDSCOPE_NONE__",
            acceptors_sel="name __MDSCOPE_NONE__", update_selections=False
        )
        if hydrogen_sel is None:
            hydrogen_sel = guesser.guess_hydrogens(selection)
        if acceptor_sel is None:
            acceptor_sel = guesser.guess_acceptors(selection)
    if not hydrogen_sel or not acceptor_sel:
        raise ValueError("No lipid donor hydrogens/acceptors inferred; supply explicit selectors")
    hydrogen_sel = f"({hydrogen_sel}) and ({selection})"
    acceptor_sel = f"({acceptor_sel}) and ({selection})"
    donor_sel = spec.get("donors_sel")
    if donor_sel is not None:
        donor_sel = f"({donor_sel}) and ({selection})"
    hba = HydrogenBondAnalysis(
        universe, donors_sel=donor_sel, hydrogens_sel=hydrogen_sel,
        acceptors_sel=acceptor_sel,
        d_h_cutoff=float(spec.get("d_h_cutoff_A", 1.2)),
        d_a_cutoff=float(spec.get("d_a_cutoff_A", 3.0)),
        d_h_a_angle_cutoff=float(spec.get("angle_cutoff_deg", 150.0)),
        update_selections=bool(spec.get("update_selections", False)),
    )
    tr = cfg["trajectory"]
    try:
        hba.run(start=tr.get("start"), stop=tr.get("stop"),
                step=int(tr.get("stride", 1)))
    except Exception as exc:
        raise RuntimeError(
            "Lipid H-bond analysis requires explicit H atoms/bond connectivity; verify selections"
        ) from exc
    return write_lipid_hydrogen_bonds(universe, cfg, repout, hba.results.hbonds, hba.frames, hba.times, atoms, molecule_ids, atom_species)


def write_lipid_hydrogen_bonds(universe, cfg, repout, raw, frames, times, atoms, molecule_ids, atom_species):
    spec = cfg["analyses"]["lipid_hydrogen_bonds"]
    times = np.asarray(times, float)
    if cfg.get("time", {}).get("dt_ps") is not None:
        times = cfg["time"].get("origin_ps", 0.) + np.arange(len(times)) * float(cfg["time"]["dt_ps"])
    raw = np.asarray(raw)
    if raw.size:
        donors = raw[:, 1].astype(int)
        acceptors = raw[:, 3].astype(int)
        same_molecule = molecule_ids[donors] == molecule_ids[acceptors]
        same_residue = universe.atoms[donors].resindices == universe.atoms[acceptors].resindices
    else:
        same_molecule = np.array([], bool)
        same_residue = np.array([], bool)
    df = tabulate_events(
        raw, universe.atoms, frames, times,
        exclude_same_residue=bool(spec.get("exclude_same_residue", True)),
        molecule_ids=molecule_ids, atom_species=atom_species,
        exclude_same_molecule=bool(spec.get("exclude_same_molecule", True)),
    )
    out = ensure_dir(Path(repout) / "interactions" / "lipid_hydrogen_bonds")
    atomic_parquet(out / "events.parquet", df)
    counts = count_timeseries(
        df, times, "lipid_hydrogen_bond", "interaction_class",
        ["donor_molecule_id", "acceptor_molecule_id", "donor", "hydrogen",
         "acceptor", "interaction_class"],
    )
    atomic_parquet(out / "count_timeseries.parquet", counts)
    audit = pd.DataFrame([{
        "raw_hbond_events": int(len(raw)),
        "same_residue_events": int(same_residue.sum()),
        "same_molecule_events": int(same_molecule.sum()),
        "cross_residue_same_molecule_events": int((same_molecule & ~same_residue).sum()),
        "intermolecular_events": int(len(df)),
        "selected_atoms": int(len(atoms)),
        "selected_molecules": int(len(np.unique(molecule_ids[atoms.indices]))),
    }])
    atomic_parquet(out / "molecule_audit.parquet", audit)
    return {
        "frames": len(times), "events": len(df),
        "species_pairs": df["species_pair"].nunique(),
        "same_molecule_removed": int((same_molecule & ~same_residue).sum()),
    }
