"""Molecule-resolved lipid neighbor matrices from bounded periodic searches.

One anchor residue per bonded molecule defines species and leaflet. Output
pair counts are undirected, repeated for AB/BA rows to expose directional
neighbors per A molecule. No atom identities are persisted.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance
from ..core.lipid_geometry import anchor_groups, aligned_leaflets
from ..core.atom_selection import select_heavy_atoms
from ..config import atomic_parquet, atomic_json, ensure_dir


def contact_matrix(species, upper, molecule_pairs, categories=None):
    """Finite-population random-mixing reference, including same-type pairs."""
    species, upper = np.asarray(species, str), np.asarray(upper, bool)
    pairs = np.asarray(molecule_pairs, int).reshape(-1, 2)
    if len(pairs):
        if pairs.min() < 0 or pairs.max() >= len(species):
            raise ValueError('Molecule pair outside anchor population')
        pairs = np.unique(np.sort(pairs[pairs[:, 0] != pairs[:, 1]], axis=1), axis=0)
    rows = []
    categories = sorted(set(species) if categories is None else categories)
    for leaflet, mask in [('upper', upper), ('lower', ~upper)]:
        kept = pairs[mask[pairs].all(axis=1)] if len(pairs) else pairs
        total_pairs = len(kept)
        n = int(mask.sum())
        possible_total = n*(n-1)/2
        for a in categories:
            na = int(np.count_nonzero(mask & (species == a)))
            for b in categories:
                nb = int(np.count_nonzero(mask & (species == b)))
                count = int(np.count_nonzero(
                    ((species[kept[:, 0]] == a) & (species[kept[:, 1]] == b))
                    | ((species[kept[:, 0]] == b) & (species[kept[:, 1]] == a)))) if len(kept) else 0
                possible = na*(na-1)/2 if a == b else na*nb
                expected = total_pairs*possible/possible_total if possible_total else np.nan
                rows.append(dict(leaflet=leaflet, species_a=a, species_b=b,
                    molecules_a=na, molecules_b=nb, pair_count_undirected=count,
                    possible_pairs=possible, total_contact_pairs=total_pairs,
                    neighbors_per_molecule=(count*(2 if a == b else 1)/na if na else np.nan),
                    enrichment=(count/expected if expected > 0 else np.nan)))
    return pd.DataFrame(rows)


class LipidContacts:
    def __init__(self, universe, spec, repout):
        self.u, self.spec = universe, spec
        self.cache = ensure_dir(Path(repout)/'cache/lipid_contacts')
        self.out = ensure_dir(Path(repout)/'extensions/lipid_contacts')
        self.species, residues, self.tracker, self.initial = anchor_groups(universe, spec)
        self.axis = {'x':0, 'y':1, 'z':2}[spec.get('axis', 'z')]
        self.cutoff = float(spec.get('cutoff_A', 5.))
        self.batch = int(spec.get('batch_atoms', 4096))
        if not np.isfinite(self.cutoff) or self.cutoff <= 0 or self.batch < 1:
            raise ValueError('lipid_contacts requires positive cutoff_A and batch_atoms')
        self.atoms = select_heavy_atoms(universe.select_atoms(spec['selection']), policy='mechanistic')
        if not len(self.atoms):
            raise ValueError('lipid_contacts selection is empty')
        by_residue = dict(zip(map(int, residues), range(len(residues))))
        atom_map = np.full(len(universe.atoms), -1, dtype=int)
        seen = set()
        try:
            fragments = self.atoms.fragments
        except Exception as exc:
            raise ValueError('lipid_contacts requires bonded molecule identities') from exc
        for fragment in fragments:
            anchors = {by_residue[int(r)] for r in fragment.resindices if int(r) in by_residue}
            if len(anchors) != 1:
                raise ValueError(f'Every selected lipid molecule needs exactly one anchor residue; found {len(anchors)} in fragment with residues {sorted(set(fragment.resnames))}')
            anchor = anchors.pop()
            if anchor in seen:
                raise ValueError('One anchor cannot identify multiple lipid molecules')
            seen.add(anchor)
            atom_map[fragment.indices] = anchor
        if len(seen) != len(residues):
            raise ValueError('Lipid contact selection does not cover every anchor molecule')
        self.owners = atom_map[self.atoms.indices]
        self.rows, self.series = [], []
        atomic_json(self.out/'selection_audit.json', dict(
            molecules=len(residues), heavy_atoms=len(self.atoms), cutoff_A=self.cutoff,
            batch_atoms=self.batch, normalization='contact share / finite-population possible-pair share',
            interpretation='geometric proximity; not hydrogen bonds or affinity'))

    def update(self, ts, frame_index, time_ps):
        box = np.asarray(self.u.dimensions, float)
        if box.shape != (6,) or not np.isfinite(box).all() or not np.allclose(box[3:], 90.):
            raise ValueError('lipid_contacts leaflet labels require an orthorhombic planar membrane')
        leaf, changes = aligned_leaflets(self.tracker.positions()[:, self.axis], box[self.axis], self.initial)
        n = len(self.species)
        encoded = set()
        xyz = self.atoms.positions
        for first in range(0, len(xyz), self.batch):
            pairs = capped_distance(xyz[first:first+self.batch], xyz,
                                    self.cutoff, box=box, return_distances=False)
            a, b = self.owners[pairs[:, 0]+first], self.owners[pairs[:, 1]]
            keep = a != b
            a, b = np.minimum(a[keep], b[keep]), np.maximum(a[keep], b[keep])
            encoded.update(np.unique(a*n+b).tolist())
        codes = np.asarray(sorted(encoded), dtype=int)
        pairs = np.column_stack((codes//n, codes % n))
        table = contact_matrix(self.species, leaf, pairs, self.spec['species_selections'])
        table.insert(0, 'time_ps', time_ps)
        self.rows.extend(table.to_dict('records'))
        self.series.append(dict(time_ps=time_ps, lipid_contact_pair_count=len(pairs),
            lipid_contact_cross_leaflet_pair_count=int(np.count_nonzero(leaf[pairs[:, 0]] != leaf[pairs[:, 1]])),
            lipid_contact_leaflet_changes=changes))

    def flush(self, part):
        for key, rows in [('matrix', self.rows), ('timeseries', self.series)]:
            if rows:
                atomic_parquet(self.cache/f'{key}-{part:06d}.parquet', pd.DataFrame(rows))
                rows.clear()

    def prune_from(self, part):
        for path in self.cache.glob('*.parquet'):
            match = re.match(r'(matrix|timeseries)-(\d+)\.parquet$', path.name)
            if match and int(match[2]) >= part:
                path.unlink()

    def finalize(self):
        for key in ('matrix', 'timeseries'):
            paths = sorted(self.cache.glob(f'{key}-*.parquet'))
            table = pd.concat([pd.read_parquet(p) for p in paths], ignore_index=True) if paths else pd.DataFrame()
            atomic_parquet(self.out/f'{key}.parquet', table)
        return {'fused':True}
