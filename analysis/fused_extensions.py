"""Frame observers for expensive analyses, without additional XTC traversals.

All observers operate on the *raw* current Universe frame, before protein
alignment. MDAnalysis AnalysisBase algorithms are stepped with their native
_single_frame implementation, retaining the same scientific detectors as
standalone .run(). Raw event chunks are checkpointed atomically.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re
import time

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

from .config import analysis_groups, atomic_parquet, ensure_dir, time_ps
from .membrane_order import _atomic_npz
from .progress import track
from .hydrogen_bonds import write_hydrogen_bonds
from .lipid_hydrogen_bonds import (
    GLYCOLIPID_CHAIN_IDS, lipid_molecule_map, write_lipid_hydrogen_bonds,
)
from .sasa import _freesasa_structure
from .profile_common import _residue_label
from .biophysical_common import heavy, box_of, residue_pairs, sequence_motifs
from MDAnalysis.lib.distances import distance_array


def enabled_names(cfg):
    specs = cfg.get("analyses", {})
    names = []
    for name in ("hydrogen_bonds", "lipid_hydrogen_bonds", "dssp", "sasa", "local_environment"):
        if not specs.get(name, {}).get("enabled", False):
            continue
        # MDAnalysis parallel DSSP runs have different worker semantics.
        # Keep their existing specialized implementation.
        if name == "dssp" and specs[name].get("backend", "serial") != "serial":
            continue
        names.append(name)
    return tuple(names)


def _part_paths(path, extension):
    return sorted(path.glob("part-*" + extension))


def _concat(path):
    parts = _part_paths(path, ".parquet")
    return pd.concat((pd.read_parquet(p) for p in parts), ignore_index=True) if parts else pd.DataFrame()


class FusedExtensions:
    def __init__(self, universe, cfg, repout):
        self.u, self.cfg, self.repout = universe, cfg, Path(repout)
        self.names = enabled_names(cfg)
        self.cache = ensure_dir(self.repout / "cache" / "fused_extensions")
        self.rows = defaultdict(list)
        self.hb = {}
        self.dssp = None
        self.sasa = None
        self.local = None
        self.frames = 0
        self.timings = defaultdict(float)
        if "hydrogen_bonds" in self.names:
            self.hb["hydrogen_bonds"] = self._setup_protein_hba()
        if "lipid_hydrogen_bonds" in self.names:
            self.hb["lipid_hydrogen_bonds"] = self._setup_lipid_hba()
        if "dssp" in self.names:
            self._setup_dssp()
        if "sasa" in self.names:
            self._setup_sasa()
        if "local_environment" in self.names:
            self._setup_local()

    def _setup_protein_hba(self):
        spec = self.cfg["analyses"]["hydrogen_bonds"]
        between = spec.get("between")
        if between is None and spec.get("mode", "protein_protein") == "protein_ligand":
            ligand = self.cfg["selections"].get("ligand")
            if not ligand:
                raise ValueError("hydrogen_bonds.mode=protein_ligand requires selections.ligand")
            between = [self.cfg["selections"]["protein"], ligand]
        elif between is None:
            between = [self.cfg["selections"]["protein"]] * 2
        hba = HydrogenBondAnalysis(
            self.u, donors_sel=spec.get("donors_sel"),
            hydrogens_sel=spec.get("hydrogens_sel"),
            acceptors_sel=spec.get("acceptors_sel"),
            between=between,
            d_h_cutoff=float(spec.get("d_h_cutoff_A", 1.2)),
            d_a_cutoff=float(spec.get("d_a_cutoff_A", 3.0)),
            d_h_a_angle_cutoff=float(spec.get("angle_cutoff_deg", 150.)),
            update_selections=bool(spec.get("update_selections", False)),
        )
        hba._prepare()
        return {"hba": hba, "frames": [], "times": []}

    def _setup_lipid_hba(self):
        spec = self.cfg["analyses"]["lipid_hydrogen_bonds"]
        selection = str(spec.get("lipid_selection") or (
            "(resname POPC POPE POPI POPS DOPC DPPC PSM CHL1) "
            "or (chainID GM1a BGAL CER160)"
        ))
        atoms = self.u.select_atoms(selection)
        if not len(atoms):
            raise ValueError("lipid_hydrogen_bonds.lipid_selection matched no atoms")
        molecule_ids, atom_species = lipid_molecule_map(self.u, atoms)
        hydrogens = spec.get("hydrogens_sel")
        acceptors = spec.get("acceptors_sel")
        if hydrogens is None or acceptors is None:
            guesser = HydrogenBondAnalysis(
                self.u, donors_sel=None, hydrogens_sel="name __MDSCOPE_NONE__",
                acceptors_sel="name __MDSCOPE_NONE__", update_selections=False
            )
            if hydrogens is None:
                hydrogens = guesser.guess_hydrogens(selection)
            if acceptors is None:
                acceptors = guesser.guess_acceptors(selection)
        if not hydrogens or not acceptors:
            raise ValueError("No lipid donor hydrogens/acceptors inferred; supply explicit selectors")
        hydrogens = f"({hydrogens}) and ({selection})"
        acceptors = f"({acceptors}) and ({selection})"
        donors = spec.get("donors_sel")
        if donors is not None:
            donors = f"({donors}) and ({selection})"
        hba = HydrogenBondAnalysis(
            self.u, donors_sel=donors, hydrogens_sel=hydrogens,
            acceptors_sel=acceptors,
            d_h_cutoff=float(spec.get("d_h_cutoff_A", 1.2)),
            d_a_cutoff=float(spec.get("d_a_cutoff_A", 3.0)),
            d_h_a_angle_cutoff=float(spec.get("angle_cutoff_deg", 150.)),
            update_selections=bool(spec.get("update_selections", False)),
        )
        hba._prepare()
        return {"hba": hba, "frames": [], "times": [],
                "atoms": atoms, "molecule_ids": molecule_ids,
                "atom_species": atom_species}

    def _setup_dssp(self):
        from MDAnalysis.analysis.dssp import DSSP
        spec = self.cfg["analyses"]["dssp"]
        selection = str(spec.get("selection") or self.cfg["selections"]["protein"])
        atoms = self.u.select_atoms(selection)
        if not len(atoms):
            raise ValueError(f"dssp.selection matched zero atoms: {selection}")
        dssp = DSSP(
            atoms, guess_hydrogens=bool(spec.get("guess_hydrogens", True)),
            heavyatom_names=tuple(spec.get("heavyatom_names", ("N", "CA", "C", "O O1 OT1"))),
            hydrogen_name=str(spec.get("hydrogen_name", "H HN HT1 HT2 HT3")),
        )
        dssp._prepare()
        residues = list(atoms.select_atoms("protein").residues)
        resids = np.asarray(dssp._heavy_atoms["CA"].resids)
        if len(residues) != len(resids):
            mapped, used = [], set()
            for resid in resids:
                candidates = [r for r in residues if r.resid == resid and r.ix not in used]
                if not candidates:
                    raise RuntimeError("Could not map DSSP residues back to topology residues")
                mapped.append(candidates[0])
                used.add(candidates[0].ix)
            residues = mapped
        self.dssp = {"analysis": dssp, "residues": residues, "times": []}

    def _setup_sasa(self):
        spec = self.cfg["analyses"]["sasa"]
        selection = str(spec.get("selection") or self.cfg["selections"]["protein"])
        atoms = self.u.select_atoms(selection)
        if not len(atoms):
            raise ValueError(f"sasa.selection matched zero atoms: {selection}")
        residues = list(atoms.residues)
        lookup = {int(res.ix): i for i, res in enumerate(residues)}
        mapping = np.fromiter(
            (lookup[int(ix)] for ix in atoms.resindices), dtype=int, count=len(atoms)
        )
        import freesasa
        params = freesasa.Parameters()
        algorithm = str(spec.get("algorithm", "shrake_rupley")).lower()
        if algorithm in {"shrake_rupley", "shrake-rupley", "sr"}:
            params.setAlgorithm(freesasa.ShrakeRupley)
            if spec.get("n_points") is not None:
                params.setNPoints(int(spec["n_points"]))
        elif algorithm in {"lee_richards", "lee-richards", "lr"}:
            params.setAlgorithm(freesasa.LeeRichards)
            if spec.get("n_slices") is not None:
                params.setNSlices(int(spec["n_slices"]))
        else:
            raise ValueError("sasa.algorithm must be shrake_rupley or lee_richards")
        if spec.get("probe_radius_A") is not None:
            params.setProbeRadius(float(spec["probe_radius_A"]))
        if spec.get("n_threads") is not None:
            params.setNThreads(int(spec["n_threads"]))
        self.sasa = {"atoms": atoms, "residues": residues,
                     "mapping": mapping, "params": params}

    def _setup_local(self):
        s = self.cfg["analyses"]["local_environment"]
        prot = heavy(self.u.select_atoms(s.get("protein_selection") or self.cfg["selections"]["protein"]))
        lip = heavy(self.u.select_atoms(s.get("lipid_selection", "resname POPC POPE POPI POPS DOPC DPPC PSM CHL1 CER160 CER1 BGAL GM1a")))
        if not len(lip) or not len(prot):
            raise ValueError("No protein/lipid atoms matched local_environment selections")
        glycan = heavy(self.u.select_atoms(s.get("glycan_selection", "resname AGLCNA AIDOA BGLCA BGAL BXYL")))
        ions = heavy(self.u.select_atoms(s.get("ion_selection", "resname SOD NA CLA CL POT K CAL CA")))
        cutoff = float(s.get("contact_cutoff_A", 5.))
        motifs = sequence_motifs(prot.residues, s.get("motifs", {
            "nxip": "N.IP", "gxxxg": "G...G", "hs_candidate": "SG"
        }))
        groups = {}
        for hit in motifs:
            key = "motif_" + hit["motif"] + "_" + hit["chain"] + "_" + str(hit["start_resid"])
            groups[key] = heavy(self.u.atoms[np.isin(self.u.atoms.resindices, hit["resindices"])])
        grouped = {str(name): lip[lip.resnames == name] for name in np.unique(lip.resnames)}
        meta = [{"motif_id": key, "sequence": "",
                 "resids": ",".join(map(str, group.residues.resids))}
                for key, group in groups.items()]
        self.local = {"prot": prot, "lip": lip, "glycan": glycan, "ions": ions,
                      "groups": groups, "grouped": grouped, "meta": meta, "cutoff": cutoff}

    def update(self, ts, frame_index):
        t = time_ps(ts, frame_index, self.cfg)
        for name, entry in self.hb.items():
            start = time.perf_counter()
            hba = entry["hba"]
            hba._ts = ts
            hba._single_frame()
            entry["frames"].append(int(ts.frame))
            entry["times"].append(t)
            self.timings[name] += time.perf_counter() - start
        if self.dssp is not None:
            start = time.perf_counter()
            self.dssp["analysis"]._single_frame()
            self.dssp["times"].append(t)
            self.timings["dssp"] += time.perf_counter() - start
        if self.sasa is not None:
            start = time.perf_counter()
            import freesasa
            state = self.sasa
            _, structure = _freesasa_structure(state["atoms"])
            result = freesasa.calc(structure, state["params"])
            areas = np.fromiter(
                (result.atomArea(i) for i in range(len(state["atoms"]))),
                dtype=float, count=len(state["atoms"])
            )
            residue_areas = np.bincount(
                state["mapping"], weights=areas,
                minlength=len(state["residues"])
            ) / 100.
            self.rows["sasa"].extend(
                (t, int(res.ix), _residue_label(res), float(area))
                for res, area in zip(state["residues"], residue_areas)
            )
            self.timings["sasa"] += time.perf_counter() - start
        if self.local is not None:
            start = time.perf_counter()
            self._update_local(t)
            self.timings["local_environment"] += time.perf_counter() - start
        self.frames += 1

    def _update_local(self, t):
        s = self.local
        prot, lip, cutoff = s["prot"], s["lip"], s["cutoff"]
        box = box_of(self.u)
        row = {"time_ps": t}
        cp = residue_pairs(prot, lip, cutoff, box)
        contacted = np.unique(cp[:, 0]) if len(cp) else np.empty(0, int)
        row["local_protein_lipid_contact_residues_count"] = len(contacted)
        if len(cp):
            a, b = np.unique(cp[:, 0], return_counts=True)
            self.rows["local_residue"].extend((t, int(resid), int(n)) for resid, n in zip(a, b))
        for species, atoms in s["grouped"].items():
            hits = residue_pairs(prot, atoms, cutoff, box)
            row["local_" + species + "_contact_residues_count"] = len(np.unique(hits[:, 0])) if len(hits) else 0
        for name, group in s["groups"].items():
            if len(group):
                d = distance_array(group.positions, lip.positions, box=box)
                row[name + "_min_lipid_distance_A"] = float(d.min())
                row[name + "_lipid_contact_fraction"] = float(np.mean(d.min(axis=1) <= cutoff))
                others = prot[~np.isin(prot.indices, group.indices)]
                if len(others):
                    d2 = distance_array(group.positions, others.positions, box=box)
                    row[name + "_protein_occlusion_fraction"] = float(np.mean(np.any(d2 < cutoff, axis=1)))
        if len(s["glycan"]):
            contacts = residue_pairs(s["glycan"], lip, cutoff, box)
            row["local_hs_lipid_contact_pairs_count"] = len(contacts)
            if len(s["ions"]):
                row["local_hs_ion_contacts_count"] = len(residue_pairs(s["glycan"], s["ions"], 3.5, box))
        self.rows["local_environment"].append(row)

    def flush(self, part):
        start = time.perf_counter()
        for name, entry in self.hb.items():
            hba = entry["hba"]
            hba._conclude()
            raw = np.asarray(hba.results.hbonds, dtype=float).reshape(-1, 6)
            _atomic_npz(
                ensure_dir(self.cache / name) / f"part-{part:06d}.npz",
                events=raw, frames=np.asarray(entry["frames"], dtype=int),
                times=np.asarray(entry["times"], dtype=float),
            )
            hba._prepare()
            entry["frames"].clear()
            entry["times"].clear()
        if self.dssp is not None:
            entry = self.dssp
            analysis = entry["analysis"]
            analysis._conclude()
            codes = np.asarray(analysis.results.dssp).astype(str)
            codes[codes == "-"] = "C"
            nframes, nres = codes.shape
            residues = entry["residues"]
            atomic_parquet(
                ensure_dir(self.cache / "dssp") / f"part-{part:06d}.parquet",
                pd.DataFrame({
                    "time_ps": np.repeat(entry["times"], nres),
                    "topology_residue_index": np.tile([int(r.ix) for r in residues], nframes),
                    "residue": np.tile([_residue_label(r) for r in residues], nframes),
                    "dssp": codes.ravel(),
                })
            )
            analysis._prepare()
            entry["times"].clear()
        for name, rows in self.rows.items():
            if rows:
                if name == "sasa":
                    df = pd.DataFrame(rows, columns=[
                        "time_ps", "topology_residue_index", "residue", "sasa_nm2"
                    ])
                elif name == "local_residue":
                    df = pd.DataFrame(rows, columns=[
                        "time_ps", "resindex", "n_lipid_residues"
                    ])
                else:
                    df = pd.DataFrame(rows)
                atomic_parquet(
                    ensure_dir(self.cache / name) / f"part-{part:06d}.parquet", df
                )
                rows.clear()
        self.timings["chunk_write"] += time.perf_counter() - start

    def prune_from(self, next_part):
        for path in self.cache.rglob("part-*"):
            match = re.search(r"part-(\d+)", path.name)
            if match and int(match.group(1)) >= next_part and path.is_file():
                path.unlink()

    def finalize(self):
        start = time.perf_counter()
        results = {}
        for name, entry in self.hb.items():
            parts = _part_paths(self.cache / name, ".npz")
            events, frames, times = [], [], []
            for path in parts:
                with np.load(path, allow_pickle=False) as data:
                    events.append(data["events"])
                    frames.append(data["frames"])
                    times.append(data["times"])
            arr = np.concatenate(events) if events else np.empty((0, 6))
            frame_arr = np.concatenate(frames) if frames else np.array([], int)
            time_arr = np.concatenate(times) if times else np.array([], float)
            if name == "hydrogen_bonds":
                results[name] = write_hydrogen_bonds(
                    self.u, self.cfg, self.repout, arr, frame_arr, time_arr
                )
            else:
                results[name] = write_lipid_hydrogen_bonds(
                    self.u, self.cfg, self.repout, arr, frame_arr, time_arr,
                    entry["atoms"], entry["molecule_ids"], entry["atom_species"]
                )
        for name, output in (
            ("dssp", self.repout / "profiles" / "dssp" / "assignments.parquet"),
            ("sasa", self.repout / "profiles" / "sasa" / "residue_sasa.parquet"),
            ("local_environment", self.repout / "extensions" / "local_environment" / "timeseries.parquet"),
        ):
            if name in self.names:
                atomic_parquet(output, _concat(self.cache / name))
                results[name] = {"frames": self.frames, "fused": True}
        if self.local is not None:
            if _part_paths(self.cache / "local_residue", ".parquet"):
                atomic_parquet(
                    self.repout / "extensions" / "local_environment" / "residue_contacts.parquet",
                    _concat(self.cache / "local_residue")
                )
            if self.local["meta"]:
                atomic_parquet(
                    self.repout / "extensions" / "local_environment" / "motifs.parquet",
                    pd.DataFrame(self.local["meta"])
                )
        self.timings["finalize"] += time.perf_counter() - start
        from .config import atomic_json
        atomic_json(self.repout / "fused_extensions_timings.json", {
            "enabled": list(self.names), "frames_processed_this_invocation": self.frames,
            "compute_seconds_by_analysis": dict(self.timings),
        })
        return results


def run_fused_extensions(rep, cfg, repout):
    """One shared fallback trajectory pass if protein streaming is disabled."""
    if not enabled_names(cfg):
        return {}
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    tr = cfg["trajectory"]
    start = 0 if tr.get("start") is None else int(tr["start"])
    u.trajectory[start]
    observer = FusedExtensions(u, cfg, repout)
    observer.prune_from(0)
    frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
    chunk = max(1, int(tr.get("chunk_frames", 1000)))
    for i, ts in enumerate(track(frames, total=len(frames),
                                 desc=f"{rep.system}/{rep.label}: fused extensions",
                                 enabled=cfg.get("progress", {}).get("enabled", True))):
        observer.update(ts, i)
        if (i + 1) % chunk == 0:
            observer.flush(i // chunk)
    if observer.frames % chunk:
        observer.flush(observer.frames // chunk)
    return observer.finalize()
