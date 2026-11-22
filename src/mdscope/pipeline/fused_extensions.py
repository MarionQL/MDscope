"""Frame observers for expensive analyses, without additional XTC traversals.

All observers operate on the *raw* current Universe frame, before protein
alignment. MDAnalysis AnalysisBase algorithms are stepped with their native
_single_frame implementation, retaining the same scientific detectors as
standalone .run(). Raw event chunks are checkpointed atomically.
"""
from __future__ import annotations

import re
import time
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

from ..config import time_ps
from ..core.persistence import atomic_parquet, ensure_dir
from ..core.profiling import _residue_label
from ..core.progress import track
from ..interactions.hydrogen_bonds import write_hydrogen_bonds
from ..interactions.lipid_hydrogen_bonds import (
    lipid_molecule_map,
    write_lipid_hydrogen_bonds,
)
from ..observables.hydration import HydrationObserver
from ..mechanisms.local_environment import (
    local_environment_frame,
    prepare_local_environment,
)
from ..observables.membrane_order import _atomic_npz
from ..observables.sasa import sasa_instances, sasa_residue_rows, write_sasa_burial
from ..observables.solvation import SolvationObserver
from .analysis_selection import extensions_enabled_names as enabled_names


def _part_paths(path, extension):
    return sorted(path.glob("part-*" + extension))


def _concat(path):
    parts = _part_paths(path, ".parquet")
    return pd.concat((pd.read_parquet(p) for p in parts), ignore_index=True) if parts else pd.DataFrame()


class FusedExtensions:
    def __init__(self, universe, cfg, repout, *, names=None, sasa_run_ids=None, cache_dir=None):
        self.u, self.cfg, self.repout = universe, cfg, Path(repout)
        self.names = enabled_names(cfg) if names is None else tuple(names)
        self.sasa_run_ids = None if sasa_run_ids is None else tuple(sasa_run_ids)
        self.cache = ensure_dir(Path(cache_dir) if cache_dir is not None
                                else self.repout / "cache" / "fused_extensions")
        self.rows = defaultdict(list)
        self.hb = {}
        self.dssp = None
        self.sasa = None
        self.local = None
        self.solvation = None
        self.hydration = None
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
        if "solvation" in self.names:
            opts = ({"cache_dir": self.cache / "solvation"}
                    if cache_dir is not None else {})
            self.solvation = SolvationObserver(universe, cfg, repout, **opts)
        if "hydration" in self.names:
            opts = ({"cache_dir": self.cache / "hydration"}
                    if cache_dir is not None else {})
            self.hydration = HydrationObserver(universe, cfg, repout, **opts)

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
        donors_sel = spec.get("donors_sel")
        hydrogens_sel = spec.get("hydrogens_sel")
        acceptors_sel = spec.get("acceptors_sel")
        # Opt-in fast path. MDAnalysis applies 'between' *after* generating
        # donor/acceptor candidates; on solvated systems this can include
        # hundreds of thousands of water atoms. Restrict candidate chemistry
        # before the distance search, retaining the original between filter.
        # Keep explicit selections and other modes unchanged.
        if (spec.get("restrict_search_to_protein", False)
                and spec.get("mode", "protein_protein") == "protein_protein"
                and spec.get("between") is None):
            protein = self.cfg["selections"]["protein"]
            if hydrogens_sel is None or acceptors_sel is None:
                guesser = HydrogenBondAnalysis(
                    self.u, donors_sel=None,
                    hydrogens_sel="name __MDSCOPE_NONE__",
                    acceptors_sel="name __MDSCOPE_NONE__",
                    update_selections=False,
                )
                if hydrogens_sel is None:
                    hydrogens_sel = guesser.guess_hydrogens(protein)
                if acceptors_sel is None:
                    acceptors_sel = guesser.guess_acceptors(protein)
            if not hydrogens_sel or not acceptors_sel:
                raise ValueError("No protein H-bond hydrogens/acceptors inferred")
            hydrogens_sel = f"({hydrogens_sel}) and ({protein})"
            acceptors_sel = f"({acceptors_sel}) and ({protein})"
            if donors_sel is not None:
                donors_sel = f"({donors_sel}) and ({protein})"
        hba = HydrogenBondAnalysis(
            self.u, donors_sel=donors_sel,
            hydrogens_sel=hydrogens_sel,
            acceptors_sel=acceptors_sel,
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
        # Validate every named run/reference against the original config,
        # but calculate and persist only the requested runs.
        states = sasa_instances(self.u, self.cfg)
        if self.sasa_run_ids is None:
            self.sasa = states
            return
        unknown = set(self.sasa_run_ids) - set(states)
        if unknown:
            raise ValueError(f"Unknown SASA run IDs: {sorted(unknown)}")
        self.sasa = {key: value for key, value in states.items()
                     if key in self.sasa_run_ids}

    def _setup_local(self):
        self.local = prepare_local_environment(self.u, self.cfg)

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
            for run_id, state in self.sasa.items():
                name = "sasa" if run_id == "default" else f"sasa__{run_id}"
                atom_sink = (self.rows[f"sasa_atoms_{run_id}"]
                             if len(state["atom_positions"]) else None)
                self.rows[name].extend(sasa_residue_rows(state, t, atom_sink))
            self.timings["sasa"] += time.perf_counter() - start
        if self.local is not None:
            start = time.perf_counter()
            self._update_local(t)
            self.timings["local_environment"] += time.perf_counter() - start
        if self.solvation is not None:
            start = time.perf_counter()
            self.solvation.update(ts, frame_index, t)
            self.timings['solvation'] += time.perf_counter() - start
        if self.hydration is not None:
            start = time.perf_counter()
            self.hydration.update(ts, frame_index, t)
            self.timings['hydration'] += time.perf_counter() - start
        self.frames += 1

    def _update_local(self, t):
        row, residue_contacts = local_environment_frame(self.local, self.u, t)
        self.rows["local_environment"].append(row)
        self.rows["local_residue"].extend(residue_contacts)

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
                if name.startswith("sasa_atoms_"):
                    df = pd.DataFrame(rows, columns=[
                        "time_ps", "atom_index", "topology_residue_index",
                        "atom_name", "sasa_nm2",
                    ])
                elif name == "sasa" or name.startswith("sasa__"):
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
        if self.solvation is not None:
            self.solvation.flush(part)
        if self.hydration is not None:
            self.hydration.flush(part)
        self.timings["chunk_write"] += time.perf_counter() - start

    def prune_from(self, next_part):
        if self.solvation is not None:
            self.solvation.prune_from(next_part)
        if self.hydration is not None:
            self.hydration.prune_from(next_part)
        active = set(self.hb)
        if self.dssp is not None:
            active.add("dssp")
        if self.local is not None:
            active.update(("local_environment", "local_residue"))
        if self.sasa is not None:
            for run_id in self.sasa:
                active.add("sasa" if run_id == "default" else f"sasa__{run_id}")
                active.add(f"sasa_atoms_{run_id}")
        for name in active:
            for path in (self.cache / name).glob("part-*"):
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
            if name in self.names and (name != "sasa" or "default" in self.sasa):
                atomic_parquet(output, _concat(self.cache / name))
                results[name] = {"frames": self.frames, "fused": True}
        if self.sasa is not None:
            for run_id in self.sasa:
                if run_id == "default":
                    continue
                name = f"sasa__{run_id}"
                output = self.repout / "profiles" / "sasa" / run_id / "residue_sasa.parquet"
                atomic_parquet(output, _concat(self.cache / name))
                results[name] = {"frames": self.frames, "fused": True}
            for run_id, state in self.sasa.items():
                if len(state["atom_positions"]):
                    output = (self.repout / "profiles" / "sasa" /
                              (run_id if run_id != "default" else "") /
                              "atom_sasa.parquet")
                    atomic_parquet(output, _concat(self.cache / f"sasa_atoms_{run_id}"))
                if state["reference_run"] is not None:
                    write_sasa_burial(self.repout, run_id, state["reference_run"])
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
        if self.solvation is not None:
            results['solvation'] = self.solvation.finalize()
        if self.hydration is not None:
            results['hydration'] = self.hydration.finalize()
        self.timings["finalize"] += time.perf_counter() - start
        from ..core.persistence import atomic_json
        atomic_json(self.repout / "fused_extensions_timings.json", {
            "enabled": list(self.names), "frames_processed_this_invocation": self.frames,
            "compute_seconds_by_analysis": dict(self.timings),
        })
        return results


def run_fused_extensions(rep, cfg, repout, *, names=None, sasa_run_ids=None):
    """Replay selected missing extensions without touching committed stream chunks."""
    if not enabled_names(cfg):
        return {}
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    tr = cfg["trajectory"]
    start = 0 if tr.get("start") is None else int(tr["start"])
    u.trajectory[start]
    chunk = max(1, int(tr.get("chunk_frames", 1000)))
    with TemporaryDirectory(dir=ensure_dir(Path(repout)), prefix=".extensions-recovery-") as scratch:
        observer = FusedExtensions(
            u, cfg, repout, names=names, sasa_run_ids=sasa_run_ids,
            cache_dir=scratch,
        )
        observer.prune_from(0)
        frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
        for i, ts in enumerate(track(frames, total=len(frames),
                                     desc=f"{rep.system}/{rep.label}: fused extensions",
                                     enabled=cfg.get("progress", {}).get("enabled", True))):
            observer.update(ts, i)
            if (i + 1) % chunk == 0:
                observer.flush(i // chunk)
        if observer.frames % chunk:
            observer.flush(observer.frames // chunk)
        return observer.finalize()
