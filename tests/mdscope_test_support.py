"""Shared tiny MD systems and CLI project builders used by unrelated test modules.

Keep only fixtures with multiple consumers. These builders must not change
scientific expectations or introduce a production dependency on test data.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import yaml
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.pipeline.fused_membrane import FusedMembrane


def membrane_universe():
    u = mda.Universe.empty(
        10, n_residues=10, atom_resindex=np.arange(10), trajectory=True
    )
    u.add_TopologyAttr("names", ["P"] * 8 + ["CA", "CA"])
    u.add_TopologyAttr("resnames", ["POPC"] * 8 + ["ALA"] * 2)
    u.add_TopologyAttr("resids", list(range(1, 11)))
    x = np.array([10., 25., 10., 25., 10., 25., 10., 25.])
    y = np.array([10., 10., 25., 25., 10., 10., 25., 25.])
    z = np.array([20., 20., 20., 20., 60., 60., 60., 60.])
    first = np.column_stack([x, y, z])
    protein = np.array([[99., 15., 40.], [1., 15., 40.]])
    coordinates = []
    for i in range(5):
        frame = np.vstack([first + [i * .1, 0, 0], protein])
        frame[8, 0] -= 100. if i % 2 else 0.
        coordinates.append(frame)
    boxes = np.tile([100., 100., 100., 90., 90., 90.], (5, 1))
    u.load_new(np.asarray(coordinates, dtype=np.float32),
               format=MemoryReader, dimensions=boxes)
    return u


def membrane_config():
    selections = {"POPC": "resname POPC and name P"}
    return {
        "trajectory": {"start": 0, "stop": 5, "stride": 1, "chunk_frames": 2},
        "time": {"dt_ps": None, "origin_ps": 0.},
        "selections": {"protein": "name CA"},
        "analyses": {
            "membrane_order": {"enabled": False},
            "bilayer_normal": {
                "enabled": True, "axis": "z",
                "species_selections": selections,
                "neighbor_radius_A": 30., "min_neighbors": 3,
                "save_lipid_angles": True,
            },
            "lipid_enrichment": {
                "enabled": True, "axis": "z",
                "species_selections": selections,
                "footprint_selection": "name CA",
                "slab_half_width_A": 25., "shell_edges_A": [0, 15, 30],
            },
        },
    }


def execute_membrane(u, cfg, root, resume=False):
    u.trajectory[0]
    observer = FusedMembrane(u, cfg, root)
    if resume:
        for k in range(3):
            observer.update(u.trajectory[k], k)
            if k == 1:
                observer.flush(0)
        observer.flush(1)  # simulate uncommitted part written before crash
        del observer
        u.trajectory[0]
        observer = FusedMembrane(u, cfg, root)
        observer.prune_from(1)
        start = 2
        next_part = 1
    else:
        observer.prune_from(0)
        start = 0
        next_part = 0
    for k in range(start, 5):
        observer.update(u.trajectory[k], k)
        if (k + 1) % 2 == 0:
            observer.flush(next_part)
            next_part += 1
    if observer.frame_count % 2:
        observer.flush(next_part)
    return observer.finalize()


def solvation_universe():
    u = mda.Universe.empty(5, n_residues=4, atom_resindex=[0, 0, 1, 2, 3],
                           trajectory=True)
    u.add_TopologyAttr("names", ["OE1", "OE2", "OW", "OW", "SOD"])
    u.add_TopologyAttr("resnames", ["GLU", "SOL", "SOL", "SOD"])
    u.add_TopologyAttr("resids", [7, 1, 2, 3])
    xyz = np.array([[1, 1, 1], [1.5, 1, 1], [3, 1, 1],
                    [19.5, 1, 1], [3, 1, 2]], np.float32)
    shifted = xyz.copy()
    shifted[2] = [12, 1, 1]
    u.load_new(np.array([xyz, shifted]), format=MemoryReader,
               dimensions=np.tile([20, 20, 20, 90, 90, 90], (2, 1)))
    return u


def invoke_cli(tmp_path: Path, *args: str, stdin: str | None = None, ok: bool = True):
    command = ["mdscope", *map(str, args)]
    result = subprocess.run(
        command, input=stdin, text=True, capture_output=True,
        cwd=tmp_path, timeout=150, check=False,
    )
    if ok:
        assert result.returncode == 0, (
            f"{command!r} failed ({result.returncode})\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    else:
        assert result.returncode != 0, f"{command!r} unexpectedly succeeded"
    return result


def tiny_cli_project(tmp_path):
    pdb = tmp_path / "tiny.pdb"
    xyz = np.array([
        [10., 10., 10.], [13., 10., 10.],
        [10., 13., 10.], [10., 10., 13.],
    ])
    lines = ["CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1"]
    for i, (x, y, z) in enumerate(xyz, 1):
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}           C"
        )
    pdb.write_text("\n".join([*lines, "TER", "END", ""]))
    universe = mda.Universe(str(pdb))
    trajectory = tmp_path / "tiny.dcd"
    with mda.Writer(str(trajectory), n_atoms=len(universe.atoms)) as writer:
        for frame in range(6):
            positions = xyz.copy()
            positions[-1, 2] += frame * 0.12
            universe.atoms.positions = positions
            universe.dimensions = np.array([50., 50., 50., 90., 90., 90.])
            writer.write(universe.atoms)

    def config(label="results", statistics=False):
        analyses = {
            name: {"enabled": False}
            for name in (
                "membrane", "membrane_order", "dssp", "sasa",
                "salt_bridges", "hydrogen_bonds", "ligand_contacts",
                "density_1d", "residue_membrane_distance",
                "bilayer_normal", "lipid_enrichment",
                "lipid_hydrogen_bonds", "sdc4_mechanisms",
                "lifetimes",
            )
        }
        analyses.update({
            "alignment": {"enabled": True},
            "rmsd": {"enabled": True},
            "rmsf": {"enabled": True},
            "radius_of_gyration": {"enabled": True},
        })
        specification = {
            "project": {"name": "cli_smoke", "output": label},
            "trajectory": {"stride": 1, "chunk_frames": 3},
            "selections": {
                "protein": "protein", "alignment": "protein and name CA",
                "features": "protein and name CA", "ligand": None,
            },
            "systems": {
                "tiny": {
                    "topology": pdb.name,
                    "replicates": [{"trajectory": trajectory.name, "label": "rep1"}],
                },
            },
            "analyses": analyses,
            "statistics": {"enabled": statistics, "bootstrap_draws": 100},
            "plots": {"enabled": False},
        }
        path = tmp_path / f"{label}.yaml"
        path.write_text(yaml.safe_dump(specification, sort_keys=False))
        return path, tmp_path / label
    return config
