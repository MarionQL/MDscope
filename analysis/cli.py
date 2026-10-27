"""
Tiny CLI for creating, validating, running, selecting one replicate, and
finalizing MDscope.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import MDAnalysis as mda
import yaml

from .config import analysis_runs, iter_replicates, load_config
from .profiles import _without_hydrogens
from .runner import finalize_project, run_all, run_one, project_root
from .report import generate_report
from .timeseries import resolved_selection_counts


def interactive_init(output="analysis.yaml"):
    print("MDscope configuration generator (analysis runs themselves are non-interactive).")
    name = input("Project name [mdscope_project]: ").strip() or "mdscope_project"
    topology = input("Topology path: ").strip()
    trajectories = [
        x.strip() for x in input("Trajectory paths, comma-separated: ").split(",") if x.strip()
    ]
    ligand = input("Ligand MDAnalysis selection (blank for none): ").strip() or None
    reps = [{"trajectory": path, "label": f"rep{i}"} for i, path in enumerate(trajectories, 1)]
    cfg = {
        "project": {"name": name, "output": "results"},
        "trajectory": {"stride": 10, "chunk_frames": 1000},
        "selections":{
            "protein": "protein",
            "alignment": "protein and name CA",
            "features": "protein and name CA",
            "ligand": ligand
        },
        "systems": {"system1": {"topology": topology, "replicates": reps}},
        "analyses":{
            "rmsd": True,
            "rmsf": True,
            "radius_of_gyration": True,
            "dssp": True,
            "salt_bridges": True,
            "sasa": {"enabled": True}
        },
    }
    if ligand:
        cfg["analyses"]["ligand_contacts"] = {"enabled": True, "cutoff_A": 4.0}
    Path(output).write_text(yaml.safe_dump(cfg, sort_keys=False))
    return Path(output)


def cmd_validate(args):
    cfg = load_config(args.config)
    reps = list(iter_replicates(cfg, args.config))
    print(f"Valid config: {len(cfg['systems'])} systems, {len(reps)} replicates")
    for rep in reps:
        print(f"  {rep.system}/{rep.label}: {rep.trajectory}")
    print("\nResolved structural selections:")
    base = Path(args.config).resolve().parent
    for system, spec in cfg["systems"].items():
        topology = Path(spec["topology"])
        topology = topology if topology.is_absolute() else (base / topology).resolve()
        universe = mda.Universe(str(topology))
        counts = resolved_selection_counts(universe, cfg)
        print(
            f"  {system}: protein={counts['protein']}, "
            f"alignment={counts['alignment']}, features={counts['features']}"
        )
        if counts["protein"] == 0 or counts["features"] == 0:
            raise ValueError(
                f"System {system!r} has a zero-atom protein/feature selection"
            )
        if (
            cfg["analyses"]["alignment"].get("enabled", True)
            and counts["alignment"] == 0
        ):
            raise ValueError(
                f"System {system!r} has a zero-atom alignment selection"
            )
        for group_id, group in counts["groups"].items():
            print(
                f"    group {group_id}: group={group['group_atoms']}, "
                f"features={group['feature_atoms']}, "
                f"alignment={group['alignment_atoms']}"
            )
            if not group["group_atoms"] or not group["feature_atoms"]:
                raise ValueError(
                    f"System {system!r}, group {group_id!r} resolves to zero atoms"
                )
            if (
                cfg["analyses"]["alignment"].get("enabled", True)
                and not group["alignment_atoms"]
            ):
                raise ValueError(
                    f"System {system!r}, group {group_id!r} has a zero-atom "
                    "alignment selection"
                )
        for kind in ("bilayer_normal","lipid_enrichment"):
            extension=cfg.get("analyses",{}).get(kind,{})
            if not extension.get("enabled",False):
                continue
            lipid_residues=set()
            matched=0
            for species,selection in extension.get("species_selections",{}).items():
                atoms=universe.select_atoms(str(selection))
                if not len(atoms):
                    if extension.get("allow_absent_species",False):
                        print(f"    {kind}: {species} absent (permitted)")
                        continue
                    raise ValueError(f"{system}: {kind} species {species} matched no atoms")
                ids=set(map(int,atoms.resindices))
                if lipid_residues.intersection(ids):
                    raise ValueError(f"{system}: {kind} lipid anchor residue selections overlap")
                lipid_residues.update(ids)
                matched+=len(ids)
                print(f"    {kind}: {species} {len(ids)} lipid anchors")
            if matched<4:
                raise ValueError(f"{system}: {kind} requires >=4 lipid anchors")
            if kind=="lipid_enrichment":
                protein=universe.select_atoms(str(
                    extension.get("protein_selection") or cfg["selections"]["protein"]))
                if not len(protein):
                    raise ValueError(f"{system}: lipid enrichment protein selection matched no atoms")
        hb=cfg.get("analyses",{}).get("lipid_hydrogen_bonds",{})
        if hb.get("enabled",False):
            lipid=universe.select_atoms(str(hb["lipid_selection"]))
            if not len(lipid):
                raise ValueError(f"{system}: lipid H-bond selection matched no atoms")
            print(f"    lipid_hydrogen_bonds: {len(lipid)} lipid atoms")
        for run in analysis_runs(cfg, "residue_membrane_distance"):
            membrane_selection = run.spec.get("membrane_selection") or "resname " + " ".join(map(
                str,
                run.spec.get("lipid_resnames", [])
            ))
            membrane = _without_hydrogens(
                universe.select_atoms(str(membrane_selection)),
                bool(run.spec.get("exclude_hydrogens", True))
            )
            if not len(membrane):
                raise ValueError(
                    f"System {system!r}, residue-membrane run {run.run_id!r} "
                    "matched zero membrane atoms"
                )
            print(
                f"    residue-membrane {run.run_id}: "
                f"{len(membrane)} membrane atoms"
            )



def main():
    parser = argparse.ArgumentParser(prog="mdscope")
    sub = parser.add_subparsers(required=True)
    p = sub.add_parser("init")
    p.add_argument("--output", default="analysis.yaml")
    p.set_defaults(func=lambda a: print(interactive_init(a.output)))
    p = sub.add_parser("validate")
    p.add_argument("config")
    p.set_defaults(func=cmd_validate)
    p = sub.add_parser("run")
    p.add_argument("config")
    p.set_defaults(func=lambda a: print(run_all(a.config)))
    p = sub.add_parser("replicate", help="run one zero-based systemxreplicate entry")
    p.add_argument("config")
    p.add_argument("--index", type=int, required=True)
    p.set_defaults(func=lambda a: print(run_one(a.config, a.index)))
    p = sub.add_parser("finalize")
    p.add_argument("config")
    p.set_defaults(func=lambda a: print(finalize_project(a.config)))
    p = sub.add_parser("report", help="Rebuild PDF/ZIP from existing finalized tables without re-running or archiving")
    p.add_argument("config")
    def _report(args):
        cfg = load_config(args.config)
        root = project_root(cfg, args.config)
        reps = list(iter_replicates(cfg, args.config))
        outputs = {
            system: [root / rep.system / rep.label for rep in reps if rep.system == system]
            for system in cfg["systems"]
        }
        print(generate_report(cfg, root, outputs, force=True))
    p.set_defaults(func=_report)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
