"""
Tiny CLI for creating, validating, running, selecting one replicate, and
finalizing MDscope.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from .config import iter_replicates, load_config
from .core.validation import validate_selections
from .pipeline.runner import finalize_project, run_all, run_one
from .reporting.report import generate_report
from .statistics.residue_regions import scan_saved_residue_regions


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
    validate_selections(cfg, Path(args.config).resolve().parent, emit=print)


def cmd_plan(args):
    """Print a no-run output plan; topology inspection is opt-in."""
    import json
    from .planning import build_plan

    cfg = load_config(args.config, check_paths=args.check_inputs)
    counts = (validate_selections(cfg, Path(args.config).resolve().parent)
              if args.check_topology else None)
    mode = ("topology_aware" if args.check_topology else
            "input_paths" if args.check_inputs else "schema_only")
    plan = build_plan(cfg, args.config, validation_mode=mode,
                      topology_counts=counts)
    if args.format == "json":
        print(json.dumps(plan, indent=2, sort_keys=True))
        return
    print(f"MDscope plan ({mode}): {len(plan['systems'])} systems, "
          f"{len(plan['replicates'])} replicas")
    for item in plan["analyses"]:
        print(f"  {item['analysis']}[{item['run_id']}] "
              f"source={item['coordinate_source']} "
              f"artifacts={item['artifact_catalog_status']}")
    print("Registered per-replica files (required when their observer completes):")
    for item in plan["registered_replica_outputs"]:
        print(f"  {item['path_within_replica']} "
              f"({item['format']}, {item['role']})")
    print("Statistics families (conditional on saved inputs):")
    for item in plan["statistics_families"]:
        print(f"  {item['family']}: {item['status']}")
    print("Catalog is partial: exact table columns, all raw products and "
          "individual figure IDs are not yet registered.")


def main():
    parser = argparse.ArgumentParser(prog="mdscope")
    sub = parser.add_subparsers(required=True)
    p = sub.add_parser("init")
    p.add_argument("--output", default="analysis.yaml")
    p.set_defaults(func=lambda a: print(interactive_init(a.output)))
    p = sub.add_parser("validate")
    p.add_argument("config")
    p.set_defaults(func=cmd_validate)
    p = sub.add_parser("plan", help="Inspect expected outputs without running trajectories")
    p.add_argument("config")
    p.add_argument("--format", choices=("text", "json"), default="text")
    p.add_argument("--check-inputs", action="store_true",
                   help="Require topology and trajectory paths to exist")
    p.add_argument("--check-topology", action="store_true",
                   help="Inspect topology selections (no trajectory frames)")
    p.set_defaults(func=cmd_plan)
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
    for phase in ("statistics", "plotting"):
        p = sub.add_parser(phase, help=f"Run only the {phase} finalization phase")
        p.add_argument("config")
        p.set_defaults(func=lambda a, phase=phase: print(finalize_project(a.config, phase=phase)))
    p = sub.add_parser("scan-residue-regions", help="Find contiguous changes in saved residue profiles")
    p.add_argument("statistics_dir", type=Path)
    p.set_defaults(func=lambda a: print(scan_saved_residue_regions(a.statistics_dir)))
    p = sub.add_parser("report", help="Summarize bootstrap findings from a completed statistics directory")
    p.add_argument("statistics_dir", type=Path)
    p.set_defaults(func=lambda a: print(generate_report(a.statistics_dir)))
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
