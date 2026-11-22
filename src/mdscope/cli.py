"""
Tiny CLI for creating, validating, running, selecting one replicate, and
finalizing MDscope.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import yaml

from . import __version__
from .config import iter_replicates, load_config, project_root
from .errors import UserInputError


def interactive_init(output="analysis.yaml"):
    print("MDscope configuration generator (analysis runs themselves are non-interactive).")
    name = input("Project name [mdscope_project]: ").strip() or "mdscope_project"
    topology = input("Topology path: ").strip()
    trajectories = [
        x.strip() for x in input("Trajectory paths, comma-separated: ").split(",") if x.strip()
    ]
    ligand = input("Ligand MDAnalysis selection (blank for none): ").strip() or None
    reps = [{"trajectory": path, "label": f"rep{i}"} for i, path in enumerate(trajectories, 1)]
    from .observables.sasa import FREESASA_INSTALL_HINT, import_freesasa
    try:
        import_freesasa()
        sasa_available = True
    except RuntimeError:
        sasa_available = False
        print("WARNING: FreeSASA is not installed, so analyses.sasa is disabled. "
              f"To enable it, run `{FREESASA_INSTALL_HINT}` and set "
              "analyses.sasa.enabled: true.")
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
            "sasa": {"enabled": sasa_available}
        },
    }
    if ligand:
        cfg["analyses"]["ligand_contacts"] = {"enabled": True, "cutoff_A": 4.0}
    Path(output).write_text(yaml.safe_dump(cfg, sort_keys=False))
    return Path(output)


def cmd_init(args):
    try:
        print(interactive_init(args.output))
    except EOFError as exc:
        raise UserInputError(
            "init reads one answer per prompt from stdin, which ended early; "
            "run it in a terminal or pipe all answers"
        ) from exc


def cmd_validate(args):
    from .pipeline.preflight import validate_selections

    cfg = load_config(args.config)
    reps = list(iter_replicates(cfg, args.config))
    print(f"Valid config: {len(cfg['systems'])} systems, {len(reps)} replicates")
    print(f"Resolved output directory: {project_root(cfg, args.config)}")
    for rep in reps:
        print(f"  {rep.system}/{rep.label}: {rep.trajectory}")
    print("\nResolved structural selections:")
    validate_selections(cfg, Path(args.config).resolve().parent, emit=print)


def cmd_plan(args):
    """Print a no-run output plan; topology inspection is opt-in."""
    import json

    from .planning import build_plan

    cfg = load_config(args.config, check_paths=args.check_inputs)
    counts = None
    if args.check_topology:
        from .pipeline.preflight import validate_selections

        counts = validate_selections(cfg, Path(args.config).resolve().parent)
    mode = ("topology_aware" if args.check_topology else
            "input_paths" if args.check_inputs else "schema_only")
    plan = build_plan(cfg, args.config, validation_mode=mode,
                      topology_counts=counts)
    if args.format == "json":
        print(json.dumps(plan, indent=2, sort_keys=True))
        return
    print(f"MDscope plan ({mode}): {len(plan['systems'])} systems, "
          f"{len(plan['replicates'])} replicas")
    print(f"Resolved output directory: {plan['resolved_output_directory']}")
    for pair in plan["inference_design"]["comparisons"]:
        if pair["resolution_kind"] == "exact":
            print(f"  Exact two-sided p floor for {pair['system_1']} vs "
                  f"{pair['system_2']}: {pair['min_attainable_p']:.4g}")
        else:
            print(f"  Monte Carlo p resolution for {pair['system_1']} vs "
                  f"{pair['system_2']}: {pair['monte_carlo_resolution']:.4g} "
                  "(not the exact design floor)")
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


def cmd_run(args):
    from .pipeline.runner import run_all

    print(run_all(args.config))


def cmd_replicate(args):
    from .pipeline.runner import run_one

    print(run_one(args.config, args.index))


def cmd_finalize(args, phase="statistics"):
    from .pipeline.runner import finalize_project

    print(finalize_project(args.config, phase=phase, results_dir=getattr(args, "results_dir", None)))


def cmd_report(args):
    from .reporting.report import generate_report

    print(generate_report(args.statistics_dir))


def cmd_scan_residue_regions(args):
    from .statistics.residue_regions import scan_saved_residue_regions

    cfg = load_config(args.config, check_paths=False) if getattr(args, "config", None) else None
    print(scan_saved_residue_regions(args.statistics_dir, cfg=cfg))


def main():
    parser = argparse.ArgumentParser(prog="mdscope")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--debug", action="store_true",
                        help="Show tracebacks for expected input errors (or set MDSCOPE_DEBUG=1)")
    sub = parser.add_subparsers(required=True)
    p = sub.add_parser("init")
    p.add_argument("--output", default="analysis.yaml")
    p.set_defaults(func=cmd_init)
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
    p.set_defaults(func=cmd_run)
    p = sub.add_parser("replicate", help="run one zero-based systemxreplicate entry")
    p.add_argument("config")
    p.add_argument("--index", type=int, required=True)
    p.set_defaults(func=cmd_replicate)
    p = sub.add_parser("finalize")
    p.add_argument("config")
    p.set_defaults(func=lambda a: cmd_finalize(a, phase="statistics"))
    for phase in ("statistics", "plotting"):
        p = sub.add_parser(phase, help=f"Run only the {phase} finalization phase")
        p.add_argument("config")
        p.add_argument("--results-dir", type=Path,
                       help="Reuse validated existing analysis products from this directory; "
                            "do not schedule trajectory analyses")
        p.set_defaults(func=lambda a, phase=phase: cmd_finalize(a, phase=phase))
    p = sub.add_parser("scan-residue-regions", help="Find contiguous changes in saved residue profiles")
    p.add_argument("statistics_dir", type=Path)
    p.add_argument("--config", help="Use validated YAML discovery settings instead of saved methods")
    p.set_defaults(func=cmd_scan_residue_regions)
    p = sub.add_parser("report", help="Summarize bootstrap findings from a completed statistics directory")
    p.add_argument("statistics_dir", type=Path)
    p.set_defaults(func=cmd_report)
    args = parser.parse_args()
    try:
        args.func(args)
    except UserInputError as exc:
        if args.debug or os.environ.get("MDSCOPE_DEBUG") == "1":
            raise
        parser.exit(2, f"{parser.prog}: error: {exc}\n")


if __name__ == "__main__":
    main()
