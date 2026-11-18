#!/usr/bin/env python3
"""Render execution and import diagrams from the current MDscope checkout.

Execution edges summarize runner/streaming/fused_all control flow; they are
maintained documentation, not runtime traces. Import edges come from the
existing AST inventory in scripts/architecture_map.py.
No MDscope package or scientific data are imported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys


FOCUS = {
    "mdscope.cli", "mdscope.config", "mdscope.planning", "mdscope.registry",
    "mdscope.output_paths", "mdscope.statistics_families", "mdscope.statistics",
    "mdscope.pipeline.runner", "mdscope.pipeline.streaming",
    "mdscope.pipeline.fused_all", "mdscope.pipeline.fused_membrane",
    "mdscope.pipeline.fused_extensions", "mdscope.pipeline.fused_distance",
    "mdscope.pipeline.replica_tasks", "mdscope.pipeline.analysis_selection",
    "mdscope.observables.lipid_enrichment", "mdscope.observables.solvation",
    "mdscope.observables.solvation_postprocess", "mdscope.observables.solvation_common",
    "mdscope.mechanisms.sdc4_mechanisms",
}

# Conditions are represented with dashed edges. The initial streaming
# trajectory loop and observer-only loop are mutually exclusive; later
# recovery and specialized work may require additional trajectory traversals.
EXECUTION = {
    "nodes": {
        "entry": "CLI: run / replicate",
        "runner": "pipeline.runner: run_replicate",
        "precheck": "Validate existing products before writes",
        "choice": "Core streaming or aligned cache required?",
        "stream": "streaming: shared primary frame pass",
        "fusedonly": "fused_all: pending instances, one shared pass",
        "membrane": "FusedMembrane: order / normal / enrichment / contacts",
        "extensions": "FusedExtensions: H-bonds / SASA / DSSP / solvation / hydration",
        "distance": "FusedDistances: residue-membrane",
        "sdc4": "SDC4Mechanisms: motif / HS / membrane",
        "core": "Core metrics: RMSD / RMSF / Rg / density / membrane",
        "audit": "Per-family missing/incomplete product checks",
        "fallback": "Selective fused-family fallback readers",
        "tasks": "Conditional specialized tasks / multi-origin MSD",
        "verify": "Validate registered artifact contracts / stage markers",
        "summary": "Saved per-replica data and summary",
        "late": "Conditional density completion / hydration profiles",
        "stats": "Offline statistics, findings and figures",
    },
    "edges": [
        ("entry", "runner", False), ("runner", "precheck", False),
        ("precheck", "choice", False),
        ("choice", "stream", True), ("choice", "fusedonly", True),
        ("stream", "core", False),
        ("stream", "membrane", False), ("stream", "extensions", False),
        ("stream", "distance", False), ("stream", "sdc4", False),
        ("fusedonly", "membrane", True), ("fusedonly", "extensions", True),
        ("fusedonly", "distance", True), ("fusedonly", "sdc4", True),
        ("core", "audit", False),
        ("membrane", "audit", False), ("extensions", "audit", False),
        ("distance", "audit", False), ("sdc4", "audit", False),
        ("choice", "audit", True),
        ("audit", "fallback", True), ("fallback", "verify", False),
        ("audit", "verify", False),
        ("verify", "tasks", True), ("tasks", "summary", False),
        ("verify", "summary", False),
        ("summary", "late", True), ("late", "stats", False),
        ("summary", "stats", False),
    ],
}



# A conceptual decision tree for artifact reuse and reconstruction.
# Fused-only replays use per-instance selection; existing corrupt products
# fail preflight rather than being silently replaced.
RECOVERY = {
    "nodes": {
        "start": "Start replicate / read output registry",
        "preflight": "Validate existing primary and required companions",
        "invalid": "Invalid existing artifact?",
        "abort": "Fail; preserve existing outputs for inspection",
        "stream": "Core stream required?",
        "checkpoint": "Stream checkpoint state?",
        "reuse_stream": "Reuse stream; validate checkpoint inventory",
        "replay_stream": "Resume or execute primary stream using checked chunks",
        "pending": "Pending fused artifacts? (analysis + run ID)",
        "skip": "No pending fused output; avoid opening trajectory",
        "shared": "Observer-only: one source-frame pass for selected instances",
        "scratch": "Isolated scratch chunks; publish selected products only",
        "fallback": "Per-family selective fallback if stream left missing products",
        "validate": "Validate registered artifacts and their companions",
        "summary": "Reuse or rebuild replica summary as required",
        "extras": "Optional specialized tasks / further passes",
        "end": "Saved replica ready for offline finalization",
    },
    "edges": [
        ("start", "preflight", ""),
        ("preflight", "invalid", ""),
        ("invalid", "abort", "yes"),
        ("invalid", "stream", "no"),
        ("stream", "checkpoint", "yes"),
        ("checkpoint", "reuse_stream", "completed + valid"),
        ("checkpoint", "replay_stream", "absent / resumable"),
        ("checkpoint", "abort", "invalid / inconsistent"),
        ("reuse_stream", "fallback", ""),
        ("replay_stream", "fallback", ""),
        ("stream", "pending", "no"),
        ("pending", "skip", "no"),
        ("pending", "shared", "yes"),
        ("shared", "scratch", ""),
        ("scratch", "validate", ""),
        ("skip", "validate", ""),
        ("fallback", "validate", ""),
        ("validate", "summary", ""),
        ("summary", "extras", ""),
        ("extras", "end", ""),
    ],
}


def quote(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def write_execution(path: Path) -> None:
    lines = [
        "digraph execution {",
        '  graph [rankdir=TB, bgcolor=white, label="MDscope execution (conceptual, conditional)", labelloc=t];',
        '  node [shape=box, style="rounded,filled", fillcolor="#e6effb", color="#7299cc"];',
    ]
    for key, label in EXECUTION["nodes"].items():
        attrs = ['label=' + quote(label)]
        if key in {"stream", "fusedonly", "membrane", "extensions", "distance", "sdc4"}:
            attrs.append('fillcolor="#dcf4e8"')
        elif key in {"choice", "audit", "fallback", "tasks", "late"}:
            attrs.append('fillcolor="#fff1d9"')
        lines.append(f'  {quote(key)} [{", ".join(attrs)}];')
    for source, target, conditional in EXECUTION["edges"]:
        style = ' [style=dashed, color="#c88a3b"]' if conditional else ""
        lines.append(f'  {quote(source)} -> {quote(target)}{style};')
    lines.append("}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_recovery(path: Path) -> None:
    lines = [
        "digraph recovery {",
        '  graph [rankdir=TB, bgcolor=white, label="MDscope artifact recovery (conceptual)", labelloc=t];',
        '  node [shape=box, style="rounded,filled", fillcolor="#e6effb", color="#7299cc"];',
    ]
    for key, label in RECOVERY["nodes"].items():
        attrs = ["label=" + quote(label)]
        if key in {"invalid", "stream", "checkpoint", "pending"}:
            attrs.append('shape=diamond')
            attrs.append('fillcolor="#fff1d9"')
        elif key in {"abort"}:
            attrs.append('fillcolor="#ffe1df"')
        elif key in {"reuse_stream", "skip", "shared", "scratch", "replay_stream"}:
            attrs.append('fillcolor="#dcf4e8"')
        lines.append(f'  {quote(key)} [{", ".join(attrs)}];')
    for source, target, decision in RECOVERY["edges"]:
        attrs = f' [label={quote(decision)}]' if decision else ""
        lines.append(f'  {quote(source)} -> {quote(target)}{attrs};')
    lines.append("}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_imports(path: Path, inventory: dict, focused: bool) -> None:
    included = (set(inventory["modules"]) & FOCUS) if focused else set(inventory["modules"])
    phases = {}
    for item in inventory["import_sites"]:
        phases.setdefault((item["source"], item["target"]), set()).add(item["phase"])
    cycle_edges = {
        (source, target)
        for group in inventory["cycles"] for source in group for target in group
        if source != target
    }
    lines = [
        "digraph imports {",
        '  graph [rankdir=LR, bgcolor=white, overlap=false, splines=true, label="MDscope static imports (current checkout)", labelloc=t];',
        '  node [shape=box, style="rounded,filled", fillcolor="#ecf2fc", color="#7299cc"];',
    ]
    for module in sorted(included):
        lines.append(f'  {quote(module)} [label={quote(module.removeprefix("mdscope."))}];')
    for source, target in inventory["imports"]:
        if source not in included or target not in included:
            continue
        deferred = "module_load" not in phases.get((source, target), set())
        attrs = ["style=dashed" if deferred else "style=solid"]
        if (source, target) in cycle_edges:
            attrs += ['color="#906bd0"', "penwidth=2.2"]
        lines.append(f'  {quote(source)} -> {quote(target)} [{", ".join(attrs)}];')
    lines.append("}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inventory", type=Path,
                        help="Existing architecture.json; generated from this checkout if omitted")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    inventory_path = args.inventory or args.output / "architecture.json"
    if args.inventory is None:
        subprocess.run(
            [sys.executable, str(root / "scripts" / "architecture_map.py"),
             "--output", str(args.output)], check=True,
        )
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    if not {"modules", "imports", "import_sites", "cycles"} <= inventory.keys():
        parser.error("Invalid architecture inventory")
    write_execution(args.output / "execution.dot")
    write_recovery(args.output / "recovery.dot")
    write_imports(args.output / "imports_focused.dot", inventory, focused=True)
    write_imports(args.output / "imports_full.dot", inventory, focused=False)
    dot = shutil.which("dot")
    if dot:
        for name in ("execution", "recovery", "imports_focused", "imports_full"):
            subprocess.run(
                [dot, "-Tsvg", str(args.output / (name + ".dot")),
                 "-o", str(args.output / (name + ".svg"))], check=True,
            )
    print(f"Graphs: {args.output} (SVG: {'yes' if dot else 'requires Graphviz dot'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
