"""Summarize persisted MDscope end-to-end profiling without reading trajectories.

Usage:
  python scripts/benchmark_workflow.py results/my_project
  python scripts/benchmark_workflow.py results/my_project --output benchmark.csv
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def summarize(root):
    root = Path(root)
    rows = []
    for path in sorted(root.glob("*/*/streaming_phase_profile.json")):
        data = json.loads(path.read_text())
        for phase in data.get("timings", []):
            rows.append({"system": path.parent.parent.name,
                         "replicate": path.parent.name,
                         "scope": "streaming_phase",
                         "stage": phase["stage"],
                         "seconds": float(phase["seconds"]),
                         "frames": int(data.get("frames_completed", 0)),
                         "seconds_per_frame": phase.get("seconds_per_frame")})
    for path in sorted(root.glob("*/*/stage_timings.json")):
        data = json.loads(path.read_text())
        for phase in data.get("stages", []):
            rows.append({"system": path.parent.parent.name,
                         "replicate": path.parent.name,
                         "scope": "replicate_stage",
                         "stage": phase["stage"],
                         "seconds": float(phase["seconds"]),
                         "frames": None, "seconds_per_frame": None})
    path = root / "project_summary" / "finalization_timings.json"
    if path.exists():
        data = json.loads(path.read_text())
        for phase in data.get("stages", []):
            rows.append({"system": "project", "replicate": "all",
                         "scope": "finalization",
                         "stage": phase["stage"],
                         "seconds": float(phase["seconds"]),
                         "frames": None, "seconds_per_frame": None})
    return sorted(rows, key=lambda row: row["seconds"], reverse=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = summarize(args.root)
    if not rows:
        parser.error("No profiling JSON found. Set MDSCOPE_PROFILE=1 for replicate jobs.")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    for scope in ("replicate_stage", "streaming_phase", "finalization"):
        group = [row for row in rows if row["scope"] == scope]
        if not group:
            continue
        print(f"\\n{scope} — slowest stages")
        for row in group[:20]:
            rate = (f" | {row['seconds_per_frame']:.5f} s/frame"
                    if row["seconds_per_frame"] is not None else "")
            print(f"  {row['system']}/{row['replicate']} "
                  f"{row['stage']}: {row['seconds']:.2f} s{rate}")
    print("\\nScopes overlap: streaming phases are contained within replicate stages. "
          "Do not sum them together or interpret the ranking as CPU profiling.")


if __name__ == "__main__":
    main()
