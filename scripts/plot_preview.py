#!/usr/bin/env python3
"""Render deterministic fake MD results through the real MDscope plotters.

Usage: python scripts/plot_preview.py --output /tmp/mdscope-plot-preview
No XTC/TPR, statistics calculation, or scientific outputs are required.
Run this *same script* with alternate PYTHONPATHs to compare renderer versions.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

from mdscope.reporting.plotting import plot_all

STAMP = "SYNTHETIC EXAMPLE — NOT MD RESULTS"


def create_fixture(root: Path):
    """Write small, deterministic examples of the *existing* saved-data schema."""
    times = np.linspace(0, 120_000, 61)
    system_repouts = {}
    for system_index, system in enumerate(("condition_A", "condition_B")):
        paths = []
        for rep_index in range(2):
            rep = f"rep{rep_index + 1}"
            repout = root / system / rep
            summary = repout / "summary"
            summary.mkdir(parents=True, exist_ok=False)
            phase = 0.22 * rep_index + 0.4 * system_index
            t = times / 1000.0
            pd.DataFrame({
                "time_ps": times,
                "rmsd_A": (0.5 + 0.15 * system_index
                           + 0.009 * t + 0.16 * np.sin(t / 11 + phase)),
                "radius_of_gyration_A": (14.0 + 0.4 * system_index
                                       + 0.15 * np.sin(t / 17 + phase)),
                "membrane_thickness_A": (31.0 + 1.2 * system_index
                                         + 0.6 * np.cos(t / 15 + phase)),
            }).to_parquet(summary / "timeseries.parquet", index=False)
            by_group = []
            for group_index, group in enumerate(("chain_A", "chain_B")):
                by_group.append(pd.DataFrame({
                    "time_ps": times,
                    "group": group,
                    "rmsd_A": (0.3 + 0.25 * group_index
                               + 0.008 * t + 0.2 * np.sin(t / 10 + phase + group_index)),
                }))
            pd.concat(by_group, ignore_index=True).to_parquet(
                summary / "group_timeseries.parquet", index=False
            )
            paths.append(repout)
        system_repouts[system] = paths
    cfg = {
        "project": {"name": "synthetic_figure_preview"},
        "systems": {
            name: {"label": f"Condition {i + 1}",
                   "replicates": [{"label": f"rep{j+1}"} for j in range(2)]}
            for i, name in enumerate(system_repouts)
        },
        "plots": {
            "enabled": True,
            "families": {"time_series": True, "profiles": False, "interactions": False},
            "distributions": False,
            "autocorrelation": False,
            "bootstrap_tables": False,
            "dpi": 120,
        },
        "statistics": {"enabled": False},
    }
    return cfg, system_repouts


def mark_synthetic(path: Path):
    """Mark the preview copy, not the underlying production renderer."""
    with Image.open(path) as original:
        frame = original.convert("RGB")
    draw = ImageDraw.Draw(frame)
    draw.rectangle((7, 7, 273, 27), fill="white", outline="black", width=1)
    draw.text((12, 12), "SYNTHETIC - NOT MD RESULTS", fill="black")
    frame.save(path)


def render_preview(root: Path):
    root = Path(root)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Refusing to replace preview or results in nonempty {root}")
    root.mkdir(parents=True, exist_ok=True)
    cfg, system_repouts = create_fixture(root)
    plot_all(cfg, root, system_repouts)  # Production renderer, no plot duplication.
    figures = sorted((root / "project_summary/figures").rglob("*.png"))
    if len(figures) < 4:
        raise RuntimeError(f"Expected at least four representative plot previews; got {len(figures)}")
    for path in figures:
        mark_synthetic(path)
    import mdscope
    record = {
        "notice": STAMP,
        "renderer_module": str(Path(mdscope.__file__).resolve()),
        "input_data_kind": "deterministic synthetic Parquet",
        "plots": [str(path.relative_to(root)) for path in figures],
    }
    (root / "preview_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = render_preview(args.output)
    print(f"Rendered {len(record['plots'])} marked synthetic previews in {args.output}")


if __name__ == "__main__":
    main()
