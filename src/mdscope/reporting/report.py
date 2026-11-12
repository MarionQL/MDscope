"""Deterministic review of saved between-condition bootstrap intervals."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyarrow.parquet import read_metadata

from ..core.persistence import atomic_parquet, atomic_text
from ..statistics.coverage import COLUMNS, _identity, _interval_status


_STATE_SOURCES = {
    "full_population_comparisons": ("state",),
    "matched_population_comparisons": ("state",),
    "state_conditioned_nxip_comparisons": ("chain", "state", "endpoint"),
}
_DIFFERENCE = "difference_system2_minus_system1"
_REPORT_FILES = {"findings.parquet", "output_inventory.parquet"}


def _shared_state_coverage(statistics_dir):
    sources = []
    for name, keys in _STATE_SOURCES.items():
        relative = f"../shared_states/{name}.parquet"
        path = statistics_dir / relative
        if not path.is_file():
            continue
        data = pd.read_parquet(path)
        if name == "state_conditioned_nxip_comparisons":
            data = data.rename(columns={
                "system_1": "system_2", "system_0": "system_1",
                "difference_system1_minus_system0": _DIFFERENCE,
                "replicate_bootstrap_ci_low": "difference_ci_low",
                "replicate_bootstrap_ci_high": "difference_ci_high",
            })
            data["ci_level"] = .95
        if data.empty:
            continue
        status = [
            saved if isinstance(saved, str) else _interval_status(low, high)
            for saved, low, high in zip(
                data.get("ci_status", [None] * len(data)),
                data["difference_ci_low"], data["difference_ci_high"],
            )
        ]
        sources.append(data.assign(
            family=f"shared_states/{name}", observable=_identity(data, keys),
            estimator="saved", status=status, source_table=relative,
            ci_low=data["difference_ci_low"], ci_high=data["difference_ci_high"],
        ).reindex(columns=COLUMNS))
    return sources


def _inventory(statistics_dir, manifest, comparison_sources):
    # Inventory numerical files by metadata. Raw trajectories and working
    # caches are intentionally outside the report input.
    roots = [statistics_dir, statistics_dir.parent / "shared_states",
             statistics_dir.parent / "diffusion", statistics_dir.parent / "coupling"]
    if statistics_dir.name == "statistics" and statistics_dir.parent.name == "project_summary":
        project = statistics_dir.parent.parent
        for entry in manifest["replicate_paths"]:
            rel = Path(entry)
            if not rel.is_absolute() and len(rel.parts) == 2 and ".." not in rel.parts:
                roots.extend(project / rel / folder for folder in (
                    "summary", "profiles", "extensions", "interactions",
                ))
    from os.path import relpath
    files = {
        path for root in roots for path in root.rglob("*")
        if path.suffix in {".parquet", ".npz", ".npy", ".json"}
        and path.name not in _REPORT_FILES and path.parent.name != "frames"
    }
    rows = []
    for path in sorted(files):
        relative = relpath(path, statistics_dir)
        rows.append({
            "path": relative,
            "rows": read_metadata(path).num_rows if path.suffix == ".parquet" else None,
            "role": "comparison source" if relative in comparison_sources else "available output",
        })
    return pd.DataFrame(rows, columns=["path", "rows", "role"])


def _cell(value):
    return str(value).replace("|", r"\|").replace("\n", " ")


def _examples(findings):
    """Give each profile or interaction kind one slot before filling the review."""
    examples = findings.copy()
    for name in ("ci_low", "ci_high"):
        examples[name] = examples[name].abs()
    examples["separation"] = examples[["ci_low", "ci_high"]].min(axis=1) / (
        examples[["ci_low", "ci_high"]].max(axis=1)
    )
    examples = examples.sort_values("separation", ascending=False, kind="stable")
    identities = examples["observable"].map(json.loads)
    examples["example_group"] = [
        identity.get("profile", identity.get("interaction_type", row))
        for identity, row in zip(identities, examples["observable"])
    ]
    return examples.drop_duplicates("example_group").head(8).index


def _markdown(coverage, findings, inventory):
    lines = [
        "# MDscope findings", "",
        "Findings have finite bootstrap intervals entirely above or below zero.",
        "These are pointwise interval exclusions. The q column reports an adjusted test where available.",
        "Differences are system 2 minus system 1 in the source units.", "",
        "## Comparison coverage", "",
        "| Family | CI excludes zero | CI includes zero | Unavailable or disabled |",
        "| --- | ---: | ---: | ---: |",
    ]
    for family, rows in coverage.groupby("family", sort=True):
        supported = int(findings.family.eq(family).sum())
        overlapping = int(rows.status.eq("includes_zero").sum())
        lines.append(f"| {_cell(family)} | {supported} | {overlapping} | {len(rows) - supported - overlapping} |")
    gaps = coverage.loc[~coverage.status.isin(["above_zero", "below_zero", "includes_zero"])]
    if not gaps.empty:
        lines.extend(["", "Unavailable and disabled entries by reason", ""])
        lines.extend(f"- {_cell(status)} ({n})" for status, n in gaps.status.value_counts().sort_index().items())
    if findings.empty:
        lines.extend(["", "No saved comparison met the interval criterion."])
    for (left, right, family), rows in findings.groupby(["system_1", "system_2", "family"], sort=True):
        lines.extend([
            "", f"## {_cell(right)} minus {_cell(left)} — {_cell(family)}", "",
            f"{len(rows)} interval exclusions. Representative rows are ordered by separation from zero relative to interval size.",
            "", "| Observable | Estimator | Difference | Bootstrap interval | q |",
            "| --- | --- | ---: | --- | ---: |",
        ])
        for row in rows.loc[_examples(rows)].itertuples(index=False):
            identity = ", ".join(f"{key}={value}" for key, value in json.loads(row.observable).items() if value != "")
            q = f"{row.q_value:.3g}" if pd.notna(row.q_value) else "—"
            level = f"{row.ci_level:.0%} " if pd.notna(row.ci_level) else ""
            lines.append(
                f"| {_cell(identity)} | {_cell(row.estimator)} | {getattr(row, _DIFFERENCE):.5g} | "
                f"{level}[{row.ci_low:.5g}, {row.ci_high:.5g}] | {q} |"
            )
    lines.extend([
        "", "## Saved evidence", "",
        f"The [output inventory](output_inventory.parquet) lists {len(inventory)} saved products.",
        "[Full findings](findings.parquet) contains every selected comparison with its source path.",
        "Use the [coverage audit](coverage/contrast_coverage.parquet) for the complete registered comparison set.",
        "Adjacent positions, alternative estimators, and related metrics can describe the same underlying effect.",
        "Check the coordinate reference of density profiles before interpreting shifts as redistribution.",
    ])
    return "\n".join(lines) + "\n"


def generate_report(statistics_dir: Path) -> Path:
    """Build a findings review without loading trajectories or rendering figures."""
    statistics_dir = Path(statistics_dir)
    manifest = json.loads((statistics_dir / "manifest.json").read_text())
    if manifest.get("status") != "complete":
        raise ValueError("The findings report requires completed statistics")
    coverage = pd.read_parquet(statistics_dir / "coverage/contrast_coverage.parquet")
    extra = _shared_state_coverage(statistics_dir)
    if extra:
        coverage = pd.concat([coverage, *extra], ignore_index=True)
    bounds = coverage[[_DIFFERENCE, "ci_low", "ci_high"]].apply(pd.to_numeric, errors="coerce")
    supported = (
        coverage.status.isin(["above_zero", "below_zero"])
        & np.isfinite(bounds.to_numpy(dtype=float)).all(axis=1)
        & bounds.ci_low.le(bounds.ci_high)
        & (bounds.ci_low.gt(0) | bounds.ci_high.lt(0))
    )
    findings = coverage.loc[supported].sort_values(
        ["family", "system_1", "system_2", "observable", "estimator"], ignore_index=True,
    )
    sources = set(coverage.loc[coverage.observable.notna(), "source_table"])
    inventory = _inventory(statistics_dir, manifest, sources)
    markdown = _markdown(coverage, findings, inventory)
    atomic_parquet(statistics_dir / "findings.parquet", findings)
    atomic_parquet(statistics_dir / "output_inventory.parquet", inventory)
    report = statistics_dir / "findings.md"
    atomic_text(report, markdown)
    return report
