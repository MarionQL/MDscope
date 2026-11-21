"""Deterministic review of saved between-condition bootstrap intervals."""
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from pyarrow.parquet import read_metadata

from ..core.persistence import atomic_parquet, atomic_text
from ..statistics.coverage import COLUMNS, SOURCES, _identity, _interval_status


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
        for entry in manifest.get("replicate_paths", []):
            rel = Path(entry)
            if not rel.is_absolute() and len(rel.parts) == 2 and ".." not in rel.parts:
                roots.extend(project / rel / folder for folder in (
                    "summary", "profiles", "extensions", "interactions",
                ))
    from os.path import relpath
    files = {
        path for root in roots for path in root.rglob("*")
        if path.suffix in {".parquet", ".npz", ".npy", ".json"}
        and path.name not in _REPORT_FILES and path.parent.name not in {"frames", "replayed_frames"}
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


def _replicate_support(statistics_dir):
    """Count saved independent replicate identities for each observable.

    The manifest lists scheduled inputs; only per-replicate statistics rows
    establish that a particular observable has data in a particular replica.
    """
    support = {}
    for family, folder, replica_file, _, keys, _ in SOURCES:
        path = statistics_dir / folder / f"{replica_file}.parquet"
        if not path.is_file():
            continue
        columns = ["system", "replicate", *keys]
        try:
            table = pd.read_parquet(path, columns=columns)
        except (KeyError, ValueError) as exc:
            raise ValueError(f"Cannot determine replicate support from {path}: {exc}") from exc
        if table.empty:
            continue
        table["observable"] = _identity(table, keys)
        for (system, observable), rows in table.groupby(["system", "observable"], dropna=False):
            support[family, observable, system] = set(rows.replicate.dropna().astype(str))
    # Matched time-series data may reuse the saved ordinary time-series
    # replicate table; do not assign counts from other metric families.
    for (family, observable, system), labels in tuple(support.items()):
        if family == "time_series":
            support["matched_time_series", observable, system] = labels
    return support


def _support_floor(left_labels, right_labels, methods):
    """Resolution for the available replica identities of one observable."""
    tests = methods.get("hypothesis_tests", {}) or {}
    design = tests.get("design", "independent")
    if left_labels is None or right_labels is None:
        return float("nan"), "unknown"
    a, b = len(left_labels), len(right_labels)
    if not a or not b:
        return float("nan"), "unknown"
    total = (2 ** len(left_labels & right_labels) if design == "paired"
             else math.comb(a + b, a))
    if total < 2:
        return float("nan"), "unknown"
    if total <= int(tests.get("exact_max_assignments", 100000)):
        return 2 / total, "exact"
    return 1 / (int(tests.get("permutation_draws", 10000)) + 1), "monte_carlo_resolution"


def _design_resolution(manifest):
    """Minimum attainable exact two-sided p, or a conservative draw-resolution note."""
    counts = Counter()
    for entry in manifest.get("replicate_paths", []):
        parts = Path(entry).parts
        if len(parts) == 2 and ".." not in parts:
            counts[parts[0]] += 1
    methods = manifest.get("methods", {})
    tests = methods.get("hypothesis_tests", {}) or {}
    design = tests.get("design", "independent")
    limit = int(tests.get("exact_max_assignments", 100000))
    draws = int(tests.get("permutation_draws", 10000))
    def floor(left, right):
        na, nb = counts[left], counts[right]
        if not na or not nb:
            return float("nan"), "unknown"
        if design == "paired":
            # Pair identities are required: count matched labels, not the smaller size.
            paths = {system: {Path(s).parts[1] for s in manifest.get("replicate_paths", [])
                              if len(Path(s).parts) == 2 and Path(s).parts[0] == system}
                     for system in (left, right)}
            n = len(paths[left] & paths[right])
            assignments = 2 ** n
            minimum = 2 / assignments if n else float("nan")
        else:
            assignments = math.comb(na + nb, na)
            minimum = 2 / assignments
        if assignments > limit:
            # Monte Carlo's +1 correction has resolution 1/(draws+1);
            # this is NOT a guaranteed attainable exact p-value.
            return 1 / (draws + 1), "monte_carlo_resolution"
        return minimum, "exact"
    return counts, floor


def _markdown(coverage, findings, inventory, *, design=None, alpha=.05):
    lines = [
        "# MDscope findings", "",
        "Findings have finite bootstrap intervals entirely above or below zero.",
        "These are pointwise interval exclusions. The q column reports an adjusted test where available.",
        "Replicate counts refer to observed per-replica metric rows, not scheduled runs.",
        "Min p denotes the exact two-sided design floor; (MC) denotes only Monte Carlo draw resolution.",
        "Differences are system 2 minus system 1 in the source units.", "",
        "## Comparison coverage", "",
        "| Family | CI excludes zero | CI includes zero | Unavailable or disabled |",
        "| --- | ---: | ---: | ---: |",
    ]
    if design is not None:
        _, floor = design
        comparisons = sorted({
            (left, right) for left, right in zip(coverage.system_1, coverage.system_2)
            if isinstance(left, str) and isinstance(right, str)
        })
        warnings = []
        for left, right in comparisons:
            value, method = floor(left, right)
            if method == "exact" and np.isfinite(value) and value > alpha:
                warnings.append(
                    f"{left} vs {right}: smallest attainable exact two-sided p = {value:.4g} "
                    f"(alpha = {alpha:.4g})."
                )
        if warnings:
            lines[5:5] = ["**Design limitation:** " + " ".join(warnings),
                          "Pointwise bootstrap exclusions do not imply BH-adjusted significance.", ""]
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
            "", "| Observable | Estimator | Replicates | Min p | Difference | Bootstrap interval | q |",
            "| --- | --- | --- | ---: | ---: | --- | ---: |",
        ])
        for row in rows.loc[_examples(rows)].itertuples(index=False):
            identity = ", ".join(f"{key}={value}" for key, value in json.loads(row.observable).items() if value != "")
            q = f"{row.q_value:.3g}" if pd.notna(row.q_value) else "—"
            level = f"{row.ci_level:.0%} " if pd.notna(row.ci_level) else ""
            rep_count = (
                f"{int(row.n_replicates_system_1)} vs {int(row.n_replicates_system_2)}"
                if pd.notna(row.n_replicates_system_1) and pd.notna(row.n_replicates_system_2)
                else "—"
            )
            minimum = f"{row.min_attainable_p:.4g}" if pd.notna(row.min_attainable_p) else "—"
            if row.p_resolution_kind == "monte_carlo_resolution":
                minimum += " (MC)"
            lines.append(
                f"| {_cell(identity)} | {_cell(row.estimator)} | {rep_count} | {minimum} | {getattr(row, _DIFFERENCE):.5g} | "
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
    from ..statistics_manifest import read_statistics_manifest
    manifest_path = statistics_dir / "manifest.json"
    facts = read_statistics_manifest(manifest_path)
    if facts.state == "missing":
        raise FileNotFoundError(
            f"The findings report requires completed statistics: missing {manifest_path}. "
            "Run 'mdscope statistics' before generating the report."
        )
    if facts.state == "invalid":
        raise ValueError(
            f"The findings report requires completed statistics: invalid {manifest_path}. "
            "Inspect the saved manifest before retrying."
        )
    if not facts.complete:
        raise ValueError(
            f"The findings report requires completed statistics: "
            f"{manifest_path} has no complete status."
        )
    manifest = facts.payload
    coverage_path = statistics_dir / "coverage/contrast_coverage.parquet"
    if not coverage_path.is_file():
        raise FileNotFoundError(
            f"Missing saved contrast coverage: {coverage_path}. "
            "Complete statistics before generating findings."
        )
    coverage = pd.read_parquet(coverage_path)
    coverage = coverage.loc[~coverage.family.map(facts.family_disabled)].copy()
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
    counts, floor = _design_resolution(manifest)
    support = _replicate_support(statistics_dir)
    methods = manifest.get("methods", {})
    for table in (coverage, findings):
        left_sets = [support.get((family, observable, left))
                     for family, observable, left in zip(table.family, table.observable, table.system_1)]
        right_sets = [support.get((family, observable, right))
                      for family, observable, right in zip(table.family, table.observable, table.system_2)]
        table["n_replicates_system_1"] = [len(labels) if labels is not None else None
                                            for labels in left_sets]
        table["n_replicates_system_2"] = [len(labels) if labels is not None else None
                                            for labels in right_sets]
        resolutions = [_support_floor(a, b, methods) for a, b in zip(left_sets, right_sets)]
        table["min_attainable_p"] = [value for value, _ in resolutions]
        table["p_resolution_kind"] = [kind for _, kind in resolutions]
    sources = set(coverage.loc[coverage.observable.notna(), "source_table"])
    inventory = _inventory(statistics_dir, manifest, sources)
    alpha = manifest.get("methods", {}).get("hypothesis_tests", {}).get("alpha", .05)
    markdown = _markdown(coverage, findings, inventory, design=(counts, floor), alpha=float(alpha))
    atomic_parquet(statistics_dir / "findings.parquet", findings)
    atomic_parquet(statistics_dir / "output_inventory.parquet", inventory)
    report = statistics_dir / "findings.md"
    atomic_text(report, markdown)
    return report
