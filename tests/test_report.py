"""Findings include only valid CI exclusions and preserve saved evidence."""
import json

import pandas as pd
import pytest

from mdscope.reporting.report import generate_report
from mdscope.statistics.coverage import COLUMNS


def test_findings_from_saved_coverage_and_shared_states_are_deterministic(tmp_path):
    project_summary = tmp_path / "project_summary"
    stats = project_summary / "statistics"
    (stats / "coverage").mkdir(parents=True)
    manifest = stats / "manifest.json"
    manifest.write_text(json.dumps({"status": "complete", "replicate_paths": ["a/rep1"]}))
    base = dict(family="profiles", estimator="saved", system_1="a", system_2="b",
                source_table="profiles/system_comparisons.parquet", ci_level=.90)
    rows = []
    for name, diff, low, high, status in [
        ("strict_positive", .4, .2, .6, "above_zero"),
        ("strict_negative", -.4, -.6, -.2, "below_zero"),
        ("at_zero", .1, 0., .2, "includes_zero"),
        ("missing", .1, None, None, "incomplete_replica_coverage"),
        ("invalid", .4, .5, .3, "above_zero"),
    ]:
        rows.append(dict(base, observable=json.dumps({"profile": name}),
                         difference_system2_minus_system1=diff, ci_low=low, ci_high=high,
                         status=status))
    coverage = stats / "coverage/contrast_coverage.parquet"
    pd.DataFrame(rows, columns=COLUMNS).to_parquet(coverage)
    original = coverage.read_bytes()
    source = stats / base["source_table"]
    source.parent.mkdir()
    pd.DataFrame({"position": [1, 2]}).to_parquet(source)
    replica_dir = tmp_path / "a/rep1/profiles"
    replica_dir.mkdir(parents=True)
    pd.DataFrame({"value": [2.]}).to_parquet(replica_dir / "saved.parquet")

    shared = project_summary / "shared_states"
    shared.mkdir()
    pd.DataFrame([dict(chain="A", state=1, endpoint="hs_contacts", system_0="a", system_1="b",
                       difference_system1_minus_system0=.25, replicate_bootstrap_ci_low=.1,
                       replicate_bootstrap_ci_high=.4)]).to_parquet(
        shared / "state_conditioned_nxip_comparisons.parquet",
    )
    result = generate_report(stats)
    first = result.read_bytes()
    assert b"strict_positive" in first and b"strict_negative" in first
    assert b"at_zero" not in first and b"profile=invalid" not in first
    assert b"incomplete_replica_coverage" in first and b"90% [0.2, 0.6]" in first
    findings = pd.read_parquet(stats / "findings.parquet")
    assert len(findings) == 3
    assert set(findings.system_2) == {"b"}
    inventory = pd.read_parquet(stats / "output_inventory.parquet").set_index("path")
    assert inventory.loc[base["source_table"], "role"] == "comparison source"
    assert inventory.loc["../../a/rep1/profiles/saved.parquet", "rows"] == 1
    assert generate_report(stats).read_bytes() == first
    assert coverage.read_bytes() == original

    manifest.write_text(json.dumps({"status": "running"}))
    with pytest.raises(ValueError, match="completed statistics"):
        generate_report(stats)
    assert result.read_bytes() == first
