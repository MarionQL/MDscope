"""Replica-weighted ceramide scaffold comparison and missingness contracts."""
import numpy as np
import pandas as pd
import pytest

from mdscope.statistics.ceramide_scaffold import (
    scaffold_condition_comparisons,
    summarize_scaffold_replicates,
    write_scaffold_condition_summary,
)

SOURCE = "extensions/sdc4_mechanisms/ceramide_scaffold_pooled.parquet"


def _write_replica(root, system, label, cohort, fractions):
    path = root / system / label
    source = path / SOURCE
    source.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame({
        "time_ps": np.arange(len(fractions), dtype=float) * 10.,
        "cohort": cohort,
        "region": "TM_A",
        "leaflet": "upper",
        "moiety": "ceramide_head",
        "contact_fraction": fractions,
        "matched_local_contact_fraction": np.asarray(fractions) / 2.,
    })
    frame.to_parquet(source, index=False)
    return path


def test_scaffold_contrast_weights_replicas_not_frames_or_molecules(tmp_path):
    systems = {
        "glyco": [
            _write_replica(tmp_path, "glyco", "rep1", "GM1a+BGAL", [.1, .3]),
            _write_replica(tmp_path, "glyco", "rep2", "GM1a+BGAL", [.3, .5]),
        ],
        "no_glyco": [
            _write_replica(tmp_path, "no_glyco", "rep1", "CER160", [.5, .7]),
            _write_replica(tmp_path, "no_glyco", "rep2", "CER160", [.7, .9]),
        ],
    }
    replicas = summarize_scaffold_replicates(systems)
    assert len(replicas) == 4
    assert set(replicas.n_valid_frames) == {2}
    comparisons = scaffold_condition_comparisons(replicas, systems)
    contact = comparisons.loc[comparisons.metric == "contact_fraction"].iloc[0]
    assert contact.status == "ok"
    assert (contact.cohort_1, contact.cohort_2) == ("GM1a+BGAL", "CER160")
    assert contact.mean_1 == pytest.approx(.3)
    assert contact.mean_2 == pytest.approx(.7)
    assert contact.difference_system2_minus_system1 == pytest.approx(.4)
    assert (contact.n_replicates_1, contact.n_replicates_2) == (2, 2)

    cfg = {"analyses": {"sdc4_mechanisms": {"ceramide_scaffold_pooling": True}}}
    assert write_scaffold_condition_summary(cfg, tmp_path, systems) == 2
    outdir = tmp_path / "project_summary/statistics/ceramide_scaffold"
    assert len(pd.read_parquet(outdir / "replicate_means.parquet")) == 4
    assert len(pd.read_parquet(outdir / "condition_comparisons.parquet")) == 2


def test_scaffold_contrast_marks_absent_replicas_unavailable(tmp_path):
    systems = {
        "glyco": [
            _write_replica(tmp_path, "glyco", "rep1", "GM1a+BGAL", [.4]),
            _write_replica(tmp_path, "glyco", "rep2", "GM1a+BGAL", [.6]),
        ],
        "no_glyco": [
            _write_replica(tmp_path, "no_glyco", "rep1", "CER160", [.8]),
            _write_replica(tmp_path, "no_glyco", "rep2", "CER160", []),
        ],
    }
    contrast = scaffold_condition_comparisons(summarize_scaffold_replicates(systems), systems)
    assert set(contrast.status) == {"incomplete_replicas"}
    assert contrast.difference_system2_minus_system1.isna().all()
    assert set(contrast.n_replicates_2) == {0}


def test_scaffold_contrast_does_not_invent_chemical_difference(tmp_path):
    systems = {
        "a": [_write_replica(tmp_path, "a", "rep1", "CER160", [.2])],
        "b": [_write_replica(tmp_path, "b", "rep1", "CER160", [.8])],
    }
    contrast = scaffold_condition_comparisons(summarize_scaffold_replicates(systems), systems)
    assert set(contrast.status) == {"no_chemical_contrast"}
    assert contrast.difference_system2_minus_system1.isna().all()


def test_scaffold_missing_source_and_duplicate_frames_fail(tmp_path):
    path = _write_replica(tmp_path, "a", "rep1", "CER160", [.1, .2])
    source = path / SOURCE
    source.unlink()
    with pytest.raises(FileNotFoundError, match="Required ceramide"):
        summarize_scaffold_replicates({"a": [path]})
    _write_replica(tmp_path, "a", "rep1", "CER160", [.1, .2])
    df = pd.read_parquet(source)
    pd.concat([df, df.iloc[[0]]], ignore_index=True).to_parquet(source, index=False)
    with pytest.raises(ValueError, match="duplicate ceramide"):
        summarize_scaffold_replicates({"a": [path]})


def test_disabled_scaffold_comparison_produces_no_artifacts(tmp_path):
    assert write_scaffold_condition_summary(
        {"analyses": {"sdc4_mechanisms": {"ceramide_scaffold_pooling": False}}},
        tmp_path, {},
    ) == 0
    assert not (tmp_path / "project_summary/statistics/ceramide_scaffold").exists()
