"""Cross-metric region-screen fixtures: missingness, boundaries, direction and CI."""
import numpy as np
import pandas as pd
import pytest

from mdscope.statistics.residue_response import residue_response_screen


def _inputs():
    reps, contrasts = [], []
    for profile, sign in (("sasa_mean__isolated", 1.), ("residue_occlusion__any", -1.)):
        for residue in range(10, 17):
            label = f"PROA:ALA{residue}"
            for system, shift in (("control", 0.), ("treated", sign * (2. if residue < 14 else -2.))):
                for replicate, noise in (("rep1", -.1), ("rep2", 0.), ("rep3", .1)):
                    reps.append(dict(system=system, replicate=replicate, profile=profile,
                                     position_label=label, value=shift + noise))
            difference = sign * (2. if residue < 14 else -2.)
            contrasts.append(dict(system_1="control", system_2="treated", profile=profile,
                                  position_label=label,
                                  difference_system2_minus_system1=difference,
                                  difference_ci_low=difference - .5,
                                  difference_ci_high=difference + .5))
    return pd.DataFrame(reps), pd.DataFrame(contrasts)


def _screen(reps, contrasts):
    return residue_response_screen(
        reps, contrasts,
        {"control": ["rep1", "rep2", "rep3"],
         "treated": ["rep1", "rep2", "rep3"]},
    )


def test_adjacent_multimetric_changes_and_direction_split():
    reps, contrasts = _inputs()
    matrix, regions = _screen(reps, contrasts)
    assert len(matrix) == 14
    assert len(regions) == 2
    assert [(int(r.start_residue), int(r.end_residue))
            for r in regions.itertuples()] == [(10, 13), (14, 16)]
    assert set(regions.n_metrics) == {2}
    assert matrix.ci_excludes_zero.all()
    assert matrix.loc[matrix.profile == "residue_occlusion__any",
                      "standardized_effect"].iloc[0] < 0


def test_missing_replica_is_missing_not_zero_or_a_discovered_region():
    reps, contrasts = _inputs()
    reps = reps.loc[~((reps.system == "control") &
                      (reps.replicate == "rep3") &
                      (reps.position_label == "PROA:ALA12"))]
    matrix, regions = _screen(reps, contrasts)
    assert not (matrix.resid == 12).any()
    assert not ((regions.start_residue <= 12) & (regions.end_residue >= 12)).any()


def test_adjacent_numbers_across_chains_are_not_joined():
    reps, contrasts = _inputs()
    reps["position_label"] = reps.position_label.str.replace("PROA:", "PROB:")
    contrasts["position_label"] = contrasts.position_label.str.replace("PROA:", "PROB:")
    _, regions = _screen(reps, contrasts)
    assert set(regions.chain) == {"PROB"}


def test_duplicate_labels_fail_ambiguously_mapped_inputs():
    reps, contrasts = _inputs()
    with pytest.raises(ValueError, match="Ambiguous"):
        _screen(pd.concat([reps, reps.iloc[[0]]], ignore_index=True), contrasts)


def test_missing_bootstrap_interval_still_descriptive():
    reps, contrasts = _inputs()
    contrasts.loc[contrasts.position_label == "PROA:ALA10",
                  ["difference_ci_low", "difference_ci_high"]] = np.nan
    matrix, regions = _screen(reps, contrasts)
    assert not matrix.loc[matrix.resid == 10, "ci_excludes_zero"].any()
    assert not regions.empty
    assert regions.exploratory.all()
