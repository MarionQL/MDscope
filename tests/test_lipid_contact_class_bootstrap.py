"""Independent arithmetic and end-to-end saved-data checks for lipid classes."""
import numpy as np
import pandas as pd
import pytest

from mdscope.observables.lipid_contacts import contact_matrix
from mdscope.statistics.lipid_contacts import (
    _classes, contact_frames, analyze_lipid_contacts,
)


CLASSES = {
    "PC": ["POPC", "DOPC"],
    "SM": ["PSM"],
    "glycolipid": ["GM1a"],
}
SPECIES = ["POPC", "POPC", "DOPC", "PSM",
           "POPC", "POPC", "DOPC", "PSM"]
UPPER = [True] * 4 + [False] * 4


def _frame(time=0, extra=True):
    upper = [(0, 1), (0, 2), (1, 2), (2, 3)] if extra else [(0, 1), (2, 3)]
    lower = [(4, 5), (4, 6)]
    result = contact_matrix(
        SPECIES, UPPER, upper + lower, categories=[*CLASSES["PC"], "PSM", "GM1a"],
    )
    result.insert(0, "time_ps", float(time))
    return result


def test_class_counts_deduplicate_species_pairs_and_recompute_denominators():
    result = _classes(_frame(), {species: group for group, members in CLASSES.items()
                                 for species in members})
    upper = result[result.leaflet == "upper"]
    pc_pc = upper[(upper.species_a == "PC") & (upper.species_b == "PC")].iloc[0]
    assert pc_pc.pair_count_undirected == 3  # 1 POPC-POPC + 2 POPC-DOPC
    assert pc_pc.neighbors_per_molecule == 2  # 2*3 / 3 PC lipids
    assert pc_pc.enrichment == pytest.approx(1.5)  # 3/(4*3/6)
    pc_sm = upper[(upper.species_a == "PC") & (upper.species_b == "SM")].iloc[0]
    sm_pc = upper[(upper.species_a == "SM") & (upper.species_b == "PC")].iloc[0]
    assert pc_sm.pair_count_undirected == sm_pc.pair_count_undirected == 1
    assert pc_sm.neighbors_per_molecule == pytest.approx(1/3)
    assert sm_pc.neighbors_per_molecule == pytest.approx(1.)
    assert pc_sm.enrichment == pytest.approx(.5)
    absent = upper[(upper.species_a == "glycolipid") &
                   (upper.species_b == "PC")].iloc[0]
    assert np.isnan(absent.pair_count_undirected)
    assert np.isnan(absent.enrichment)
    lower = result[(result.leaflet == "lower") &
                   (result.species_a == "PC") & (result.species_b == "PC")].iloc[0]
    assert lower.pair_count_undirected == 2


def test_all_species_and_class_counts_remain_symmetric():
    combined = contact_frames(_frame(), CLASSES)
    for level, group in combined.groupby("level"):
        for _, row in group.iterrows():
            reverse = group[(group.leaflet == row.leaflet) &
                            (group.species_a == row.species_b) &
                            (group.species_b == row.species_a)]
            assert len(reverse) == 1
            assert (np.isnan(row.pair_count_undirected) and
                    np.isnan(reverse.iloc[0].pair_count_undirected) or
                    row.pair_count_undirected == reverse.iloc[0].pair_count_undirected)


def test_overlapping_or_unmapped_lipid_classes_fail():
    with pytest.raises(ValueError, match="multiple classes"):
        contact_frames(_frame(), {"one": ["POPC", "DOPC"],
                                  "two": ["POPC", "PSM"]})
    with pytest.raises(ValueError, match="missing species"):
        contact_frames(_frame(), {"only": ["POPC"]})


def test_corrupt_directional_counts_fail():
    table = _frame()
    row = table[(table.leaflet == "upper") &
                (table.species_a == "POPC") &
                (table.species_b == "DOPC")].index[0]
    table.loc[row, "pair_count_undirected"] += 10
    with pytest.raises(ValueError, match="AB/BA"):
        contact_frames(table, CLASSES)


def test_bootstrap_exports_complete_replicas_leaflets_and_condition_differences(tmp_path):
    reps = {}
    for system, richer in (("control", False), ("treated", True)):
        reps[system] = []
        for i in range(3):
            out = tmp_path / system / f"rep{i + 1}"
            (out / "extensions/lipid_contacts").mkdir(parents=True)
            raw = pd.concat([_frame(time=100 * t, extra=richer)
                             for t in range(8)], ignore_index=True)
            raw.to_parquet(out / "extensions/lipid_contacts/matrix.parquet")
            reps[system].append(out)
    cfg = {"analyses": {"lipid_contacts": {"enabled": True, "classes": CLASSES}},
           "statistics": {"enabled": True, "bootstrap_draws": 64, "seed": 17,
                          "ci": .95,
                          "block": {"length": "auto", "tau_multiplier": 2},
                          "hierarchy": {"resample_replicates": True},
                          "hypothesis_tests": {"design": "independent"}}}
    analyze_lipid_contacts(cfg, tmp_path, reps)
    out = tmp_path / "project_summary/statistics/lipid_contacts"
    replica = pd.read_parquet(out / "replicate_statistics.parquet")
    systems = pd.read_parquet(out / "system_statistics.parquet")
    comp = pd.read_parquet(out / "system_comparisons.parquet")
    assert set(replica.level) == {"class", "species"}
    assert set(comp.leaflet) == {"upper", "lower"}
    assert systems.n_replicates.eq(3).all()
    selected = comp[(comp.level == "class") & (comp.leaflet == "upper") &
                    (comp.species_a == "PC") & (comp.species_b == "PC") &
                    (comp.metric == "pair_count_undirected")].iloc[0]
    assert selected.difference_system2_minus_system1 == 2
    assert selected.difference_ci_low == 2
    assert selected.difference_ci_high == 2
    assert selected.ci_excludes_zero
    lower = comp[(comp.level == "class") & (comp.leaflet == "lower") &
                 (comp.species_a == "PC") & (comp.species_b == "PC") &
                 (comp.metric == "pair_count_undirected")].iloc[0]
    assert lower.difference_system2_minus_system1 == 0
    assert not lower.ci_excludes_zero
    missing = comp[(comp.level == "species") &
                   (comp.species_a == "GM1a") &
                   (comp.species_b == "POPC")]
    assert missing.difference_ci_low.isna().all()
