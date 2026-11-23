"""Raw residue associations: exact joins, availability and no bootstrap."""
import numpy as np
import pandas as pd
import pytest

from mdscope.statistics.residue_correlations import correlate_replica
from mdscope.config.validators.statistics import residue_correlations_settings


def _write(path, frame):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)


def test_raw_correlations_join_by_exact_time_and_topology_identity(tmp_path):
    times = np.arange(20, dtype=float) * 10.
    idx = [3] * 20 + [9] * 20
    labels = ["PROA:LEU3"] * 20 + ["PROB:LEU9"] * 20
    water = np.r_[np.arange(20), np.arange(20)[::-1]]
    base = pd.DataFrame({"time_ps": np.tile(times, 2),
                         "topology_residue_index": idx, "residue": labels,
                         "water_count": water,
                         "protein_fraction": water / 20.,
                         "hs_fraction": np.zeros(40),
                         "lipid_fraction": np.zeros(40)})
    _write(tmp_path / "extensions/sdc4_mechanisms/residue_occlusion.parquet", base)
    # A deliberately shuffled source with the same times must be joined by time,
    # not by row number; the one missing frame must stay missing.
    sasa = base.loc[(base.time_ps != 30) | (base.topology_residue_index != 3),
                    ["time_ps", "topology_residue_index", "residue"]].copy()
    sasa["sasa_nm2"] = base.loc[sasa.index, "water_count"].to_numpy() * 2.
    _write(tmp_path / "profiles/sasa/isolated/residue_sasa.parquet",
           sasa.sample(frac=1, random_state=3))
    result = correlate_replica(tmp_path, [
        ("water_count", "sasa_isolated"), ("water_count", "hs_fraction")
    ], min_frames=5)
    first = result.loc[(result.topology_residue_index == 3)
                       & (result.metric_y == "sasa_isolated")].iloc[0]
    assert first.n_paired_frames == 19
    assert first.pearson_r == pytest.approx(1.0)
    assert first.spearman_rho == pytest.approx(1.0)
    flat = result.loc[result.metric_y == "hs_fraction"]
    assert set(flat.status) == {"constant_series"}
    assert flat.pearson_r.isna().all()


def test_missing_metric_does_not_turn_into_zero(tmp_path):
    frame = pd.DataFrame({"time_ps": [0., 10., 20.],
                          "topology_residue_index": [2, 2, 2],
                          "residue": ["PROA:SER2"] * 3,
                          "water_count": [1, 2, 3]})
    _write(tmp_path / "extensions/sdc4_mechanisms/residue_occlusion.parquet", frame)
    result = correlate_replica(tmp_path, [("water_count", "sasa_environment")],
                               min_frames=3)
    assert len(result) == 1
    assert result.status.iloc[0] == "missing_metric"
    assert np.isnan(result.pearson_r.iloc[0])


def test_duplicate_frames_or_conflicting_labels_fail(tmp_path):
    src = tmp_path / "extensions/sdc4_mechanisms/residue_occlusion.parquet"
    frame = pd.DataFrame({"time_ps": [0., 0.],
                          "topology_residue_index": [2, 2],
                          "residue": ["PROA:SER2"] * 2,
                          "water_count": [1., 2.]})
    _write(src, frame)
    with pytest.raises(ValueError, match="duplicate time/residue"):
        correlate_replica(tmp_path, [("water_count", "sasa_isolated")])


def test_correlation_config_rejects_unsupported_pairs():
    cfg = {"statistics": {"residue_correlations": {
        "enabled": True, "pairs": [["water_count", "unknown"]], "min_frames": 20
    }}}
    with pytest.raises(ValueError, match="supported metric pairs"):
        residue_correlations_settings(cfg)
    cfg["statistics"]["residue_correlations"]["pairs"] = [["water_count", "protein_fraction"]]
    assert residue_correlations_settings(cfg)["enabled"] is True


def test_raw_water_correlations_exclude_residues_without_polar_sites(tmp_path):
    frame = pd.DataFrame({
        "time_ps": [0., 10., 20.] * 2,
        "topology_residue_index": [2] * 3 + [3] * 3,
        "residue": ["PROA:GLY2"] * 3 + ["PROA:SER3"] * 3,
        "water_count": [0, 0, 0, 0, 1, 2],
        "polar_atom_count": [0] * 3 + [1] * 3,
        "protein_fraction": [0., .5, 1.] * 2,
    })
    _write(tmp_path / "extensions/sdc4_mechanisms/residue_occlusion.parquet", frame)
    result = correlate_replica(
        tmp_path, [("water_count", "protein_fraction")], min_frames=3,
    ).set_index("topology_residue_index")
    assert result.loc[2, "n_paired_frames"] == 0
    assert result.loc[2, "status"] == "insufficient_paired_frames"
    assert result.loc[3, "n_paired_frames"] == 3
    assert result.loc[3, "pearson_r"] == pytest.approx(1.)
