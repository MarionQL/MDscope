from types import SimpleNamespace

import numpy as np
import pandas as pd

from mdscope.mechanisms.sdc4_mechanisms import _add_acidic_patch_sasa
from mdscope.observables.multi_origin_pipeline import _sampled_window


def test_msd_window_preserves_sampled_time_index():
    trajectory = [SimpleNamespace(time=frame * 4) for frame in range(24)]
    cfg = {"trajectory": {"start": 3, "stop": 20, "stride": 2}}
    indices, offset = _sampled_window(trajectory, cfg, [19, 54])
    assert list(indices) == [5, 7, 9, 11, 13]
    assert offset == 1
    # Configured times use the sampled index rather than the raw frame index.
    cfg["time"] = {"dt_ps": 8, "origin_ps": 12}
    assert _sampled_window(trajectory, cfg, [19, 54]) == (indices, offset)


def test_acidic_patch_sasa_uses_matched_atom_and_residue_totals(tmp_path):
    isolated, environment = tmp_path / "isolated", tmp_path / "environment"
    isolated.mkdir()
    environment.mkdir()
    pd.DataFrame({
        "time_ps": [0.0, 0.0], "topology_residue_index": [10, 11],
        "environment_sasa_nm2": [0.5, 1.5], "buried_sasa_nm2": [0.3, 0.2],
    }).to_parquet(environment / "residue_burial.parquet")
    for folder, area in ((isolated, 0.8), (environment, 0.5)):
        pd.DataFrame({
            "time_ps": [0.0, 0.0], "atom_index": [100, 101],
            "topology_residue_index": [10, 10], "atom_name": ["OD1", "CA"],
            "sasa_nm2": [area, 5.0],
        }).to_parquet(folder / "atom_sasa.parquet")
    result = _add_acidic_patch_sasa(
        pd.DataFrame({"time_ps": [0.0]}), {"acidic_patch_A": {"resindices": [10, 11]}},
        isolated, environment,
    )
    prefix = "sdc4_motif_acidic_patch_A_"
    np.testing.assert_allclose(result[prefix + "environment_sasa_nm2"], 2)
    np.testing.assert_allclose(result[prefix + "buried_sasa_nm2"], 0.5)
    np.testing.assert_allclose(result[prefix + "acidic_oxygen_buried_sasa_nm2"], 0.3)
