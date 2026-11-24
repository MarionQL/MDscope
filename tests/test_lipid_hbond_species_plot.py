import numpy as np
import pandas as pd

from mdscope.reporting.plotting import _lipid_hbond_pair_matrix


def test_lipid_hbond_pair_matrix():
    def table(pairs, occupancies):
        return pd.DataFrame({"interaction_class": pairs, "occupancy": occupancies})

    first = table(["GM1a__PSM", "GM1a__PSM", "PSM__PSM"], [0.8, 0.7, 0.2])
    second = table(["GM1a__GM1a", "PSM__PSM"], [0.4, 0.3])
    third = table(["CER160__PSM"], [0.5])
    pairs, matrix = _lipid_hbond_pair_matrix([first, second, third, None])
    assert pairs == ["CER160__PSM", "GM1a__GM1a", "GM1a__PSM", "PSM__PSM"]
    np.testing.assert_allclose(matrix, [
        [np.nan, np.nan, 0.5, np.nan],
        [0, 0.4, np.nan, np.nan],
        [1.5, 0, np.nan, np.nan],
        [0.2, 0.3, 0, np.nan],
    ], equal_nan=True)
