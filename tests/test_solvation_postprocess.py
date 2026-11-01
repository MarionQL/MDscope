"""Independent checks for cutoff suggestions, SASA reuse and survival."""
import numpy as np
import pandas as pd
from mdscope.observables.solvation_postprocess import (
    suggest_shell_boundaries, kaplan_meier, join_existing_sasa)


def test_cutoff_suggestions_are_explicit_not_assumed():
    edges = np.arange(0, 10, .5)
    centers = edges[:-1] + .25
    g = np.ones(len(centers))
    g[4:11] = [1.2, 2.4, 3.1, 2.2, 1.2, .7, 1.1]
    rdf = pd.DataFrame({"solvent":"water","target_class":"protein",
                         "site_class":"carboxylate", "r_inner_A":edges[:-1],
                         "r_outer_A":edges[1:], "g_bulk":g})
    result = suggest_shell_boundaries(rdf, smoothing_bins=1)
    assert result.status.item() == "candidate_only"
    assert abs(result.suggested_cutoff_A.item() - 4.75) < 1e-9
    flat = rdf.assign(g_bulk=1.)
    assert suggest_shell_boundaries(flat).suggested_cutoff_A.isna().all()


def test_product_limit_censors_and_known_origin():
    episodes = pd.DataFrame([
        ("water", "site", 1., False, False),
        ("water", "site", 2., False, True),
        ("water", "site", 3., False, False),
        ("water", "site", 99., True, False),
    ], columns=["solvent","target_id","duration_ps","left_censored","right_censored"])
    km = kaplan_meier(episodes)
    np.testing.assert_allclose(km.survival_probability, [2/3, 2/3, 0])
    assert km.at_risk.tolist() == [3,2,1]
    assert km.right_censored.tolist() == [0,1,0]


def test_sasa_exact_frame_join_and_missing_residue():
    counts = pd.DataFrame({"time_ps":[0.,5.],"target_id":["patch","patch"],
                           "n_solvent":[2,4]})
    mapping = pd.DataFrame({"target_id":["patch","patch"],
                            "topology_residue_index":[1,2]})
    sasa = pd.DataFrame({"time_ps":[0.,0.,5.],
                         "topology_residue_index":[1,2,1],
                         "sasa_nm2":[.2,.3,.8]})
    result = join_existing_sasa(counts,mapping,sasa)
    assert result.sasa_nm2.iloc[0] == .5
    assert np.isnan(result.sasa_nm2.iloc[1])
    assert result.sasa_coverage.tolist() == [1.,.5]
