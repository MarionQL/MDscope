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


def test_chunk_shell_calibration_reports_stable_vs_ambiguous():
    from mdscope.observables.solvation_postprocess import calibration_stability
    edges = np.arange(0, 10, .5)
    # Construct PAIR COUNTS from a known RDF, including r² shell volume.
    # Flat counts are not a flat RDF and would spuriously remove its peak.
    g = np.ones(len(edges)-1)
    g[4:14] = [1.2, 2.4, 3.1, 2.2, 1.2, .7, .9, 1., 1., 1.]
    volume = 4*np.pi/3*np.diff(edges**3)
    exposure = 1000.
    counts = np.rint(g*volume*exposure).astype(int)
    parts = [{"solvent": "water", "target_class": "protein",
              "site_class": "carboxylate", "exposure": exposure,
              **{f"bin_{j}": int(x) for j, x in enumerate(counts)}}
             for _ in range(4)]
    stable = calibration_stability(pd.DataFrame(parts), edges)
    assert stable.status.item() == "stable"
    assert stable.chunks_resolved.item() == 4
    parts[-1].update({f"bin_{j}": int(x) for j, x in
                      enumerate(np.rint(volume*exposure).astype(int))})
    unstable = calibration_stability(pd.DataFrame(parts), edges)
    assert unstable.status.item() == "unresolved_chunks"
    assert unstable.chunks_resolved.item() == 3
