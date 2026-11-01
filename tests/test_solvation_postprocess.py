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


def test_coordination_aggregation_equal_weight_targets():
    from mdscope.observables.solvation_postprocess import aggregate_coordination
    values = pd.DataFrame({
        "solvent": ["water"]*4,
        "target_class": ["lipid"]*4,
        "site_class": ["GM1a"]*4,
        "target_id": ["A", "B", "A", "B"],
        "time_ps": [0., 0., 5., 5.],
        "n_solvent": [2, 6, 0, 4]})
    by_class, by_target = aggregate_coordination(values)
    assert by_class.mean_n_solvent.tolist() == [4., 2.]
    assert by_class.target_count.tolist() == [2, 2]
    assert by_target.set_index("target_id").loc["A", "hydrated_fraction"] == .5
    assert by_target.set_index("target_id").loc["B", "mean_n_solvent"] == 5.


def test_chunked_coordination_summary_matches_whole_dataframe(tmp_path):
    from mdscope.observables.solvation_postprocess import (
        aggregate_coordination, aggregate_coordination_parts)
    rng = np.random.default_rng(134)
    frames = []
    for step in range(38):
        for target in range(5):
            frames.append((float(step), "water", "protein", "acidic",
                           f"p{target}", int(rng.integers(0, 7))))
    df = pd.DataFrame(frames, columns=[
        "time_ps", "solvent", "target_class", "site_class",
        "target_id", "n_solvent"])
    paths = []
    for i, frame_indices in enumerate(np.array_split(np.arange(38), 6)):
        path = tmp_path / f"counts-{i:06d}.parquet"
        df[df.time_ps.isin(frame_indices)].to_parquet(path, index=False)
        paths.append(path)
    # A frame can span chunks only if a caller explicitly splits it;
    # the fused flush machinery always writes entire frames.
    whole_class, whole_target = aggregate_coordination(df)
    partial_class, partial_target = aggregate_coordination_parts(paths)
    for actual, expected in (
        (partial_class, whole_class), (partial_target, whole_target)):
        keys = [k for k in ("time_ps", "target_id", "solvent")
                if k in actual.columns]
        actual = actual.sort_values(keys).reset_index(drop=True)
        expected = expected.sort_values(keys).reset_index(drop=True)
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False,
                                      rtol=1e-10, atol=1e-10)


def test_chunked_coordination_summary_empty_inputs():
    from mdscope.observables.solvation_postprocess import aggregate_coordination_parts
    by_class, by_target = aggregate_coordination_parts([])
    assert by_class.empty and by_target.empty
