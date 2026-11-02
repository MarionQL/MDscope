"""Independent, interpretable numerical oracles for MDscope physical contracts.

These tests use hand-countable configurations and formulas that do not call
the implementation to compute their expected answers. Passing them does not
establish converged simulation estimates or transferable binding affinities.
"""
import numpy as np
import pandas as pd
import pytest
import MDAnalysis as mda

from mdscope.mechanisms.motif_residence import (
    summarize_lipid_residence, summarize_intermittent_residence,
)
from mdscope.mechanisms.sdc4_mechanisms import (
    _composition_enrichment, _contacts, _accessible_area_grid,
)


def _events(times):
    return pd.DataFrame({
        "time_ps": np.asarray(times, float),
        "motif": ["NxIP"] * len(times),
        "molecule_id": [7] * len(times),
        "species": ["GM1a"] * len(times),
    })


def test_composition_reference_two_species_and_missing_denominator():
    # One GM1a versus three POPC molecules in the sampled leaflet.
    # Both have one contacting molecule. Contact shares are 1/2 each;
    # abundance shares are 1/4 and 3/4 respectively.
    table = pd.DataFrame({
        "time_ps": [0.] * 5,
        "region": ["TM_A"] * 5,
        "leaflet": ["upper"] * 5,
        "species": ["GM1a", "POPC", "GM1a", "BGAL", "BGAL"],
        "moiety": ["any", "any", "glycan", "any", "ceramide_tail"],
        "contacted_molecules": [1, 1, 1, 0, 0],
        "available_molecules": [1, 3, 1, 0, 0],
        "locally_available_molecules": [1, 3, 1, 0, 0],
    })
    for column in ("available_molecules", "locally_available_molecules"):
        result = _composition_enrichment(table, column)
        assert result[0] == pytest.approx(2.0)
        assert result[1] == pytest.approx(2 / 3)
        assert np.isnan(result[2])  # no double-counting moiety subclasses
        assert np.isnan(result[3])  # BGAL unavailable, not zero-enrichment
        assert np.isnan(result[4])


def test_local_population_mismatch_reference_is_detectable():
    # Contacting GM1a is beyond the radius used for local anchors.
    # This is a deliberately diagnostic counterexample: 1 contact / 0
    # locally available is not a physically valid local occupancy.
    table = pd.DataFrame({
        "time_ps": [0., 0.], "region": ["TM_A"] * 2,
        "leaflet": ["upper"] * 2, "species": ["GM1a", "POPC"],
        "moiety": ["any", "any"], "contacted_molecules": [1, 1],
        "available_molecules": [1, 3],
        "locally_available_molecules": [0, 3],
    })
    result = _composition_enrichment(table, "locally_available_molecules")
    assert np.isnan(result[0])
    assert result[1] == pytest.approx(.5)  # contact share 1/2, local share 3/3
    assert table.loc[0, "contacted_molecules"] > table.loc[0, "locally_available_molecules"]
    # The oracle exposes a population mismatch; it does not silently
    # claim local enrichment for the unavailable GM1a.


def test_orthorhombic_minimum_image_contact_oracle():
    # Across a 10 Å periodic boundary: 9.7 -> 0.2 is 0.5 Å,
    # while 9.7 -> 4.2 is 4.5 Å. Test strict cutoff behavior too.
    class Atoms:
        def __init__(self, coords):
            self.positions = np.asarray(coords, dtype=np.float32)
        def __len__(self):
            return len(self.positions)
    box = np.array([10., 10., 10., 90., 90., 90.])
    left = Atoms([[9.7, 5., 5.]])
    right = Atoms([[.2, 5., 5.], [4.2, 5., 5.]])
    actual = _contacts(left, right, 1., box)
    assert actual.tolist() == [[0, 0]]
    assert _contacts(left, right, .4, box).size == 0
    # Rigid translation by one box length must preserve contacts.
    translated = Atoms([[19.7, 5., 5.]])
    np.testing.assert_array_equal(_contacts(translated, right, 1., box), actual)


def test_projected_area_grid_periodicity_and_slab_oracle():
    grid = np.array([[x, y] for x in (-1., 0., 1.)
                     for y in (-1., 0., 1.)], float)
    box = np.array([10., 10., 10., 90., 90., 90.])
    # One protein atom at periodic x=-0.2 blocks only (0,0) on a
    # 1 Å grid when exclusion radius is 0.4 Å.
    wrapped = _accessible_area_grid([0., 0.], 0., [[9.8, 0., 0.]],
                                    box, grid, 1., .4, 1.)
    direct = _accessible_area_grid([0., 0.], 0., [[-.2, 0., 0.]],
                                   box, grid, 1., .4, 1.)
    outside_slab = _accessible_area_grid([0., 0.], 0., [[9.8, 0., 2.]],
                                         box, grid, 1., .4, 1.)
    assert wrapped == pytest.approx(8.)
    assert direct == pytest.approx(8.)
    assert outside_slab == pytest.approx(9.)


def test_residence_stride_and_gap_tolerance_are_observation_dependent():
    # Same hypothetical event: contact at 0, 2 and 4 ps, absent
    # at 1 and 3. A 2-ps stride cannot detect the two missing frames.
    dense_times = np.arange(5., dtype=float)
    sparse_times = np.array([0., 2., 4.])
    events = _events([0., 2., 4.])
    dense = summarize_lipid_residence(events, dense_times)
    sparse = summarize_lipid_residence(events, sparse_times)
    assert dense.n_frames.tolist() == [1, 1, 1]
    assert sparse.n_frames.tolist() == [3]
    assert dense.occupancy_fraction.tolist() == pytest.approx([3 / 5] * 3)
    assert sparse.occupancy_fraction.tolist() == pytest.approx([1.])
    bridged = summarize_intermittent_residence(events, dense_times, gap_frames=1)
    assert bridged.n_frames.tolist() == [3]
    assert bridged.duration_ps.tolist() == pytest.approx([5.])


def test_residence_censoring_and_irregular_last_interval_oracle():
    times = np.array([0., 3., 7., 12.])
    events = _events([0., 3., 12.])
    result = summarize_lipid_residence(events, times)
    assert len(result) == 2
    assert result.left_censored.tolist() == [True, False]
    assert result.right_censored.tolist() == [False, True]
    assert result.n_frames.tolist() == [2, 1]
    # The first dwell includes its next sampled interval (3->7).
    assert result.duration_ps.iloc[0] == pytest.approx(7.)
    # The terminal episode has no observed next interval. Current
    # convention extrapolates the median spacing (3,4,5)->4.
    assert result.duration_ps.iloc[1] == pytest.approx(4.)
    assert result.occupancy_fraction.tolist() == pytest.approx([.75, .75])


def test_ceramide_head_tail_glycan_chemistry_selection_oracle():
    # A synthetic topology has no force-field chemistry inference:
    # label identities and atom names are explicit, checkable inputs.
    u = mda.Universe.empty(
        7, n_residues=3,
        atom_resindex=[0, 0, 0, 0, 1, 1, 2],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["NF", "O3", "C4S", "C2F",
                                   "C1", "O5", "O1"])
    u.add_TopologyAttr("resnames", ["CER160", "BGLCA", "CER160"])
    u.add_TopologyAttr("chainIDs", ["GM1a"] * 6 + ["CER160"])
    head = u.select_atoms(
        "chainID GM1a CER160 and resname CER160 and name NF O3 O1"
    )
    tail = u.select_atoms(
        "chainID GM1a CER160 and resname CER160 and name C4S C2F"
    )
    glycan = u.select_atoms(
        "chainID GM1a and not resname CER160 and not name H*"
    )
    assert set(head.indices.tolist()) == {0, 1, 6}
    assert set(tail.indices.tolist()) == {2, 3}
    assert set(glycan.indices.tolist()) == {4, 5}
    assert not (set(head.indices) & set(tail.indices))
    assert not (set(glycan.indices) & (set(head.indices) | set(tail.indices)))
    assert "O1" in set(head.names)  # free ceramide hydroxyl in control


def test_matched_local_contact_population_independent_reference():
    # Whole leaflet: one GM1 and three POPC, with one contact each.
    # Local disk: GM1 anchor is outside, POPC has 3 available, 1 contacting.
    # Historical local enrichment compares 1 GM1 contact to zero local GM1.
    # Matched local enrichment conditions on the disk, so POPC alone has
    # contact share 1 and availability share 1, yielding enrichment 1.
    frame = pd.DataFrame({
        "time_ps": [0., 0.], "region": ["TM_A"] * 2,
        "leaflet": ["upper"] * 2, "species": ["GM1a", "POPC"],
        "moiety": ["any"] * 2,
        "contacted_molecules": [1, 1],
        "locally_contacted_molecules": [0, 1],
        "available_molecules": [1, 3],
        "locally_available_molecules": [0, 3],
    })
    numerator = frame.loc[:, "locally_contacted_molecules"].to_numpy()
    denominator = frame.loc[:, "locally_available_molecules"].to_numpy()
    assert np.all(numerator <= denominator)
    fractions = np.divide(numerator, denominator,
                          out=np.full(2, np.nan), where=denominator > 0)
    assert np.isnan(fractions[0])
    assert fractions[1] == pytest.approx(1 / 3)
    matched = _composition_enrichment(
        frame, "locally_available_molecules", "locally_contacted_molecules"
    )
    assert np.isnan(matched[0])
    assert matched[1] == pytest.approx(1.)
    historic = _composition_enrichment(frame, "locally_available_molecules")
    assert historic[1] == pytest.approx(.5)


def test_residence_duration_and_occupancy_have_separate_contracts():
    times = np.array([0., 3., 7., 12.])
    result = summarize_lipid_residence(_events([0., 3., 12.]), times)
    assert result.observed_span_ps.tolist() == pytest.approx([3., 0.])
    assert result.episode_contact_fraction.tolist() == pytest.approx([.5, .25])
    assert result.occupancy_fraction.tolist() == pytest.approx([.75, .75])
    assert result.duration_ps.tolist() == pytest.approx([7., 4.])
    assert set(result.duration_convention) == {
        "next_sample_interval_or_terminal_median"
    }


def test_solvation_and_lipid_residence_share_explicit_duration_contract():
    from mdscope.observables.solvation import contact_episodes
    observed = [(0, 0.), (1, 3.), (2, 7.), (3, 12.)]
    events = [(0, "water", "patch", 91),
              (1, "water", "patch", 91),
              (3, "water", "patch", 91)]
    water = contact_episodes(events, observed)
    lipid = summarize_lipid_residence(_events([0., 3., 12.]),
                                      np.array([0., 3., 7., 12.]))
    assert water.duration_ps.tolist() == pytest.approx([3., 0.])
    assert water.sampled_duration_ps.tolist() == pytest.approx([7., 4.])
    assert water.sampled_duration_ps.tolist() == pytest.approx(
        lipid.duration_ps.tolist()
    )
    assert water.duration_convention.tolist() == lipid.duration_convention.tolist()
    with pytest.raises(ValueError, match="outside observation"):
        contact_episodes([(99, "water", "patch", 91)], observed)



@pytest.mark.parametrize("times", [[0., np.nan], [np.inf], [1., 1.], [2., 1.]])
def test_episode_time_axis_rejects_unknown_or_nonincreasing_samples(times):
    from mdscope.core.residence import sampled_episode_durations
    from mdscope.observables.solvation import contact_episodes
    with pytest.raises(ValueError, match="times"):
        sampled_episode_durations(times, 0, 0)
    with pytest.raises(ValueError, match="times"):
        contact_episodes([], list(enumerate(times)))


def test_batched_duration_oracle_and_single_observation():
    from mdscope.core.residence import sampled_episode_durations
    # Intervals: 3, 4, 5 ps; terminal median: 4 ps.
    span, imputed = sampled_episode_durations(
        [0., 3., 7., 12.], np.array([0, 1, 3]), np.array([1, 2, 3]))
    np.testing.assert_allclose(span, [3., 4., 0.])
    np.testing.assert_allclose(imputed, [7., 9., 4.])
    span, imputed = sampled_episode_durations([17.], 0, 0)
    assert span == 0. and np.isnan(imputed)
    for first, last in ((-1, 0), (1, 0), (0, 4), (.5, 1)):
        with pytest.raises(ValueError, match="indices"):
            sampled_episode_durations([0., 3., 7., 12.], first, last)


def test_corrected_enrichment_figures_and_inventory(tmp_path):
    from mdscope.reporting.sdc4_figures import plot_sdc4_mechanisms
    from mdscope.reporting.figure_inventory import figure_inventory
    dest = tmp_path / "carb/rep1/extensions/sdc4_mechanisms"
    dest.mkdir(parents=True)
    pd.DataFrame({
        "region": ["TM_A", "TM_A"], "species": ["GM1a", "GM1a"],
        "leaflet": ["upper", "upper"], "moiety": ["any", "any"],
        "contact_fraction": [.1, .2], "local_composition_enrichment": [2., 4.],
        "matched_local_contact_fraction": [.25, .5],
        "matched_local_composition_enrichment": [1., np.nan],
    }).to_parquet(dest / "lipid_region_specificity.parquet")
    paths = {"carb": [tmp_path / "carb/rep1"]}
    plot_sdc4_mechanisms({"analyses": {"sdc4_mechanisms": {"enabled": True}}}, tmp_path, paths)
    summary = pd.read_csv(dest / "lipid_region_replicate_summary.csv")
    assert summary.mean_matched_local_enrichment.item() == 1.
    assert summary.mean_local_enrichment.item() == 3.
    assert summary.n_matched_enrichment_frames.item() == 1
    inventory = figure_inventory(tmp_path, paths).set_index("figure")
    prefix = "carb/rep1/extensions/sdc4_mechanisms/lipid_region_"
    assert inventory.loc[prefix + "mean_matched_local_enrichment.png", "kind"] == "primary"
    assert inventory.loc[prefix + "mean_local_enrichment.png", "kind"] == "supplement"


def test_rdf_excluded_volume_and_cell_density_reference_oracle():
    from mdscope.observables.solvation import rdf_from_counts
    # A hemisphere inaccessible to solvent halves pair counts even when
    # solvent in the accessible half has exactly the reference density.
    # If aqueous volume is half the cell, full-cell density is also half
    # the aqueous density: the two effects cancel, falsely suggesting g=1.
    shells = np.array([4*np.pi/3, 28*np.pi/3])
    rho_aqueous = .04
    hemisphere_counts = .5 * shells * rho_aqueous
    np.testing.assert_allclose(rdf_from_counts(hemisphere_counts, [0, 1, 2], rho_aqueous), .5)
    np.testing.assert_allclose(rdf_from_counts(hemisphere_counts, [0, 1, 2], rho_aqueous/2), 1.)


def test_diffusion_unit_conversion_independent_linear_curve():
    from mdscope.observables.multi_origin_msd import fit_multi_origin_curve
    # D=.5 nm²/ns = .05 Å²/ps; 2D lateral MSD = 4Dt = .2t Å².
    lag_ps = np.arange(21)*10.
    fit = fit_multi_origin_curve(lag_ps, .2*lag_ps)
    assert fit["diffusion_nm2_ns"] == pytest.approx(.5)
    assert fit["slope_A2_ps"] == pytest.approx(.2)


def test_real_freesasa_sphere_area_and_spherical_cap_occlusion():
    pytest.importorskip("freesasa")
    from mdscope.observables.sasa import sasa_instances, sasa_residue_rows
    u = mda.Universe.empty(2, n_residues=2, atom_resindex=[0, 1], trajectory=True)
    u.add_TopologyAttr("names", ["CA", "CA"])
    u.add_TopologyAttr("elements", ["C", "C"])
    u.add_TopologyAttr("resnames", ["ALA", "ALA"])
    u.add_TopologyAttr("resids", [1, 2])
    u.atoms.positions = [[0., 0., 0.], [30., 0., 0.]]
    cfg = {"selections": {"protein": "all"}, "analyses": {"sasa": {
        "enabled": True, "target_selection": "index 0", "occluder_selection": "all",
        "algorithm": "shrake_rupley", "n_points": 10000, "probe_radius_A": 1.4,
    }}}
    state = sasa_instances(u, cfg)["default"]
    sphere_nm2 = 4*np.pi*(1.7+1.4)**2/100.
    assert sasa_residue_rows(state, 0.)[0][-1] == pytest.approx(sphere_nm2, rel=.01)
    u.atoms[1].position = [3.1, 0., 0.]
    # Equal inflated spheres at separation R: cap height R/2; exposed area
    # of the target = 4pi R² - 2pi R(R/2) = 3pi R².
    assert sasa_residue_rows(state, 0.)[0][-1] == pytest.approx(.75*sphere_nm2, rel=.02)
    # Also exercise the native-crash case of isolated + connected components.
    v = mda.Universe.empty(3, n_residues=3, atom_resindex=[0, 1, 2], trajectory=True)
    for key, values in (("names", ["CA"]*3), ("elements", ["C"]*3),
                        ("resnames", ["ALA"]*3), ("resids", [1, 2, 3])):
        v.add_TopologyAttr(key, values)
    v.atoms.positions = [[0., 0., 0.], [3.1, 0., 0.], [30., 0., 0.]]
    cfg["analyses"]["sasa"]["target_selection"] = "all"
    mixed = sasa_instances(v, cfg)["default"]
    np.testing.assert_allclose([row[-1] for row in sasa_residue_rows(mixed, 0.)],
                               [.75*sphere_nm2, .75*sphere_nm2, sphere_nm2], rtol=.02)


def test_real_sasa_periodic_occluder_and_target_only_reduction():
    pytest.importorskip('freesasa')
    from mdscope.observables.sasa import sasa_instances, sasa_residue_rows
    u = mda.Universe.empty(2, n_residues=2, atom_resindex=[0,1], trajectory=True)
    for key, values in [('names',['CA','C1']),('elements',['C','C']),('resnames',['ALA','LIP']),('resids',[7,99])]:
        u.add_TopologyAttr(key, values)
    u.dimensions = [20,20,20,90,90,90]
    u.atoms.positions = [[1,0,0],[19,0,0]]
    cfg = {'selections': {'protein':'index 0'}, 'analyses': {'sasa': {'enabled':True, 'target_selection':'index 0', 'occluder_selection':'all','periodic':True,'algorithm':'shrake_rupley','n_points':10000}}}
    state = sasa_instances(u,cfg)['default']
    rows = sasa_residue_rows(state,0.)
    assert len(rows)==1 and rows[0][1]==0
    radius=3.1
    expected=4*np.pi*radius**2*(.5+2/(4*radius))/100
    assert rows[0][-1]==pytest.approx(expected,rel=.02)
    u.atoms[0].position += [20,0,0]
    assert sasa_residue_rows(state,0.)[0][-1]==pytest.approx(expected,rel=.02)


def test_sasa_burial_subtraction_and_exact_reference_identity(tmp_path):
    from mdscope.observables.sasa import write_sasa_burial
    root=tmp_path/'profiles/sasa';(root/'environment').mkdir(parents=True)
    baseline=pd.DataFrame({'time_ps':[0.,10.],'topology_residue_index':[7,7],
                           'residue':['ALA7','ALA7'],'sasa_nm2':[1.2,1.]})
    environment=baseline.assign(sasa_nm2=[.3,.25])
    baseline.to_parquet(root/'residue_sasa.parquet')
    environment.to_parquet(root/'environment/residue_sasa.parquet')
    write_sasa_burial(tmp_path,'environment','default')
    result=pd.read_parquet(root/'environment/residue_burial.parquet')
    np.testing.assert_allclose(result.buried_sasa_nm2,[.9,.75])
    np.testing.assert_allclose(result.buried_fraction,[.75,.75])
    environment.assign(time_ps=[0.,20.]).to_parquet(root/'environment/residue_sasa.parquet')
    with pytest.raises(ValueError,match='frames/residues differ'):
        write_sasa_burial(tmp_path,'environment','default')
    environment.assign(residue=['ALA7','GLY7']).to_parquet(root/'environment/residue_sasa.parquet')
    with pytest.raises(ValueError,match='identity mismatch'):
        write_sasa_burial(tmp_path,'environment','default')
