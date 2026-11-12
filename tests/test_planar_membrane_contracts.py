"""Independent planar geometry facts shared across membrane observers."""
import numpy as np
import pandas as pd
import pytest
from MDAnalysis.lib.mdamath import triclinic_vectors

from mdscope.core.lipid_geometry import anchor_groups, aligned_leaflets
from mdscope.pipeline.fused_membrane import FusedMembrane
from mdscope.observables.solvation import SolvationObserver
from mdscope.observables.bilayer_normal import fit_local_normals
from test_fused_membrane import _universe, _cfg


@pytest.mark.parametrize('axis,angles', [
    ('x', [90., 75., 90.]), ('y', [90., 90., 75.]), ('z', [75., 90., 90.]),
])
def test_leaflet_geometry_rejects_lattice_shear_through_normal(axis, angles):
    u = _universe()
    u.dimensions = [100., 100., 100., *angles]
    with pytest.raises(ValueError, match='normal'):
        anchor_groups(u, dict(axis=axis, species_selections={'POPC':'name P'}))


@pytest.mark.parametrize('family', ['bilayer_normal', 'lipid_enrichment', 'lipid_contacts'])
def test_box_becoming_incompatible_fails_on_current_frame(tmp_path, family):
    u, cfg = _universe(), _cfg()
    u.add_TopologyAttr('bonds', [])
    cfg['analyses'] = {family: cfg['analyses'].get(family, dict(
        enabled=True, selection='resname POPC', species_selections={'POPC':'name P'}))}
    observer = FusedMembrane(u, cfg, tmp_path)
    observer.update(u.trajectory.ts, 0)
    u.dimensions = [100., 100., 100., 75., 90., 90.]
    with pytest.raises(ValueError, match='normal'):
        observer.update(u.trajectory.ts, 1)
    assert observer.frame_count == 1


@pytest.mark.parametrize('gamma', [90., 70.])
def test_shared_observers_inplane_skew_and_periodic_crossing(tmp_path, gamma):
    roots = []
    for wrap in (False, True):
        u, cfg = _universe(), _cfg()
        u.add_TopologyAttr('bonds', [])
        if gamma != 90.:
            cfg['analyses']['lipid_enrichment']['enabled'] = False
        cfg['analyses']['lipid_contacts'] = dict(enabled=True, selection='resname POPC',
            species_selections={'POPC':'name P'}, cutoff_A=16.)
        cfg['analyses']['solvation'] = dict(enabled=True, water_selection='resname POPC',
            targets={'protein': {'CA':'name CA'}}, search_radius_A=25., shell_cutoffs_A={'water':25.})
        boxes = u.trajectory.dimensions_array
        boxes[:, 5] = gamma
        lattice = triclinic_vectors(boxes[0]).astype(float)
        if wrap:
            # Drift crosses the periodic normal boundary. Independent lattice
            # shifts must preserve distances, physical leaves and contact graph.
            xyz = u.trajectory.coordinate_array
            xyz += [47., -23., 47.]
            xyz += (((np.arange(30).reshape(10,3) % 3)-1) @ lattice)[None]
        u.trajectory[0]
        root = tmp_path/str(wrap)
        observer = FusedMembrane(u, cfg, root)
        solvent = SolvationObserver(u, cfg, root)
        before = u.trajectory.coordinate_array.copy()
        for i, ts in enumerate(u.trajectory):
            observer.update(ts, i)
            solvent.update(ts, i, float(i))
        observer.flush(0); solvent.flush(0)
        observer.finalize(); solvent.finalize()
        np.testing.assert_array_equal(u.trajectory.coordinate_array, before)
        roots.append(root)
    paths = ['bilayer_normal/timeseries.parquet', 'lipid_contacts/matrix.parquet',
             'solvation/counts.parquet', 'solvation/rdf.parquet']
    if gamma == 90.:
        paths.append('lipid_footprint/counts.parquet')
    for relative in paths:
        left, right = [pd.read_parquet(root/'extensions'/relative) for root in roots]
        pd.testing.assert_frame_equal(left, right, check_exact=False, atol=1e-4, rtol=1e-5)
    matrix = pd.read_parquet(roots[0]/'extensions/lipid_contacts/matrix.parquet')
    assert matrix.pair_count_undirected.eq(4).all()  # square edges, not diagonals
    assert matrix.molecules_a.eq(4).all()
    normals = pd.read_parquet(roots[0]/'extensions/bilayer_normal/timeseries.parquet')
    assert normals.bilayer_normal_upper_mean_deg.eq(0).all()
    assert normals.bilayer_normal_valid_fraction.eq(1).all()


def test_leaflet_exchange_and_label_reversal_have_same_physical_membership():
    initial = np.array([20.]*4 + [60.]*4)
    reference = np.array([False]*4 + [True]*4)
    crossed = initial.copy()
    crossed[[0, 4]] = crossed[[4, 0]]
    expected = reference.copy()
    expected[[0, 4]] = ~expected[[0, 4]]
    for labels in (reference, ~reference):
        current, changes = aligned_leaflets(crossed + 47., 100., labels)
        np.testing.assert_array_equal(current, expected if labels[0] == reference[0] else ~expected)
        assert changes == 2


def test_tilted_plane_analytic_angle_and_sparse_unavailable():
    x, y = np.meshgrid(np.arange(4)*5., np.arange(4)*5.)
    plane = np.column_stack([x.ravel(), y.ravel(), (.12*x + .05*y + 20.).ravel()])
    box = [100., 100., 100., 90., 90., 70.]
    labels = np.zeros(len(plane), bool)
    angles, _ = fit_local_normals(plane, labels, box, radius=24., min_neighbors=3)
    np.testing.assert_allclose(angles, np.degrees(np.arctan(.13)), atol=1e-5)
    angles, neighbors = fit_local_normals(plane, labels, box, radius=1., min_neighbors=3)
    assert np.isnan(angles).all()
    assert not neighbors.any()


@pytest.mark.parametrize('shape', ['tilted', 'curved'])
@pytest.mark.parametrize('wrap', [False, True])
def test_normals_and_local_surfaces_share_known_height_fields(tmp_path, shape, wrap):
    from test_sdc4_geometry_references import _fixture
    from mdscope.mechanisms.sdc4_mechanisms import SDC4Mechanisms

    u, cfg = _fixture()
    anchors = u.select_atoms('name P')
    xyz = anchors.positions.copy()
    x, y = xyz[:, 0], xyz[:, 1]
    # Binary-exact plane heights remain exact in float32 after image shifts.
    xyz[:, 2] += .125*x + .0625*y if shape == 'tilted' else .002*x*x + .003*y*y
    anchors.positions = xyz
    if wrap:
        xyz = u.atoms.positions.astype(float) + [47., -23., 47.]
        xyz += ((np.arange(len(xyz)*3).reshape(-1, 3) % 3)-1)*100.
        u.atoms.positions = xyz
    cfg.update(time=_cfg()['time'], trajectory=_cfg()['trajectory'])
    cfg['analyses']['bilayer_normal'] = dict(enabled=True, species_selections={'POPC':'name P'},
        neighbor_radius_A=35., min_neighbors=3, save_lipid_angles=True)
    normals = FusedMembrane(u, cfg, tmp_path/'normals')
    mechanism = SDC4Mechanisms(u, cfg, tmp_path/'mechanism')
    before = u.atoms.positions.copy()
    normals.update(u.trajectory.ts, 0)
    mechanism.update(u.trajectory.ts, 0)
    row = mechanism.rows['timeseries'][0]
    assert row['sdc4_local_thickness_A'] == pytest.approx(32., abs=1e-5)
    for leaf in ('upper', 'lower'):
        assert row[f'sdc4_local_{leaf}_curvature_Ainv'] == pytest.approx(
            0. if shape == 'tilted' else .005, abs=1e-7)
    angles = pd.DataFrame(normals.rows['angles'])
    if shape == 'tilted':
        np.testing.assert_allclose(angles.normal_deviation_deg,
                                   np.degrees(np.arctan(np.sqrt(5)/16)), atol=1e-5)
    else:
        # At x=y=0 a symmetric quadratic has zero slope in either leaf.
        assert angles.set_index('resindex').loc[[28, 37], 'normal_deviation_deg'].max() < 1e-5
    np.testing.assert_array_equal(u.atoms.positions, before)


@pytest.mark.parametrize('gamma', [90., 70.])
def test_msd_keeps_initial_cohorts_under_independent_lattice_wraps(tmp_path, gamma):
    from mdscope.observables.membrane import initialize_membrane
    from mdscope.observables.msd import msd_metrics

    u, cfg = _universe(), _cfg()
    u.add_TopologyAttr('masses', np.ones(len(u.atoms)))
    u.trajectory.dimensions_array[:, 5] = gamma
    lattice = triclinic_vectors(u.trajectory.dimensions_array[0]).astype(float)
    # Initial cohorts are deliberately established on a resolved, centered
    # first frame. Subsequent image choices must not change their MSD.
    u.trajectory.coordinate_array[1:] += (((np.arange(30).reshape(10,3) % 3)-1) @ lattice)[None]
    u.trajectory[0]
    cfg['analyses']['lipid_enrichment']['enabled'] = False
    cfg['analyses']['membrane'] = dict(enabled=True, headgroup_selection='name P',
        area_per_lipid={'enabled':False}, thickness={'enabled':False}, msd={'enabled':True})
    state = initialize_membrane(u, cfg)
    observer = FusedMembrane(u, cfg, tmp_path)
    for i, ts in enumerate(u.trajectory):
        observer.update(ts, i)
        metrics = msd_metrics(u, state)
        for metric in ('membrane_msd_lateral_A2', 'membrane_msd_lateral_upper_A2',
                       'membrane_msd_lateral_lower_A2'):
            assert metrics[metric] == pytest.approx((.1*i)**2, abs=2e-5)
        assert observer.rows['bilayer_normal'][-1]['bilayer_normal_valid_fraction'] == 1.


@pytest.mark.parametrize('normal', [0, 1, 2])
def test_inplane_skew_preserves_known_leaflets_for_each_axis(normal):
    u = _universe()
    xyz = u.atoms.positions.copy()
    xyz[:, [normal, 2]] = xyz[:, [2, normal]]
    u.atoms.positions = xyz
    angles = [90., 90., 90.]
    angles[normal] = 75.  # opposite angle: shear confined to lateral plane
    u.dimensions = [100., 100., 100., *angles]
    _, _, _, upper = anchor_groups(u, dict(axis='xyz'[normal],
        species_selections={'POPC':'name P'}))
    np.testing.assert_array_equal(upper, [False]*4 + [True]*4)
