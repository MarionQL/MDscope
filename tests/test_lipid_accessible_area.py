"""Independent disk/plane references for the optional lipid area probe model."""
import json
import numpy as np
import pandas as pd
import pytest
from mdscope.observables.lipid_enrichment import accessible_plane_areas, enrichment_frame
from test_fused_membrane import _universe, _cfg, _execute

BOX = np.array([40., 40., 60., 90., 90., 90.])
CENTER = np.array([20., 20., 30.])


def test_empty_plane_and_circle_area_against_analytic_reference():
    full, shell = accessible_plane_areas([], CENTER, 30., BOX, 0., 10., 2, 6., 5.)
    np.testing.assert_allclose([full,shell], [1600., 100*np.pi])
    # The retained model removes a radius-6 disk for protein inside the slab.
    # The annular quadrature is uniform in r^2, not r.
    expected = [1600.-36*np.pi, 64*np.pi]
    for n in (64,128):
        actual = accessible_plane_areas([CENTER], CENTER, 30., BOX, 0., 10., 2, 6., 5., n)
        np.testing.assert_allclose(actual, expected, atol=4., rtol=0)
    # Within the slab, dz=4 retains that disk; no spherical cross-section is implied.
    actual = accessible_plane_areas([CENTER], CENTER, 34., BOX, 0., 10., 2, 6., 5.,128)
    np.testing.assert_allclose(actual, [1600.-36*np.pi, 64*np.pi], atol=4., rtol=0)
    # An atom outside the slab does not occlude the other leaflet.
    np.testing.assert_allclose(accessible_plane_areas([CENTER], CENTER, 40., BOX,
        0.,10.,2,6.,5.), [1600.,100*np.pi])


def test_union_not_sum_and_pbc_translation_invariance():
    first = accessible_plane_areas([CENTER], CENTER, 30., BOX, 0., 10., 2, 6., 5.)
    duplicate = accessible_plane_areas([CENTER,CENTER], CENTER, 30., BOX, 0., 10., 2, 6., 5.)
    np.testing.assert_array_equal(first, duplicate)
    offset = np.array([33.,-29.,41.])
    wrapped = accessible_plane_areas([CENTER+offset+BOX[:3]], CENTER+offset,
        30.+offset[2], BOX, 0., 10., 2, 6., 5.)
    np.testing.assert_allclose(wrapped, first, atol=1e-10)


@pytest.mark.parametrize('box,outer,radius,n', [
    ([40.,40.,60.,90.,90.,80.],10.,6.,64),
    (BOX,21.,6.,64), (BOX,10.,0.,64), (BOX,10.,6.,15), (BOX,10.,6.,64.5),
])
def test_unsupported_geometry_and_parameters_fail(box,outer,radius,n):
    with pytest.raises(ValueError):
        accessible_plane_areas([CENTER],CENTER,30.,box,0.,outer,2,radius,5.,n)


def test_species_absence_and_density_definition_keep_composition_distinct():
    species = np.array(['A','B','A','B','A','B','A','B'])
    upper = np.array([True]*4+[False]*4)
    xyz = np.array([[22,20,30],[23,20,30],[35,20,30],[36,20,30],
                    [22,20,10],[23,20,10],[35,20,10],[36,20,10]],float)
    old = enrichment_frame(species,upper,xyz,CENTER,BOX,0.,10.,2)
    new = enrichment_frame(species,upper,xyz,CENTER,BOX,0.,10.,2,
        categories=['A','B','absent'],obstacles=np.empty((0,3)),
        area_reference={'exclusion_radius_A':2.,'slab_half_width_A':5.})
    for row in old:
        actual = next(r for r in new if r['leaflet']==row['leaflet'] and r['species']==row['species'])
        for key,value in row.items():
            assert actual[key] == value
        assert actual['enrichment_fold'] == 1.
        assert actual['area_density_ratio'] == pytest.approx(.5*1600/(100*np.pi))
    absent = [r for r in new if r['species']=='absent']
    assert all(np.isnan(r['enrichment_fold']) and np.isnan(r['area_density_ratio']) for r in absent)
    assert all(r['area_reference_status']=='unavailable_absent_species' for r in absent)


def test_fused_area_outputs_survive_replay_and_preserve_default_columns(tmp_path):
    default = tmp_path/'default'
    _execute(_universe(),_cfg(),default)
    cfg = _cfg()
    cfg['analyses']['lipid_enrichment']['area_reference'] = {
        'exclusion_radius_A':3.,'slab_half_width_A':5.,'quadrature_order':32}
    cfg['analyses']['lipid_enrichment']['species_selections'] = dict(cfg['analyses']['lipid_enrichment']['species_selections'])
    cfg['analyses']['lipid_enrichment']['allow_absent_species'] = True
    cfg['analyses']['lipid_enrichment']['species_selections']['absent'] = 'resname NONE'
    for label,replay in [('clean',False),('resumed',True)]:
        root=tmp_path/label
        _execute(_universe(),cfg,root,resume=replay)
    relative='extensions/lipid_enrichment'
    clean = pd.read_parquet(tmp_path/'clean'/relative/'counts.parquet')
    replay = pd.read_parquet(tmp_path/'resumed'/relative/'counts.parquet')
    pd.testing.assert_frame_equal(clean,replay)
    original = pd.read_parquet(default/relative/'counts.parquet')
    pd.testing.assert_frame_equal(original,clean.loc[clean.species.ne('absent'),original.columns].reset_index(drop=True))
    absent=clean.loc[clean.species.eq('absent')]
    assert absent.count_bulk.eq(0).all() and absent.area_density_ratio.isna().all()
    metadata=json.loads((tmp_path/'clean'/relative/'area_reference.json').read_text())
    assert metadata['exclusion_radius_A']==3. and metadata['quadrature_order']==32


def test_opposite_periodic_edges_occlude_both_samples():
    from mdscope.core.lipid_geometry import accessible_area_grid
    # Nearest distances are 0.4 and 0.0 A, both within the 0.5-A disk.
    # A nonperiodic KD-tree sees 9.6 A for the first point and misses it.
    assert accessible_area_grid([0.,0.],0.,[[-1e-16,0.,0.]],
        [10.,10.,10.,90.,90.,90.],[[0.,0.]],1.,.5,1.) == 0.
    assert accessible_area_grid([0.,0.],0.,[[4.8,0.,0.]],
        [10.,10.,10.,90.,90.,90.],np.array([[-4.8,0.],[4.8,0.]]),1.,.5,1.) == 0.


def test_uniform_accessible_anchor_population_gives_density_ratio_one():
    # Uniform square midpoint lattice, with exactly the modeled disk removed.
    x,y=np.meshgrid((np.arange(64)+.5)/64*40.,(np.arange(64)+.5)/64*40.)
    xy=np.column_stack([x.ravel(),y.ravel()])
    xy=xy[np.linalg.norm(xy-CENTER[:2],axis=1)>6.]
    n=len(xy)
    positions=np.vstack([np.column_stack([xy,np.full(n,30.)]),np.column_stack([xy,np.full(n,10.)])])
    species=np.tile(np.where(np.arange(n)%2,'A','B'),2)
    upper=np.r_[np.ones(n,bool),np.zeros(n,bool)]
    rows=enrichment_frame(species,upper,positions,CENTER,BOX,0.,10.,2,
        obstacles=np.array([CENTER,[20.,20.,10.]]),
        area_reference=dict(exclusion_radius_A=6.,slab_half_width_A=5.,quadrature_order=128))
    np.testing.assert_allclose([r['area_density_ratio'] for r in rows],1.,atol=.025,rtol=0)


def test_axis_permutation_preserves_planar_area():
    expected=accessible_plane_areas([CENTER],CENTER,30.,BOX,0.,10.,2,6.,5.)
    swapped=np.array([2,1,0])
    np.testing.assert_allclose(accessible_plane_areas([CENTER[swapped]],CENTER[swapped],30.,
        np.r_[BOX[:3][swapped],BOX[3:]],0.,10.,0,6.,5.),expected)


def test_offline_composition_inventory_distinguishes_absent_species(tmp_path):
    from mdscope.observables.lipid_enrichment import summarize_lipid_composition
    cfg=_cfg()
    cfg['analyses']['lipid_enrichment']['species_selections']=dict(
        cfg['analyses']['lipid_enrichment']['species_selections'], absent='resname NONE')
    cfg['analyses']['lipid_enrichment']['allow_absent_species']=True
    _execute(_universe(),cfg,tmp_path)
    summarize_lipid_composition(tmp_path)
    result=pd.read_parquet(tmp_path/'extensions/lipid_enrichment/composition_summary.parquet')
    assert result.loc[result.species.eq('absent'),'fraction_frames_available'].eq(0).all()
    assert result.loc[result.species.eq('POPC'),'mean_molecules'].eq(4).all()
    assert result.sampled_frames.eq(5).all()


def test_standalone_and_fused_area_paths_agree(tmp_path):
    from types import SimpleNamespace
    import MDAnalysis as mda
    from mdscope.observables.lipid_enrichment import run_lipid_enrichment
    u=_universe();top=tmp_path/'input.pdb';traj=tmp_path/'input.dcd'
    u.atoms.write(str(top))
    with mda.Writer(str(traj),n_atoms=len(u.atoms)) as writer:
        for ts in u.trajectory:
            writer.write(u.atoms)
    cfg=_cfg()
    cfg['analyses']['lipid_enrichment']['area_reference']=dict(
        exclusion_radius_A=3.,slab_half_width_A=5.,quadrature_order=32)
    _execute(_universe(),cfg,tmp_path/'fused')
    rep=SimpleNamespace(system='test',label='rep1',topology=top,trajectory=traj)
    run_lipid_enrichment(rep,cfg,tmp_path/'standalone')
    for name in ('timeseries.parquet','counts.parquet','selection_audit.parquet'):
        path='extensions/lipid_enrichment/'+name
        pd.testing.assert_frame_equal(pd.read_parquet(tmp_path/'fused'/path),
            pd.read_parquet(tmp_path/'standalone'/path),atol=1e-5,rtol=1e-6)


def test_area_and_composition_heatmaps_remain_distinct(tmp_path):
    from mdscope.reporting.membrane_figures import plot_lipid_enrichment_matrix
    path=tmp_path/'project_summary/statistics/time_series'
    path.mkdir(parents=True)
    pd.DataFrame({'system':['control']*2,'replicate':['rep1']*2,
        'metric':['lipid_enrichment_POPC_upper_fold','lipid_enrichment_POPC_upper_area_density_ratio'],
        'mean':[1.,2.]}).to_parquet(path/'replicate_statistics.parquet')
    plot_lipid_enrichment_matrix({'plots':{'dpi':30}},tmp_path)
    assert (tmp_path/'project_summary/figures/membrane/lipid_enrichment_heatmap.png').is_file()
    assert (tmp_path/'project_summary/figures/membrane/lipid_area_density_heatmap.png').is_file()


def test_fully_blocked_shell_is_unavailable_not_infinite():
    species=np.array(['A']*4)
    upper=np.array([True,True,False,False])
    positions=np.array([[21,20,30],[22,20,30],[21,20,10],[22,20,10]],float)
    rows=enrichment_frame(species,upper,positions,CENTER,BOX,0.,10.,2,
        obstacles=np.array([CENTER,[20.,20.,10.]]),
        area_reference=dict(exclusion_radius_A=15.,slab_half_width_A=5.,quadrature_order=32))
    assert all(row['count_local']==2 and row['shell_accessible_area_A2']==0 for row in rows)
    assert all(np.isnan(row['area_density_ratio']) and
        row['area_reference_status']=='unavailable_zero_accessible_area' for row in rows)
