"""Independent geometric density references; no bulk or binding inference."""
import json
from copy import deepcopy
import numpy as np
import pandas as pd
import pytest
from mdscope.observables.solvation import aqueous_slab_density, SolvationObserver, rdf_from_counts
from test_solvation import _universe


def test_slab_density_periodic_variable_volume_and_excluded_hemisphere():
    box = np.array([20., 20., 20., 90., 90., 90.])
    membrane = np.array([[0., 0., 19.], [0., 0., 1.]])
    solvent = np.array([[0., 0., 6.], [0., 0., 10.], [0., 0., 14.], [0., 0., 2.]])
    density, count, volume, _ = aqueous_slab_density(solvent, membrane, box, 5.)
    assert count == 3 and volume == 4000. and density == 3/4000.
    shifted = solvent + [20., -40., 60.]
    np.testing.assert_allclose(aqueous_slab_density(shifted, membrane+box[:3], box, 5.)[:3],
                               [density, count, volume])
    bigger = box.copy(); bigger[0] = 40.
    assert aqueous_slab_density(solvent, membrane, bigger, 5.)[0] == density/2
    # Aqueous density changes the reference only. A half-accessible shell
    # still yields 0.5 at that reference density, rather than g=1.
    shell_volume = 4*np.pi/3
    np.testing.assert_allclose(rdf_from_counts([.5*density*shell_volume], [0., 1.], density), [.5])
    np.testing.assert_allclose(rdf_from_counts([.5*density*shell_volume], [0., 1.], density/2), [1.])


@pytest.mark.parametrize('half_width', [0., 10., np.nan])
def test_invalid_or_empty_regions(half_width):
    with pytest.raises(ValueError, match='half_width'):
        aqueous_slab_density([[0.,0.,0.]], [[0.,0.,0.]], [20.,20.,20.,90.,90.,90.], half_width)


def test_empty_reference_and_unsupported_geometry():
    density, count, _, _ = aqueous_slab_density([[0.,0.,1.]], [[0.,0.,0.]], [20.,20.,20.,90.,90.,90.], 5.)
    assert count == 0 and np.isnan(density)
    with pytest.raises(ValueError, match='orthorhombic'):
        aqueous_slab_density([[0.,0.,1.]], [[0.,0.,0.]], [20.,20.,20.,90.,90.,80.], 5.)


def _run(root, *, reference, empty_second=False):
    u = _universe()
    # Slab around the target is a mathematical fixture, not a production
    # membrane selection. Count reference solvent at |z-1| > 4 in a 20-A cell.
    coordinates = u.trajectory.coordinate_array
    coordinates[:, 2, 2] = 5.5
    coordinates[:, 3, 2] = 10.
    if empty_second:
        coordinates[1, 2:4, 2] = 1.
    spec = dict(enabled=True, targets={'protein': {'site': 'resname GLU'}},
                water_selection='resname SOL and name OW', search_radius_A=5.)
    if reference:
        spec['aqueous_reference'] = dict(membrane_selection='resname GLU', half_width_A=4.)
    cfg = {'analyses': {'solvation': spec}}
    observer = SolvationObserver(u, cfg, root)
    for i, ts in enumerate(u.trajectory):
        observer.update(ts, i, float(i))
        observer.flush(i)
    observer.finalize()
    return root/'extensions/solvation'


def test_observer_preserves_legacy_curve_and_records_same_selection_density(tmp_path):
    old = _run(tmp_path/'old', reference=False)
    new = _run(tmp_path/'new', reference=True)
    legacy, updated = pd.read_parquet(old/'rdf.parquet'), pd.read_parquet(new/'rdf.parquet')
    pd.testing.assert_frame_equal(legacy, updated[legacy.columns])
    # Two selected waters in the 20*20*(20-8) reference region.
    # Full-cell density is 2/8000; new exposure is 2/4800.
    assert updated.pair_count.sum() > 0
    np.testing.assert_allclose(updated.g_aqueous_reference, updated.g_bulk*.6)
    audit = pd.read_parquet(new/'aqueous_reference.parquet')
    assert audit.n_reference.tolist() == [2,2]
    assert audit.reference_volume_A3.tolist() == [4800.,4800.]
    assert json.loads((new/'metadata.json').read_text())['aqueous_reference']['half_width_A'] == 4.
    empty = _run(tmp_path/'empty', reference=True, empty_second=True)
    assert pd.read_parquet(empty/'rdf.parquet').g_aqueous_reference.isna().all()
    assert pd.read_parquet(empty/'rdf.parquet').g_bulk.notna().all()


def test_aqueous_statistics_preserve_missing_replica_coverage(tmp_path):
    from mdscope.config import DEFAULTS
    from mdscope.statistics.solvation import analyze_solvation_rdfs
    paths = [_run(tmp_path/'system/rep1', reference=True),
             _run(tmp_path/'system/rep2', reference=False)]
    reps = [path.parents[1] for path in paths]
    cfg = deepcopy(DEFAULTS)
    analyze_solvation_rdfs(cfg, tmp_path, {'system': reps})
    system = pd.read_parquet(tmp_path/'project_summary/statistics/solvation/system_rdf.parquet')
    assert system.n_aqueous_replicates.eq(1).all()
    assert not system.complete_aqueous_replica_coverage.any()
    assert system.complete_replica_coverage.all()
