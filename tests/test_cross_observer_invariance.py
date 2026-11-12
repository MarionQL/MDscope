"""One source-space trajectory exercises geometry and contact observers together."""
import numpy as np
import pandas as pd
import pytest
from mdscope.pipeline.fused_membrane import FusedMembrane
from mdscope.observables.solvation import SolvationObserver
from test_fused_membrane import _universe, _cfg


@pytest.mark.parametrize('transform', ['translate_wrap', 'rotate_z'])
def test_cross_observer_source_geometry(tmp_path, transform):
    results = []
    for changed in (False, True):
        u, cfg = _universe(), _cfg()
        u.add_TopologyAttr("bonds", [])
        cfg['analyses']['lipid_contacts'] = dict(enabled=True, selection='resname POPC',
            species_selections={'POPC':'resname POPC and name P'}, cutoff_A=16., batch_atoms=3)
        cfg['analyses']['solvation'] = dict(enabled=True, water_selection='resname POPC',
            targets={'protein': {'CA': 'name CA'}}, search_radius_A=25., shell_cutoffs_A={'water':25.})
        if changed:
            xyz = u.trajectory.coordinate_array
            if transform == 'translate_wrap':
                xyz += [47., -23., 47.]
                xyz += (((np.arange(30).reshape(10,3) % 3)-1)*100.)[None]
            else:
                original = xyz.copy()
                xyz[...,0], xyz[...,1] = -original[...,1], original[...,0]
        root = tmp_path/str(changed)
        membrane = FusedMembrane(u, cfg, root)
        solvent = SolvationObserver(u, cfg, root)
        before = u.trajectory.coordinate_array.copy()
        for i, ts in enumerate(u.trajectory):
            membrane.update(ts, i)
            solvent.update(ts, i, float(i))
        membrane.flush(0); solvent.flush(0)
        membrane.finalize(); solvent.finalize()
        # Observer execution must leave source coordinates intact.
        np.testing.assert_array_equal(u.trajectory.coordinate_array, before)
        results.append(root)
    for relative in ('extensions/bilayer_normal/timeseries.parquet',
                     'extensions/lipid_footprint/counts.parquet',
                     'extensions/lipid_contacts/matrix.parquet',
                     'extensions/solvation/counts.parquet', 'extensions/solvation/rdf.parquet'):
        pd.testing.assert_frame_equal(pd.read_parquet(results[0]/relative),
            pd.read_parquet(results[1]/relative), check_exact=False, atol=1e-4, rtol=1e-5)
    # Four square edges per leaflet at 15 A; no diagonals or cross-leaflet
    # pairs meet the 16-A cutoff. Molecules are counted once, not twice.
    matrix = pd.read_parquet(results[0]/'extensions/lipid_contacts/matrix.parquet')
    assert matrix.pair_count_undirected.eq(4).all()
