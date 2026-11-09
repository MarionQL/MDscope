"""Independent same-entity and exhaustive small circular-shift references."""
import numpy as np
import pandas as pd
import pytest
from mdscope.observables.solvation_postprocess import bridge_pair_evidence


def _events(path, rows):
    pd.DataFrame(rows,columns=['frame','solvent','target_id','solvent_id']).to_parquet(path)
    return [path]


def test_same_entity_trace_and_all_nonzero_shift_reference(tmp_path):
    paths=_events(tmp_path/'events.parquet',[(i,'ions/K',target,7)
        for i in (0,1) for target in ('acid','HS')]+[(0,'ions/K','acid',7)])
    trace,null,summary=bridge_pair_evidence(paths,[(i,float(i)) for i in range(4)],
        [('acid','HS')],['ions/K'])
    assert trace.bridge_present_fraction.tolist()==[1,1,0,0]
    assert trace.bridge_entity_count.tolist()==[1,1,0,0]
    np.testing.assert_allclose(null.sampled_bridge_fraction,[.25,0.,.25])
    assert summary.n_offsets.item()==3
    assert summary.shift_reference_mean.item()==pytest.approx(1/6)
    assert summary.observed_minus_shift_mean.item()==pytest.approx(1/3)
    assert summary.status.item()=='descriptive_shift_reference'


def test_different_entities_never_form_a_bridge_and_unobserved_pairs_stay_zero(tmp_path):
    paths=_events(tmp_path/'events.parquet',[(0,'ions/K','acid',7),(0,'ions/K','HS',8)])
    trace,null,summary=bridge_pair_evidence(paths,[(i,float(i)) for i in range(4)],
        [('acid','HS'),('acid','HS2')],['ions/K','ions/Na'])
    assert trace.bridge_present_fraction.eq(0).all()
    assert null.sampled_bridge_fraction.eq(0).all()
    assert summary.observed_minus_shift_mean.eq(0).all()


@pytest.mark.parametrize('times,status', [([0.,1.,3.,4.],'unavailable_irregular_sampling'),
                                         ([0.,1.,2.],'unavailable_short_window')])
def test_null_unavailable_for_irregular_or_short_windows(tmp_path,times,status):
    paths=_events(tmp_path/'events.parquet',[(0,'ions/K','acid',7),(0,'ions/K','HS',7)])
    trace,null,summary=bridge_pair_evidence(paths,list(enumerate(times)),[('acid','HS')],['ions/K'])
    assert trace.bridge_present_fraction.iloc[0]==1
    assert null.empty
    assert summary.status.item()==status and summary.shift_reference_mean.isna().all()


def test_multiple_same_entities_count_once_each_and_shift_together(tmp_path):
    paths=_events(tmp_path/'events.parquet',[(0,'water',target,entity)
        for target in ('acid','HS') for entity in (7,8)])
    trace,null,_=bridge_pair_evidence(paths,[(i,float(i)) for i in range(4)],
        [('acid','HS')],['water'])
    assert trace.bridge_present_fraction.tolist()==[1,0,0,0]
    assert trace.bridge_entity_count.tolist()==[2,0,0,0]
    assert null.sampled_bridge_fraction.eq(0).all()


def test_long_window_offsets_are_bounded_and_reproducible(tmp_path):
    paths=_events(tmp_path/'events.parquet',[(0,'ions/K','acid',7),(1,'ions/K','HS',7)])
    observations=[(i,float(i)) for i in range(300)]
    first=bridge_pair_evidence(paths,observations,[('acid','HS')],['ions/K'])
    second=bridge_pair_evidence(paths,observations,[('acid','HS')],['ions/K'])
    for a,b in zip(first,second):pd.testing.assert_frame_equal(a,b)
    assert first[2].n_offsets.item()==100
    assert first[2].observed_bridge_fraction.item()==0.
    assert first[2].shift_reference_mean.item()==pytest.approx(1/30000)


def test_observer_writes_site_audit_and_preserves_pair_summary(tmp_path):
    from mdscope.observables.solvation import SolvationObserver
    from test_solvation import _universe
    import json
    u=_universe()
    cfg={'analyses':{'solvation':dict(enabled=True,preset='kinetics',
        targets={'protein':{'first':'index 0','second':'index 1'}},
        water_selection='resname SOL and name OW',search_radius_A=5.,
        options={'bridge_target_pairs':[['protein/first/0','protein/second/0']]})}}
    observer=SolvationObserver(u,cfg,tmp_path)
    for i,ts in enumerate(u.trajectory):observer.update(ts,i,float(i))
    observer.flush(0);observer.finalize()
    directory=tmp_path/'extensions/solvation'
    trace=pd.read_parquet(directory/'bridge_pair_timeseries.parquet')
    old=pd.read_parquet(directory/'bridge_pair_summary.parquet')
    assert trace.bridge_present_fraction.mean()==old.sampled_bridge_fraction.item()
    assert trace.bridge_entity_count.mean()==old.mean_simultaneous_bridging_ions.item()
    audit=pd.read_parquet(directory/'target_atoms.parquet')
    assert audit.atom_name.tolist()==['OE1','OE2']
    assert audit.charge_e.isna().all() # Missing topology charge is never guessed.
    assert json.loads((directory/'metadata.json').read_text())['solvent_selection_audit']['water']['n_entities']==2


def test_disagreement_between_bridge_and_contact_records_is_detected(tmp_path):
    from mdscope.observables.solvation import SolvationObserver
    from test_solvation import _universe
    u=_universe()
    cfg={'analyses':{'solvation':dict(enabled=True,preset='kinetics',
        targets={'protein':{'first':'index 0','second':'index 1'}},
        water_selection='resname SOL and name OW',search_radius_A=5.,
        options={'bridge_target_pairs':[['protein/first/0','protein/second/0']]})}}
    observer=SolvationObserver(u,cfg,tmp_path)
    for i,ts in enumerate(u.trajectory):observer.update(ts,i,float(i))
    observer.flush(0)
    source=tmp_path/'cache/solvation/bridges-000000.parquet'
    data=pd.read_parquet(source)
    data.iloc[1:].to_parquet(source,index=False) # Remove one real same-entity contact.
    with pytest.raises(ValueError,match='disagree'):
        observer.finalize()


def test_spatial_reference_preserves_height_and_uses_periodic_same_ion(tmp_path):
    from mdscope.observables.solvation import SolvationObserver
    from mdscope.statistics.time_series import collect_time_series
    from test_solvation import _universe

    u = _universe()
    cfg = {"analyses": {"solvation": dict(
        enabled=True, preset="kinetics", search_radius_A=5.,
        targets={"protein": {"first": "index 0", "second": "index 1"}},
        ion_selections={"K": "index 4"}, options={
            "bridge_target_pairs": [["protein/first/0", "protein/second/0"]],
            "bridge_spatial_ion": "ions/K",
        })}}
    observer = SolvationObserver(u, cfg, tmp_path)
    first = u.trajectory[0]
    observer.update(first, 0, 0.)
    observer.flush(0)
    second = u.trajectory[1]
    u.atoms[4].position = [3 + 20/3, 1, 2]
    observer.update(second, 1, 5.)
    observer.flush(1)
    observer.finalize()

    directory = tmp_path / "extensions/solvation"
    reference = pd.read_parquet(directory / "bridge_spatial_reference.parquet")
    positive = reference[reference.bridge_entity_count > 0]
    assert positive.frame.tolist() == [1]
    assert positive.offset_u_fraction.item() == pytest.approx(2/3)
    assert positive.offset_v_fraction.item() == 0
    trace = pd.read_parquet(directory / "bridge_pair_timeseries.parquet")
    assert trace.bridge_present_fraction.tolist() == [1, 0]
    assert trace.observed_minus_spatial_bridge_fraction.tolist() == [1, -1/8]
    assert any(metric == "observed_minus_spatial_bridge_fraction"
               for _, metric in collect_time_series(tmp_path))
