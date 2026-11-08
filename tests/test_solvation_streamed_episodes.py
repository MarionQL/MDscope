"""Entry-point memory regression and independent sampled-episode references."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from mdscope.observables.solvation import SolvationObserver, contact_episodes
from test_solvation import _universe


def test_finalization_never_materializes_the_full_event_table(tmp_path, monkeypatch):
    cfg={'analyses':{'solvation':dict(enabled=True,preset='kinetics',
        targets={'protein':{'site':'resname GLU'}},water_selection='resname SOL and name OW',
        search_radius_A=5.)}}
    u=_universe();observer=SolvationObserver(u,cfg,tmp_path)
    for i,ts in enumerate(u.trajectory):
        observer.update(ts,i,float(i));observer.flush(i)
    original=pd.read_parquet
    def bounded_read(path,*args,**kwargs):
        if Path(path).name.startswith('events-') or Path(path).name=='events.parquet':
            raise AssertionError('full contact-event table materialized')
        return original(path,*args,**kwargs)
    monkeypatch.setattr(pd,'read_parquet',bounded_read)
    observer.finalize()
    episodes=original(tmp_path/'extensions/solvation/episodes.parquet')
    assert len(episodes)==2
    assert episodes.duration_ps.tolist()==[0.,1.]


def _write_events(path, records):
    pd.DataFrame(records,columns=['frame','solvent','target_id','solvent_id']).to_parquet(path)
    return path


@pytest.mark.parametrize('frames', [[(0,0.),(2,2.),(4,7.),(6,11.)],[(0,0.)]])
def test_streamed_bounds_match_independent_in_memory_oracle(tmp_path,frames):
    from mdscope.observables.solvation_postprocess import write_contact_episodes
    ids=[frame for frame,_ in frames]
    events=[]
    for i,frame in enumerate(ids):
        events.extend([(frame,'water','A',7)]*(1+(i%2)))
        if i%2==0:events.append((frame,'ions/K','B',8))
    split=max(1,len(events)//2)
    paths=[_write_events(tmp_path/f'events-{i}.parquet',part)
           for i,part in enumerate([events[:split],events[split:]])]
    output=tmp_path/'episodes.parquet'
    write_contact_episodes(paths,frames,output)
    expected=contact_episodes(events,frames)
    actual=pd.read_parquet(output)
    keys=['solvent','target_id','solvent_id','start_frame']
    pd.testing.assert_frame_equal(actual.sort_values(keys).reset_index(drop=True),
                                 expected.sort_values(keys).reset_index(drop=True))
    if len(frames)==1:
        assert actual.duration_ps.eq(0).all() and actual.sampled_duration_ps.isna().all()


def test_unsorted_or_out_of_range_input_preserves_previous_output(tmp_path):
    from mdscope.observables.solvation_postprocess import write_contact_episodes
    output=tmp_path/'episodes.parquet';output.write_bytes(b'previous committed product')
    path=tmp_path/'events.parquet'
    for records,diagnostic in [([(1,'water','A',7),(0,'water','A',7)],'ordered'),
                               ([(9,'water','A',7)],'outside'),
                               ([(0,'water','A',entity) for entity in range(5000)]+
                                [(2,'water','A',7),(1,'water','A',7)],'ordered')]:
        _write_events(path,records)
        with pytest.raises(ValueError,match=diagnostic):
            write_contact_episodes([path],[(0,0.),(1,2.),(2,5.)],output)
        assert output.read_bytes()==b'previous committed product'
    assert len(list(tmp_path.glob('*.parquet')))==2


def test_multiple_output_batches_preserve_irregular_terminal_imputation(tmp_path):
    from mdscope.observables.solvation_postprocess import write_contact_episodes
    events=[(frame,'water','A',entity) for frame in (0,2) for entity in range(5000)]
    path=_write_events(tmp_path/'events.parquet',events)
    write_contact_episodes([path],[(0,0.),(1,3.),(2,10.)],tmp_path/'episodes.parquet')
    actual=pd.read_parquet(tmp_path/'episodes.parquet')
    assert len(actual)==10000
    assert actual.duration_ps.eq(0.).all()
    assert actual.loc[actual.start_frame.eq(0),'sampled_duration_ps'].eq(3.).all()
    assert actual.loc[actual.start_frame.eq(2),'sampled_duration_ps'].eq(5.).all()
