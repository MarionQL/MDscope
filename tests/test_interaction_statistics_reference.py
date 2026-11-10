"""Six-replica nonconstant interaction statistics with explicit absent events."""
from copy import deepcopy
import numpy as np
import pandas as pd
import pytest
from mdscope.config import DEFAULTS
from mdscope.statistics.interactions import analyze_interactions


def fixture(root, design='independent'):
    paths={'a':[],'b':[]}
    for system in paths:
        for rep,frames in enumerate((8,11,7),1):
            dest=root/system/f'rep{rep}'
            (dest/'summary').mkdir(parents=True)
            times=np.arange(frames+(system=='b'))*10.
            pd.DataFrame({'time_ps':times}).to_parquet(dest/'summary/timeseries.parquet')
            data=dest/'interactions/salt_bridges';data.mkdir(parents=True)
            present=[];life=[]
            for identity in range(3):
                # Last interaction is absent in replicas 1 and 2.
                active=[1,2] if identity==0 else [3] if identity==1 else [4] if rep==3 else []
                keys={'residue_1':f'ASP{identity}', 'residue_2':'LYS7','interaction_class':'interchain'}
                present.extend(dict(keys,time_ps=float(times[i])) for i in active)
                if active:
                    life.append(dict(keys,start_ps=float(times[active[0]]),lifetime_ps=len(active)*10.,left_censored=False,completed=(identity==0)))
            pd.DataFrame(present).to_parquet(data/'presence.parquet')
            pd.DataFrame(life).to_parquet(data/'lifetimes.parquet')
            paths[system].append(dest)
    cfg=deepcopy(DEFAULTS)
    cfg['statistics']['bootstrap_draws']=100
    cfg['statistics']['hypothesis_tests']['design']=design
    cfg['progress']['enabled']=False
    return cfg,paths


@pytest.mark.parametrize('design',['independent','paired'])
def test_interaction_points_missing_vs_zero_and_reproducibility(tmp_path,design):
    cfg,paths=fixture(tmp_path/'input',design)
    analyze_interactions(cfg,tmp_path/'first',paths)
    analyze_interactions(cfg,tmp_path/'second',paths)
    output='project_summary/statistics/interactions'
    first=tmp_path/'first'/output;second=tmp_path/'second'/output
    for path in first.glob('*.parquet'):
        pd.testing.assert_frame_equal(pd.read_parquet(path),pd.read_parquet(second/path.name),check_exact=True)
    rep=pd.read_parquet(first/'replicate_interaction_metrics.parquet')
    absent=rep[(rep.replicate=='rep1')&rep.interaction.str.startswith('ASP2')]
    assert absent.loc[absent.metric.eq('occupancy'),'value'].eq(0).all()
    assert absent.loc[absent.metric.eq('mean_lifetime_ps'),'value'].isna().all()
    known=rep[(rep.system=='a')&(rep.replicate=='rep1')&rep.interaction.str.startswith('ASP0')]
    assert known.loc[known.metric.eq('occupancy'),'value'].item()==.25
    assert known.loc[known.metric.eq('mean_lifetime_ps'),'value'].item()==20.


def test_occupancy_screen_keeps_boundary_and_all_replicas(tmp_path):
    cfg, paths = fixture(tmp_path / "input")
    cfg["statistics"]["interaction_min_occupancy"] = 0.20
    for system, replicas in paths.items():
        for rep in replicas:
            path = rep / "interactions/salt_bridges/presence.parquet"
            present = pd.read_parquet(path)
            if system == "a" and rep.name == "rep1":
                pd.DataFrame({"time_ps": np.arange(10) * 10.}).to_parquet(
                    rep / "summary/timeseries.parquet"
                )
            else:
                present = present[~(present.residue_1.eq("ASP0") & present.time_ps.eq(20.))]
            # Repeating rows for a single frame cannot increase occupancy.
            duplicate = present[present.residue_1.eq("ASP1")]
            pd.concat([present, duplicate, duplicate]).to_parquet(path, index=False)

    analyze_interactions(cfg, tmp_path / "output", paths)
    table = pd.read_parquet(
        tmp_path / "output/project_summary/statistics/interactions/replicate_interaction_metrics.parquet"
    )
    assert table.interaction.str.startswith("ASP0").all()
    occupancy = table[table.metric.eq("occupancy")].set_index(["system", "replicate"]).value
    assert len(occupancy) == 6
    assert occupancy.loc[("a", "rep1")] == pytest.approx(0.20)
    assert occupancy.drop(("a", "rep1")).lt(0.20).all()
