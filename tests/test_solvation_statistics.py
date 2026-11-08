"""Solvation statistics use existing bootstrap machinery, not a second estimator."""
from pathlib import Path

import numpy as np
import pandas as pd

from mdscope.statistics.solvation import analyze_solvation_rdfs
from mdscope.statistics.time_series import collect_time_series, analyze_time_series
from mdscope.observables.solvation_postprocess import aggregate_coordination


def test_hydration_frame_fraction_and_mean():
    data = pd.DataFrame({
        "time_ps": [0.,0.,5.,5.], "solvent": ["water"]*4,
        "target_class":["hs"]*4, "site_class":["oxygen"]*4,
        "target_id":["a","b","a","b"], "n_solvent":[0,4,3,0]})
    frame, targets = aggregate_coordination(data)
    np.testing.assert_allclose(frame.hydrated_fraction, [.5,.5])
    np.testing.assert_allclose(frame.mean_n_solvent, [2.,1.5])
    assert targets.set_index("target_id").loc["a","hydrated_fraction"] == .5


def test_collector_reads_solvation_as_separate_named_groups(tmp_path):
    rep = tmp_path/"rep1"
    path = rep/"extensions"/"solvation"/"species_timeseries.parquet"
    path.parent.mkdir(parents=True)
    pd.DataFrame({
        "time_ps":[0.,5.,10.,0.,5.,10.],
        "solvent":["water"]*3+["ions/K"]*3,
        "target_class":["protein"]*6,"site_class":["carboxylate"]*6,
        "mean_n_solvent":[3.,4.,2.,0.,1.,0.],
        "hydrated_fraction":[1.,1.,1.,0.,.5,0.]
    }).to_parquet(path)
    series=collect_time_series(rep)
    group="solvation/water/protein/carboxylate"
    assert (group,"mean_n_solvent") in series
    assert (group,"hydrated_fraction") in series
    assert ("solvation/ions/K/protein/carboxylate","mean_n_solvent") in series
    assert series[(group,"mean_n_solvent")].value.tolist()==[3.,4.,2.]


def test_existing_block_bootstrap_system_comparison_for_solvation(tmp_path):
    outputs={}
    for system, offset in (("glyco",1.),("nonglyco",0.)):
        paths=[]
        for j in range(3):
            rep=tmp_path/system/f"rep{j+1}"
            source=rep/"extensions"/"solvation"/"species_timeseries.parquet"
            source.parent.mkdir(parents=True)
            # Correlated but deterministic time series; the same existing
            # moving block and hierarchical draws are applied unchanged.
            t=np.arange(32,dtype=float)*5
            pd.DataFrame({
                "time_ps":t, "solvent":"water", "target_class":"protein",
                "site_class":"acidic", "mean_n_solvent":
                    2+offset+.1*j+.05*np.sin(np.arange(32)/3),
                "hydrated_fraction":.4+.1*offset+.01*j
            }).to_parquet(source)
            (rep/"summary").mkdir(parents=True)
            paths.append(rep)
        outputs[system]=paths
    cfg={"statistics":{
        "enabled":True,"bootstrap_draws":100,"ci":.95,"seed":10,
        "block":{"length":4,"max_acf_lag":12},
        "hierarchy":{"block_within_replicate":True,"resample_replicates":True},
        "hypothesis_tests":{"enabled":False,"design":"independent",
                            "multiple_testing":"none"}
    }}
    analyze_time_series(cfg,tmp_path,outputs)
    directory=tmp_path/"project_summary"/"statistics"/"time_series"
    per_rep=pd.read_parquet(directory/"replicate_statistics.parquet")
    per_system=pd.read_parquet(directory/"system_statistics.parquet")
    comparison=pd.read_parquet(directory/"system_comparisons.parquet")
    mask=(per_rep.group=="solvation/water/protein/acidic") & (
        per_rep.metric=="mean_n_solvent")
    assert mask.sum()==6
    assert (per_rep.loc[mask,"block_length_frames"]==4).all()
    assert len(per_system.loc[
        (per_system.group=="solvation/water/protein/acidic") &
        (per_system.metric=="mean_n_solvent")])==2
    delta=comparison.loc[
        (comparison.group=="solvation/water/protein/acidic") &
        (comparison.metric=="mean_n_solvent")].iloc[0]
    assert np.isfinite(delta.mean_difference_ci_low)
    assert np.isfinite(delta.mean_difference_ci_high)


def test_rdf_across_replicates_is_descriptive(tmp_path):
    outputs={"glyco":[]}
    for j in range(3):
        rep=tmp_path/"glyco"/f"rep{j+1}"
        source=rep/"extensions"/"solvation"/"rdf.parquet"
        source.parent.mkdir(parents=True)
        pd.DataFrame({
            "solvent":["water"]*3,"target_class":["protein"]*3,
            "site_class":["acidic"]*3, "r_inner_A":[0.,1.,2.],
            "r_outer_A":[1.,2.,3.],
            "g_bulk":[1.+j,2.+j,3.+j]
        }).to_parquet(source)
        outputs["glyco"].append(rep)
    analyze_solvation_rdfs({"statistics":{"enabled":True}},tmp_path,outputs)
    result=pd.read_parquet(tmp_path/"project_summary"/"statistics"/"solvation"/"system_rdf.parquet")
    assert result.n_replicates.tolist()==[3]*3
    np.testing.assert_allclose(result.mean_g_bulk,[2.,3.,4.])
    np.testing.assert_allclose(result.sd_g_bulk,[1.,1.,1.])
    assert result.complete_replica_coverage.all()



def test_rdf_coverage_counts_entirely_missing_replica(tmp_path):
    source = tmp_path / "system/rep1/extensions/solvation/rdf.parquet"
    source.parent.mkdir(parents=True)
    pd.DataFrame({"solvent": ["water"], "target_class": ["protein"],
                  "site_class": ["acidic"], "r_inner_A": [0.],
                  "r_outer_A": [1.], "g_bulk": [2.]}).to_parquet(source)
    analyze_solvation_rdfs({}, tmp_path,
                          {"system": [tmp_path/"system/rep1", tmp_path/"system/rep2"]})
    data = pd.read_parquet(tmp_path/"project_summary/statistics/solvation/system_rdf.parquet")
    assert data.n_replicates.item() == 1
    assert data.expected_replicates.item() == 2
    assert not data.complete_replica_coverage.item()


def test_bridge_traces_reuse_replicate_statistics_without_pooling_entities(tmp_path):
    outputs={}
    for system,hit in [('control',0),('bridged',1)]:
        outputs[system]=[]
        for i in range(3):
            rep=tmp_path/system/f'rep{i+1}';outputs[system].append(rep)
            source=rep/'extensions/solvation/bridge_pair_timeseries.parquet'
            source.parent.mkdir(parents=True);(rep/'summary').mkdir()
            pd.DataFrame({'time_ps':np.arange(12,dtype=float), 'solvent':'ions/K',
                'target_a':'HS/site/2','target_b':'protein/acid/1',
                'bridge_entity_count':hit*3,'bridge_present_fraction':hit}).to_parquet(source)
    cfg={'statistics':{'enabled':True,'bootstrap_draws':20,'ci':.95,'seed':10,
        'block':{'length':2},'hypothesis_tests':{'enabled':False}}}
    analyze_time_series(cfg,tmp_path,outputs)
    dest=tmp_path/'project_summary/statistics/time_series'
    reps=pd.read_parquet(dest/'replicate_statistics.parquet')
    assert len(reps)==12 # two metrics x three independent replicas x two systems
    assert reps.loc[reps.system.eq('bridged') & reps.metric.eq('bridge_entity_count'),'mean'].eq(3.).all()
    differences=pd.read_parquet(dest/'system_comparisons.parquet')
    row=differences.loc[differences.metric.eq('bridge_present_fraction')].iloc[0]
    assert abs(row.mean_difference_system2_minus_system1)==1.
    assert pd.isna(row.mean_q_value) # No shift-null or frame-level significance is invented.
