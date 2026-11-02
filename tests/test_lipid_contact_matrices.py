"""Hand-counted undirected graphs, molecule identity, periodic geometry and replay."""
import numpy as np
import pandas as pd
import pytest
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from mdscope.observables.lipid_contacts import contact_matrix
from mdscope.pipeline.fused_membrane import FusedMembrane
from mdscope.observables.solvation_postprocess import bridge_residence


def test_complete_graph_has_random_mixing_enrichment_one():
    pairs = [(a,b) for a in range(4) for b in range(a+1,4)]
    table = contact_matrix(['A','A','B','B'], [True]*4, pairs, ['A','B','absent'])
    upper = table[table.leaflet.eq('upper')]
    observed = upper[(upper.molecules_a>0)&(upper.molecules_b>0)]
    assert observed.enrichment.eq(1).all()
    assert upper.query("species_a == 'A' and species_b == 'B'").neighbors_per_molecule.item()==2
    assert upper.query("species_a == 'A' and species_b == 'A'").neighbors_per_molecule.item()==1
    assert upper[upper.molecules_a.eq(0)].enrichment.isna().all()


def _system():
    u=mda.Universe.empty(7,n_residues=5,atom_resindex=[0,0,1,2,2,3,4],trajectory=True)
    for key, vals in [('names',['P','C','P','P','C','P','P']),('resnames',['AAA','AAA','BBB','BBB','BBB']),('resids',[7,9,13,21,22]),('bonds',[(0,1),(3,4)])]:
        u.add_TopologyAttr(key, vals)
    xyz=np.array([[.25,1,15],[.4,1,15],[19.75,1,15],[10,1,15],[10.1,1,15],[10,1,5],[16,1,5]],np.float32)
    frames=np.repeat(xyz[None],4,axis=0)
    frames[2,0,0]+=20 # independent atom wrapping preserves molecular contact
    u.load_new(frames,format=MemoryReader,dimensions=np.tile([20,20,20,90,90,90],(4,1)))
    cfg={'trajectory':{'stride':1},'time':{'dt_ps':10.},'analyses':{'lipid_contacts':{'enabled':True,'selection':'all','species_selections':{'A':'resname AAA and name P','B':'resname BBB and name P'},'cutoff_A':1.,'batch_atoms':1}}}
    return u,cfg


def _run(root,batch,chunk,resume=False):
    u,cfg=_system();cfg['analyses']['lipid_contacts']['batch_atoms']=batch
    o=FusedMembrane(u,cfg,root)
    if resume:
        for i in range(3):
            o.update(u.trajectory[i],i)
            if i==1:o.flush(0)
        o.flush(1)
        u.trajectory[0];o=FusedMembrane(u,cfg,root);o.prune_from(1)
        start,part=2,1
    else:start,part=0,0
    for i in range(start,4):
        o.update(u.trajectory[i],i)
        if (i+1)%chunk==0:o.flush(part);part+=1
    o.finalize()
    return pd.read_parquet(root/'extensions/lipid_contacts/matrix.parquet')


def test_periodic_molecule_contacts_batch_chunk_and_resume(tmp_path):
    reference=_run(tmp_path/'reference',100,4)
    for label,batch,chunk,resume in [('batched',1,1,False),('resumed',2,2,True)]:
        actual=_run(tmp_path/label,batch,chunk,resume)
        pd.testing.assert_frame_equal(actual,reference)
    aa=reference.query("leaflet == 'upper' and species_a == 'A' and species_b == 'A'")
    assert aa.pair_count_undirected.eq(1).all() # atom multiplicity does not count twice
    assert aa.neighbors_per_molecule.eq(1).all()
    assert aa.enrichment.eq(3).all() # one of 3 possible upper-leaflet pairs
    ab=reference.query("leaflet == 'upper' and species_a == 'A' and species_b == 'B'")
    assert ab.pair_count_undirected.eq(0).all() and ab.enrichment.eq(0).all()


def test_same_ion_bridge_episode_identity_and_intermittent_gaps():
    bridges=pd.DataFrame({'time_ps':[0.,10.,10.,30.], 'solvent':['ions/K']*4,
        'solvent_id':[91]*4,'target_a':['acidic_patch']*4,'target_b':['HS_A']*4})
    times=np.arange(4)*10.
    continuous=bridge_residence(bridges,times)
    assert continuous.duration_ps.tolist()==[20.,10.]
    assert continuous.n_frames.tolist()==[2,1]
    assert continuous.left_censored.tolist()==[True,False]
    assert continuous.right_censored.tolist()==[False,True]
    intermittent=bridge_residence(bridges,times,gap_frames=1)
    assert intermittent.duration_ps.item()==40.
    assert intermittent.occupancy_fraction.item()==.75
    assert intermittent.molecule_id.item()==91


def test_matrix_plot_family_disabled_before_parquet_read(tmp_path, monkeypatch):
    from mdscope.reporting.membrane_figures import plot_lipid_contacts
    rep = tmp_path/'condition/rep1'
    source = rep/'extensions/lipid_contacts/matrix.parquet'
    source.parent.mkdir(parents=True);source.write_bytes(b'disabled data must not be read')
    def forbidden(*args, **kwargs):
        raise AssertionError('disabled matrix plot read source')
    monkeypatch.setattr(pd,'read_parquet',forbidden)
    plot_lipid_contacts({'plots':{'families':{'lipid_contacts':False}},
                         'analyses':{'lipid_contacts':{'enabled':True}}}, tmp_path, {'condition':[rep]})
