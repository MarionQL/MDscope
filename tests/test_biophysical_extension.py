"""Synthetic regression checks for geometry and diffusion without trajectory files."""
import numpy as np
from mdscope.analysis.biophysical_common import sequence_motifs
from mdscope.analysis.msd import fit_lateral_diffusion
from mdscope.analysis.coupling_screen import correlation,block_means

class Residue:
    def __init__(self,name,resid,ix,chain):
        self.resname=name;self.resid=resid;self.ix=ix;self.chainID=chain
def test_motif_matches_and_chain_boundaries():
    seq=[Residue(n,i,i,"A") for i,n in enumerate(["ASN","HIS","ILE","PRO"],1)]
    seq += [Residue(n,i+4,i+4,"B") for i,n in enumerate(["ASN","HIS","ILE","PRO"],1)]
    hits=sequence_motifs(seq,{"nxip":"N.IP"})
    assert len(hits)==2
    assert {m["chain"] for m in hits}=={"A","B"}

def test_rmsd_separate_from_rmsf():
    from mdscope.analysis.rmsd import compute_rmsd
    from mdscope.analysis.rmsf import Welford
    x=np.array([[0.,0.,0.],[2.,0.,0.]])
    assert np.isclose(compute_rmsd(x,x),0)
    w=Welford(x.shape)
    w.update(x)
    w.update(x+1)
    np.testing.assert_allclose(w.rmsf(),np.sqrt(3)/2*np.ones(2))

def test_diffusion_units():
    # Lateral Brownian slope of 0.4 Å²/ps = 1 nm²/ns
    t=np.arange(0,2000,20,dtype=float)
    m=.4*t
    assert np.isclose(fit_lateral_diffusion(t,m)["diffusion_nm2_ns"],1.)

def test_correlation_absent_when_no_variation():
    assert np.isnan(correlation(np.ones(9),np.arange(9)))
    np.testing.assert_equal(block_means(np.arange(10),2),[.5,2.5,4.5,6.5,8.5])
