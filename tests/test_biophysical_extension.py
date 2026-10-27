"""Synthetic regression checks for geometry and diffusion without trajectory files."""
import numpy as np
from mdscope.analysis.biophysical_common import linear_diffusion,sequence_motifs
from mdscope.analysis.lipid_diffusion import diffusion_fit,multi_origin_msd
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

def test_multiorigin_linear_motion():
    x=np.zeros((11,2,2))
    x[:,:,0]=np.arange(11)[:,None]
    result=multi_origin_msd(x,np.array([1,2,4]))
    np.testing.assert_allclose(result,[1,4,16])

def test_diffusion_units():
    # Lateral Brownian slope of 0.4 Å²/ps = 1 nm²/ns
    t=np.arange(0,2000,20,dtype=float)
    m=.4*t
    assert np.isclose(diffusion_fit(t,m),1.)
    assert np.isclose(linear_diffusion(t,m),1.)

def test_correlation_absent_when_no_variation():
    assert np.isnan(correlation(np.ones(9),np.arange(9)))
    np.testing.assert_equal(block_means(np.arange(10),2),[.5,2.5,4.5,6.5,8.5])
