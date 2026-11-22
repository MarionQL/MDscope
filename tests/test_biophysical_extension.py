"""Synthetic regression checks for geometry and diffusion without trajectory files."""
import numpy as np

from mdscope.observables.msd import fit_lateral_diffusion
from mdscope.statistics.coupling import block_means, correlation


def test_rmsd_separate_from_rmsf():
    from mdscope.observables.rmsd import compute_rmsd
    from mdscope.observables.rmsf import Welford
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
