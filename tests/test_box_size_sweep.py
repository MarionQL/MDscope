"""Verify candidate sampling classification without assuming 'clear' is box approval."""
import unittest
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from box_size_sweep import evaluate_candidates,placement_shift

class CandidateSizingTests(unittest.TestCase):
    def test_fragment_reimage_is_whole_cell_shift(self):
        box=np.array([20.,20.,24.])
        p=np.array([[21.,-1.,5.],[22.,0.,5.]])
        d=placement_shift(p,box)
        np.testing.assert_array_equal(d,[-20.,20.,0.])
    def test_candidate_lateral_clearance_changes_with_box(self):
        # Separate two assemblies; at L=10, cross copies are close, at L=20 clear.
        a=np.array([[1.,1.,2.]])
        b=np.array([[9.,1.,2.]])
        membrane=np.array([[1.,1.,0.]])
        # z nearest = TARGET_Z - 2, comfortably above cutoff.
        x=evaluate_candidates([a,b],membrane,[10.,20.],[])
        assert len(x)==2
        self.assertAlmostEqual(x[0]["min_cross_assembly_image_nm"],2.)
        self.assertGreater(x[1]["min_cross_assembly_image_nm"],5.)
    def test_vertical_image_clearance_changes(self):
        a=np.array([[1.,1.,12.]])
        b=np.array([[4.,1.,12.]])
        membrane=np.array([[1.,1.,0.]])
        x=evaluate_candidates([a,b],membrane,[],[14.,30.])
        self.assertAlmostEqual(x[0]["min_z_membrane_image_nm"],2.)
        self.assertGreater(x[1]["min_z_membrane_image_nm"],10.)
if __name__=="__main__":unittest.main()
