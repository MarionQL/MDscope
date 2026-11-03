"""Independent brute-image references for bounding-box pruning."""
import sys
from itertools import product
from pathlib import Path
import unittest
import numpy as np
from scipy.spatial.distance import cdist
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from release_box_clearance import image_minimum

class ImageMinimumTests(unittest.TestCase):
    def test_random_image_search_exact_and_swapped(self):
        rng=np.random.default_rng(2026)
        box=np.array([17.,21.,25.])
        for _ in range(5):
            a=rng.uniform(-4.,23.,(9,3))
            b=rng.uniform(-7.,24.,(17,3))
            got=image_minimum(a,b,box)
            expected=min(float(cdist(a,b+np.array(t)*box).min()) for t in product((-2,-1,0,1,2),repeat=3) if t!=(0,0,0))
            self.assertAlmostEqual(got["distance_nm"],expected,places=10)
            swapped=image_minimum(b,a,box)
            self.assertAlmostEqual(swapped["distance_nm"],expected,places=10)
    def test_translated_z_membrane_and_unequal_group_sizes(self):
        box=np.array([20.,20.,25.])
        a=np.array([[2.,4.,24.],[2.,4.,23.]])
        b=np.array([[2.,4.,0.],[17.,15.,1.],[18.,15.,1.]])
        result=image_minimum(a,b,box,z_image_only=True)
        self.assertAlmostEqual(result["distance_nm"],1.)
        self.assertNotEqual(result["shift"][2],0)
if __name__=="__main__":unittest.main()
