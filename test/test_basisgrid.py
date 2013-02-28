
import numpy as np
from numpy.testing import assert_array_almost_equal
from autogeom import cspad

class TestBasisGrid(object):
    
    def setup(self):
        self.p = np.array([0.0, 0.0, 1.0])
        self.s = np.array([1.0, 0.0, 0.0])
        self.f = np.array([0.0, 2.0, 0.0])
        self.shape = (10, 10)
        self.grid_list = [(self.p, self.s, self.f, self.shape)]
        self.bg = cspad.BasisGrid(self.grid_list)
        
    def test_add_grid(self):
        nbg = cspad.BasisGrid()
        nbg.add_grid(*self.grid_list[0])
        assert_array_almost_equal(nbg.to_explicit(), self.bg.to_explicit())
    
    def test_add_using_center(self):
        center = np.array([4.5, 9, 1.0])
        nbg = cspad.BasisGrid()
        nbg.add_grid_using_center(center, self.s, self.f, self.shape)
        assert_array_almost_equal(nbg.to_explicit(), self.bg.to_explicit())
    
    def test_get_grid(self):
        assert self.bg.get_grid(0) == self.grid_list[0]
        
    def test_to_explicit(self):
        ref = np.zeros((100,3))
        mg = np.mgrid[0:9:10j,0:18:10j]
        ref[:,0] = mg[0].flatten()
        ref[:,1] = mg[1].flatten()
        ref[:,2] = 1.0
        assert_array_almost_equal(self.bg.to_explicit(), ref)
        
    def test_grid_as_explicit(self):
        ref = np.zeros((10,10,3))
        mg = np.mgrid[0:9:10j,0:18:10j]
        ref[:,:,0] = mg[0]
        ref[:,:,1] = mg[1]
        ref[:,:,2] = 1.0
        assert_array_almost_equal(self.bg.grid_as_explicit(0), ref)