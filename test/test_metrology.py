
"""
This simple test assesses how close our "CSPad" object reproduces the precise
xyz positions of each pixel. The difference should be negligable.
"""

import numpy as np

from autogeom import cspad

class TestMetrology(object):
    
    def setup(self):
        self.met = cspad.Metrology('examples/cspad_2011-08-10-Metrology.txt')
        
    def test_cspad_metrology_consistency(self):
        
        csp = self.met.to_cspad()
        mp = self.met.pixel_positions
        cp = csp.pixel_positions
        
        assert mp.shape == cp.shape
        
        sum_err = np.sum( np.abs( mp - cp ) )
        
        tolerance = 1.0
        if sum_err > tolerance:
            raise ValueError("CSPad does not reproduce Metrology to an abs-sum "
                             "precision of %f. Error: %f" % (tolerance, sum_err))