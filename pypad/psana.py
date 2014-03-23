
"""
Utilities for interfacing with psana
"""


import os
import numpy as np

import cspad



def get_geometry(experiment, run, camera, energy, distance, dilation,
                 distance_correction=0.0, calibration_dilation=0.0):
    """
    Obtain a geometry & mask from TJL's standard location:

    /reg/d/psdm/xxx/xxxYYYYY/res/geometries/
    /reg/d/psdm/xxx/xxxYYYYY/res/masks/

    Parameters
    ----------
    experiment : str
        The 8 charater experiment identifier (e.g. cxii0114)

    run : int
        The run to acquire the geometry for

    camera : str
        The camera identifier (e.g. 'ds1')

    Returns
    -------
    dtc : thor.Detector
        A Thor detector object

    geom : pypad.CSPad
        A CSPAD object

    mask : np.ndarray
        A (4, 16, 185, 194) shape boolean mask (one = keep)
    """

    prefix = '/reg/d/psdm/%s/%s/res/geometries/%s/' % (experiment[:3], experiment, camera)
    geom = cspad.CSPad.load( os.path.join(prefix, 'current.cspad') )
    mask = np.load( os.path.join(prefix, 'mask.npy') )

    corr_dist  = distance + distance_correction
    delta_dilt = dilation - calibration_dilation

    geom.dilate(delta_dilt)

    dtc = geom.to_thor(energy, corr_dist) # eV, mm
    #print camera, corr_dist, delta_dilt, energy, dtc.q_max

    return dtc, geom, mask

