
"""
Utilities for interfacing with psana
"""


import os
import numpy as np

import cspad
from read import enforce_raw_img_shape


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


def interp_to_polar(raw_image, q_values, dtc, num_phi=1024, mask=None):
    """
    Interpolate a single image into polar coordiantes. Requires thor.

    Parameters
    ----------
    raw_image : np.ndarray
        The raw pixel intensities, in shape (4, 16, 185, 194)

    q_values : np.ndarray
        A 1d list of the q-values to interpolate into

    dtc : thor.Detector
        A detector object

    Optional Parameters
    -------------------
    num_phi : int
        The number of phi slices to interpolate (more = higher res)

    mask : np.ndarray
        A pixel mask to apply

    Notes
    -----
    Not elegant code, but it works
    """
    try:
        import thor
    except ImportError as e:
        print e
        raise ImportError('Thor must be installed to employ the interp_to_polar'
                          'function...')

    if mask == None:
        mask = np.ones((4, 16, 185, 194))

    thor_format = enforce_raw_img_shape(raw_image).flatten()
    ss = thor.Shotset(thor_format, dtc, mask=mask.flatten())
    r = ss.to_rings(q_values, num_phi=num_phi)

    return r.polar_intensities[0] * r.polar_mask


