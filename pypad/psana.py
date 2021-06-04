
"""
Utilities for interfacing with psana
"""


import os
import numpy as np

from . import cspad
from .mask import PadMask
from .read import enforce_raw_img_shape


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

    try:
        padmask = PadMask.load(os.path.join(prefix, 'current.mask'))
        m = padmask.mask

        asic_pix = 185 * 194
        mask = np.zeros(4 * 16 * asic_pix, dtype=np.int8)

        for i in range(4):
            for j in range(16):
                s = (i*16 + j) * asic_pix
                e = (i*16 + j + 1) * asic_pix
                mask[s:e] = m[i,j,:,:].flatten()

    except:
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
        print(e)
        raise ImportError('Thor must be installed to employ the interp_to_polar'
                          'function...')

    if mask == None:
        mask = np.ones((4, 16, 185, 194))

    thor_format = enforce_raw_img_shape(raw_image).flatten()
    ss = thor.Shotset(thor_format, dtc, mask=mask.flatten())
    r = ss.to_rings(q_values, num_phi=num_phi)

    return r.polar_intensities[0] * r.polar_mask


def thor_to_psana(thor_fmt_intensities):

    if thor_fmt_intensities.shape == (4, 16, 185, 194):
        ii = thor_fmt_intensities
    elif thor_fmt_intensities.shape == (2296960,):
        ii = thor_fmt_intensities.reshape((4, 16, 185, 194))
    else:
        raise ValueError('did not understand intensity shape: '
                         '%s' % str(thor_fmt_intensities.shape))

    ix = np.zeros((32, 185, 388), dtype=thor_fmt_intensities.dtype)
    for i in range(4):
        for j in range(8):
            a = i * 8 + j
            ix[a,:,:] = np.hstack(( ii[i,j*2,:,:], ii[i,j*2+1,:,:] ))

    return ix


