#!/usr/bin/env python

"""
"""

import sys
import os
import numpy as np

import matplotlib.pyplot as plt

import PyCSPadImage.CalibParsDefault   as cald
import PyCSPadImage.CalibPars          as calp
import PyCSPadImage.CalibParsEvaluated as cpe
import PyCSPadImage.CSPadConfigPars    as ccp
import PyCSPadImage.CSPadImageProducer as cip
import PyCSPadImage.GlobalMethods      as gm # getCSPadArrayFromFile for pedestal subtraction 

import PyCSPadImage.GlobalGraphics     as gg # For test purpose in main only
import PyCSPadImage.HDF5Methods        as hm # For test purpose in main only


def get_event_from_hdf5(hd5_image_file, calibration_path, dsname, event):
    
    print 'Getting raw CSPad event %d from file %s \ndataset %s' % (event, fname, dsname)
    ds1ev = hm.getAverageCSPadEvent( fname, dsname, event, nevents=5 )
    if not ds1ev.shape == (32, 185, 388):
        print 'WARNING: ds1ev.shape =', ds1ev.shape, "should be (32, 185, 388)"
        
    return ds1ev
    
    
def get_event_from_npz(npz_image_file):
    return np.load(npz_image_file)['arr_0']
    

def assemble_image(raw_image, calibration_path, run_number):
    """
    This function takes a 'raw' hd5 image file and returns an array representing
    the assembled image.
    
    Parameters
    ----------
    raw_image : np.ndarray
        A shape=(32, 185, 388) array of the raw 1x1s
    
    Returns
    -------
    image : np.ndarray
        A two-dimensional image representing the intensities on the detector.
    """
    
    print 'Loading calibration parameters from: %s' % calibration_path 
    calp.calibpars.setCalibParsForPath( run=run_number, path=calibration_path )
    cpe.cpeval.printCalibParsEvaluatedAll() 

    print 'Constructing the CSPad image from raw array'
    cspadimg = cip.CSPadImageProducer(rotation=0, tiltIsOn=True)#, mirror=True)
    image = cspadimg.getCSPadImage( raw_image )


    return image
    
    
def plot_assembled_image(image):
    
    plt.imshow(image.T)
    plt.show()
        
    return


def main():
    """
    Do a basic test on one of our files. Uses Mikails geometry.
    """
    
    # change the paths below to run this script
    calibration_path = 'mikhail_geom'
    run_number = 58
    
    raw_image = get_event_from_npz('npz_examples/cxi64813_r58_evt1.npz')
    assembled_image = assemble_image(raw_image, calibration_path, run_number)
    plot_assembled_image(assembled_image)

if __name__ == "__main__":
    main()

