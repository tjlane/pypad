#!/usr/bin/env python

"""
A high-level interface to pyana for assembling images from psana geometry
parameters. Main function to pay attention to is `assemble_image`.
"""

import sys
import os
from os.path import join as pjoin
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


# these are the parameters pyana looks for to assemble the image
# values here are the size the nparrays representing these parameters

# this is a minimal set that will probably work
# essential_params = {'center':                 (12, 8),
#                     'center_corr':            (12, 8),
#                     'marg_gap_shift':         (3, 4),
#                     'offset':                 (3, 4),
#                     'offset_corr':            (3, 4),
#                     'quad_rotation':      (4,),
#                     'quad_tilt':          (4,),
#                     'rotation':           (4, 8),
#                     'tilt':               (4, 8)}
                    
                    
# this is a full list of the parameters
essential_params = {'center' : (12, 8),
                    'center_corr' : (12, 8),
                    'common_mode' : (3,),
                    'filter' : (3,),
                    'marg_gap_shift' : (3, 4),
                    'offset' : (3, 4),
                    'offset_corr' : (3, 4),
                    'pedestals' : (5920, 388),
                    'pixel_status' : (5920, 388),
                    'quad_rotation' : (4,),
                    'quad_tilt' : (4,),
                    'rotation' : (4, 8),
                    'tilt' : (4, 8)}


def check_param_dict(param_dict):
    """
    Does a sanity check on a parameter dictionary.
    """
    
    if not param_dict.keys().sort() == essential_params.keys().sort():
        raise ValueError('keys of `param_dict` do not match expected parameter names')
        
    for key in param_dict.keys():
        if not param_dict[key].shape == essential_params[key]:
            raise ValueError('value of %s does not match expected shape: %s' % (key, str(essential_params[key])))
        

    return

def assemble_image_from_dir(raw_image, calibration_path, run_number=0):
    """
    This function takes a 'raw' hd5 image file and returns an array representing
    the assembled image. Uses parameters stored in a directory to do this.
    
    Parameters
    ----------
    raw_image : np.ndarray
        A shape=(32, 185, 388) array of the raw 1x1s
    
    calibration_path : str
        The path to a directory containing all the calibration parameters,
        which are in flat text files named runX-runY.data.
        
    run_number : int
        The run number of the image. This indicates which parameters in
        `calibration_path` should be used in the case that parameters changed
        over the course of the runs. Defaults to run 0.
    
    Returns
    -------
    image : np.ndarray
        A two-dimensional image representing the intensities on the detector.
    """
    
    print 'Loading calibration parameters from: %s' % calibration_path 
    calp.calibpars.setCalibParsForPath( run=run_number, path=calibration_path )
    calp.calibpars.printCalibPars()
    cpe.cpeval.printCalibParsEvaluatedAll() 

    print 'Constructing the CSPad image from raw array'
    cspadimg = cip.CSPadImageProducer(rotation=0, tiltIsOn=True )
    image = cspadimg.getCSPadImage( raw_image )


    return image
    
    
def assemble_image_from_params(raw_image, param_dict):
    """
    This function takes a 'raw' hd5 image file and returns an array representing
    the assembled image. Uses parameters passed as a dictionary to do this.
    
    Parameters
    ----------
    raw_image : np.ndarray
        A shape=(32, 185, 388) array of the raw 1x1s
    
    param_dict : dict
        A dictionary containing arrays of all the parameter types. Must have
        keys as indicated below, with corresponding values numpy arrays of
        the size indicated
        
        key                 np array shape
        ---                 --------------
        center 		        (12, 8)
        center_corr 		(12, 8)
        marg_gap_shift 		(3, 4)
        offset 		        (3, 4)
        offset_corr 		(3, 4)
        quad_rotation 		(4,)
        quad_tilt 		    (4,)
        rotation 		    (4, 8)
        tilt 		        (4, 8)
        
    Returns
    -------
    image : np.ndarray
        A two-dimensional image representing the intensities on the detector.
    """
    
    check_param_dict(param_dict)
    

    calp.calibpars.cpars.update(param_dict)
    
    # print so hopefully the user can see if something went wrong
    calp.calibpars.printCalibPars()
    cpe.cpeval.printCalibParsEvaluatedAll() 
    
    cspadimg = cip.CSPadImageProducer(rotation=0, tiltIsOn=True )
    image = cspadimg.getCSPadImage( raw_image )
    
    return image
    
    
def write_params_to_dir(param_dict, dir_name):
    """
    Write the parameters in `param_dict` to disk in pyana's dir layout/format.
    
    Parameters
    ----------
    param_dict : dict
        A dictionary containing arrays of all the parameter types. Must have
        keys as indicated below, with corresponding values numpy arrays of
        the size indicated
        
        key                 np array shape
        ---                 --------------
        center 		        (12, 8)
        center_corr 		(12, 8)
        marg_gap_shift 		(3, 4)
        offset 		        (3, 4)
        offset_corr 		(3, 4)
        quad_rotation 		(4,)
        quad_tilt 		    (4,)
        rotation 		    (4, 8)
        tilt 		        (4, 8)
        
    dir_name : str
        The name of the parent dir (parameter set) to place these parameters
        under
    """
    
    check_param_dict(param_dict)
    
    if os.path.exists(dir_name):
        print "WARNING: rm'ing %s..." % dir_name
        os.system('rm -r %s' % dir_name)
    
    os.mkdir(dir_name)
    
    for key in param_dict.keys():
        os.mkdir( pjoin( dir_name, key ))
        fname = pjoin( dir_name, key, '0-end.data' )
        if key == 'center' or key == 'center_corr':
            
            v = param_dict[key]
            v = v.reshape(3,4,8)
            
            f = open(fname, 'w')
            for i in range(3):
                for j in range(4):
                    f.write( '\t'.join([ str(x) for x in v[i,j,:] ]) + '\n' )
                f.write('\n')
            
        else:
            np.savetxt(fname, param_dict[key], fmt='%.2f')
        print "Wrote: %s" % fname
    
    return
    

def get_avg_from_hdf5(hd5_image_file, calibration_path, dsname, start_event,
                      n_events_to_avg=5):
    """
    Extract an average image from a psana hdf5 file.    
    """

    print 'Getting raw CSPad event %d from file %s \ndataset %s' % (event, fname, dsname)
    ds1ev = hm.getAverageCSPadEvent( fname, dsname, start_event, nevents=n_events_to_avg )
    if not ds1ev.shape == (32, 185, 388):
        print 'WARNING: ds1ev.shape =', ds1ev.shape, "should be (32, 185, 388)"

    return ds1ev


def get_event_from_npz(npz_image_file):
    return np.load(npz_image_file)['arr_0']
    
    
def plot_assembled_image(image):
    
    plt.imshow(image.T)
    plt.show()
        
    return


def test1():
    """ test assembly from disk """
    
    # change the paths below to run this script
    calibration_path = 'example_calibration_dir' # Mikhail's hand optimized geom
    run_number = 58
    
    raw_image = get_event_from_npz('npz_examples/cxi64813_r58_evt1.npz')
    assembled_image = assemble_image_from_dir(raw_image, calibration_path, run_number)
    plot_assembled_image(assembled_image)
    
    return
    
    
def test2():
    """ test assembly from memory """
    
    param_dict = {}
    raw_image = get_event_from_npz('npz_examples/cxi64813_r58_evt1.npz')
    
    for key in essential_params.keys():
        param_dict[key] = np.ones(essential_params[key])
    
    ai = assemble_image_from_params(raw_image, param_dict)
    plot_assembled_image(ai)
    
    return
    
    
def test3():
    """ test parameter write to disk """
    
    param_dict = {}
    raw_image = get_event_from_npz('npz_examples/cxi64813_r58_evt1.npz')
    
    for key in essential_params.keys():
        if key == 'rotation_index':
            param_dict[key] = np.ones(essential_params[key], dtype=np.int)
        elif key == 'quad_rotation':
            param_dict[key] = np.array([180,   90,    0,   270])
        else:
            param_dict[key] = np.ones(essential_params[key])
    
    write_params_to_dir(param_dict, 'dumb_test_params')
    ai = assemble_image_from_dir(raw_image, 'dumb_test_params')
    plot_assembled_image(ai)
    
    return
    

if __name__ == "__main__":
    #test1()
    #test2()
    test3()
