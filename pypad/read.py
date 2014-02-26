# THIS FILE IS PART OF PyPad, AND IS GOVERENED BY A PERMISSIBILITY LICENSE 
# GOVERNING ITS USE AND DISTRIBUTION. YOU SHOULD HAVE RECIEVED A COPY OF THIS
# LICENSE WITH THE SOFTWARE; IF NOT PROVIDED, WRITE TO <tjlane@stanford.edu>.
#
# AUTHORS:
# TJ Lane <tjlane@stanford.edu>
# Jonas Sellberg <sellberg@slac.stanford.edu>
#
# Apr 30, 2013

"""
read.py

Functions for reading image files from disk, from a variety of formats.
"""

import os
import h5py
import numpy as np
from scipy import io as spio


def load_raw_image(filename, image_in_file=0):
    """
    Attempts to be a general and forgiving file-loading platform for raw images.
    
    Currently supported formats:
        -- psana hdf5
        -- cheetah hdf5 (old & CXIdb)
        -- npz : numpy-z compression
        -- txt : flat text
    
    Parameters
    ----------
    filename : str
        The file to load.
        
    image_in_file : int
        The image contained in that file to load.
        
    Returns
    -------
    image : np.ndarray
        A numpy array of the image.
    """
    
    if not os.path.exists(filename):
        raise IOError('File: %s does not exist' % filename)
    
    print "Loading: %s" % filename
    
    if filename.endswith('.h5'):
        f = h5py.File(filename, 'r')

        # psana format
        if ('data%d' % image_in_file) in f.keys():
            try:
                raw_image = np.array( f[('/data%d/raw' % image_in_file)] )
            except:
                raise IOError('Was expecting psana format, but: /dataX/raw not found!')
            
        # old cheetah format
        elif 'data' in f.keys():
            try:
                raw_image = np.array( f['/data/data'] )
            except:
                raise IOError('Was expecting old cheetah format, but: /data/data not found!')
            if raw_image.shape != (1480,1552):
                try:
                    raw_image = np.array( f['/data/correcteddata'] )
                except:
                    raise IOError('/data/correcteddata not found!')
        else:
            raise IOError('Could not safely interpret hdf5 file. Can only read'
                          'Cheetah and Psana h5 images.')
        
        f.close()

    elif filename.endswith('.cxi'):
        f = h5py.File(filename, 'r')
        try:
            ds = f['/entry_1/instrument_1/detector_1/data']
            raw_image = ds[image_in_file]
        except:
             raise IOError('Error reading image %d in cheetah/CXIdb file' % image_in_file)
        f.close()
        
    elif filename.endswith('.npz'):
        raw_image = np.load(filename)['arr_%d' % image_in_file]
        
    elif filename.endswith('txt'):
        raw_image = np.loadtxt(filename)
        
    elif filename.endswith('.shot'):
        try:
            from odin import xray
        except ImportError as e:
            raise ImportError('To read `.shot` files, you must first install '
                              'odin: https://github.com/tjlane/odin')
        ss = xray.Shotset.load(filename)
        raw_image = ss.average_intensity
        
    elif filename.endswith('.mat'):
        mat = spio.loadmat(filename)
        if not 'img_avg' in mat.keys():
            raise IOError('Could not find `img_avg` key in .mat file: %s!' % filename)
        else:
            raw_image = mat['img_avg']
        
    else:
        raise ValueError('Cannot understand format of file: %s' % filename)
    
    raw_image = enforce_raw_img_shape(raw_image)
    
    return raw_image

    
def enforce_raw_img_shape(raw_image):
    """
    Make sure that the `raw_image` has shape: (4,16,185,194).
    
    Which is (quad, 2x1, fast, slow). This function will attempt to get
    the image into that form, and if it can't throw an error.
    
    Parameters
    ----------
    raw_image : np.ndarray
        The raw image, in XtcExplorer, PyCSPad, or Cheetah format
        
    Returns
    -------
    raw_image : np.ndarray
        The same image reshaped to be (4,8,185,388)
    """
    
    # XtcExporter format
    if raw_image.shape == (4,16,185,194):
        new_image = raw_image # we're good
        
    # PyCSPad format
    elif raw_image.shape == (32,185,388):
        #print "Converting: psana image --> pypad image"
        new_image = np.zeros((4,16,185,194), dtype=raw_image.dtype)
        for i in range(8):
            for j in range(4):
                psind = i + j * 8
                
                sec1, sec2 = np.hsplit(raw_image[psind,:,:], 2)
                new_image[j,i*2,:,:]   = sec1
                new_image[j,i*2+1,:,:] = sec2
        
    # Cheetah format
    # raw data format: 1480 rows x 1552 cols,
    # origin is lower left corner in doc/cspad_arrangement.pdf
    elif raw_image.shape == (1480, 1552):
        #print "Converting: cheetah image --> pypad image"
        
        new_image = np.zeros((4,16,185,194), dtype=raw_image.dtype)
        
        for q in range(4):
            for twoXone in range(8):
                
                x_start = 388 * q
                x_stop  = 388 * (q+1)
                
                y_start = 185 * twoXone
                y_stop  = 185 * (twoXone + 1)
                
                sec1, sec2 = np.hsplit(raw_image[y_start:y_stop,x_start:x_stop], 2)
                
                new_image[q,twoXone*2,:,:]   = sec1
                new_image[q,twoXone*2+1,:,:] = sec2
                
                
    # ODIN format
    elif raw_image.shape == (2296960,):
        
        new_image = np.zeros((4,16,185,194), dtype=raw_image.dtype)
        spacing = 185 * 194
        
        for i in range(4):
            for j in range(16):
                ASIC = i * 16 + j
                new_image[i,j,:,:] = raw_image[spacing*ASIC:spacing*(ASIC+1)].reshape(185, 194)
            
    
    else:
        raise ValueError("Cannot understand `raw_image`: does not have any"
                         " known convertable shape. Image shape: %s" % \
                         str(raw_image.shape))
    
    return new_image

