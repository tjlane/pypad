

import h5py
import numpy as np


def load_raw_image(filename, image_in_file=0):
    """
    Attempts to be a general and forgiving file-loading platform for raw images.
    
    Currently supported formats:
        -- psana hdf5
        -- cheetah hdf5
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
    
    print "Loading: %s" % filename
    
    if filename.endswith('.h5'):
        f = h5py.File(filename)
        try:
            # psana format
            raw_image = np.array( f[('/data%d/raw' % image_in_file)] )
        except:
            # cheetah format
            raw_image = np.array( f['/data/data'] )
        finally:
            f.close()
        
    elif filename.endswith('.npz'):
        raw_image = np.load(filename)['arr_%d' % image_in_file]
        
    elif filename.endswith('txt'):
        raw_image = np.loadtxt(filename)
        
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
    
    else:
        raise ValueError("Cannot understand `raw_image`: does not have any"
                         " known dimension structure")
    
    return new_image