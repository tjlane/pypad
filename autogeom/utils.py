
"""
utils.py

Various utility functions. Filtering, image processing, etc.
"""

import tables

import numpy as np
from scipy.ndimage import filters


def find_rings(image, threshold=0.0025, minf_size=1, medf_size=8, sobel=True):
    """
    Applies an edge filter followed by a noise reduction filter. Very good
    at locating powder rings and filtering everything else out.
    
    Parameters
    ----------
    image : ndarray
        An image to find the edges of
        
    Returns
    -------
    binary_image : ndarray, np.bool
        A binary image, with "1" where there are powder rings/strong edges
    """
    
    image = image.astype(np.float64)
    image = np.abs(filters.sobel(image, 0)) + np.abs(filters.sobel(image, 1))
    
    # do a "common-mode"-esque normalization (seems to help)
    if image.shape == (32, 185, 388):
        for i in range(32):
            image[i,:,:] -= np.min(image[i,:,:])
    elif image.shape == (4, 8, 185, 388):
        for i in range(4):
            for j in range(8):
                image[i,j,:,:] -= np.min(image[i,j,:,:])
    else:
        image -= image.min()
    
    assert image.min() == 0
    assert image.max() > 0
    
    image = (image > (image.max() * threshold)).astype(np.bool)
    
    image = filters.minimum_filter(image, size=minf_size)
    image = filters.median_filter(image, size=medf_size)
    
    return image.astype(np.bool)
    
    
def smooth(x, beta=10.0, window_size=11):
    """
    Apply a Kaiser window smoothing convolution.
    
    Parameters
    ----------
    x : ndarray, float
        The array to smooth.
    
    Optional Parameters
    -------------------
    beta : float
        Parameter controlling the strength of the smoothing -- bigger beta 
        results in a smoother function.
    window_size : int
        The size of the Kaiser window to apply, i.e. the number of neighboring
        points used in the smoothing.
    
    Returns
    -------
    smoothed : ndarray, float
        A smoothed version of `x`.
    """
    
    # make sure the window size is odd
    if window_size % 2 == 0:
        window_size += 1
    
    # apply the smoothing function
    s = np.r_[x[window_size-1:0:-1], x, x[-1:-window_size:-1]]
    w = np.kaiser(window_size, beta)
    y = np.convolve( w/w.sum(), s, mode='valid' )
    
    # remove the extra array length convolve adds
    b = (window_size-1) / 2
    smoothed = y[b:len(y)-b]
    
    return smoothed


def radial_profile(image, center):
    """
    Compute the radial intensity profile of `image`.
    
    Parameters
    ----------
    image : np.ndarray
        The image to perform the radial profile on.
        
    center : tuple of floats
        The center to use, in pixel units
        
    Returns
    -------
    bin_centers : np.ndarray
        The radial position (x-coordinate)
    
    bin_values
        The intensities corresponding to `bin_centers`.
    """
    
    # compute the radii
    x = np.arange(image.shape[0])
    y = np.arange(image.shape[1])
    
    XX, YY = np.meshgrid(y, x)

    dx = np.power( XX - center[0], 2 )
    dy = np.power( YY - center[1], 2 )
    r = np.sqrt( dx + dy )

    assert r.shape == image.shape
    
    # histogram the intensities
    n_bins = max(image.shape) / 2
    
    if image.dtype == np.bool:
        bin_values, bin_edges = np.histogram( r * image, bins=n_bins )
    else:
        bin_values, bin_edges = np.histogram( r, weights=image, bins=n_bins )
    
    bin_values = bin_values[1:]
    bin_centers = bin_edges[1:-1] + np.abs(bin_edges[2] - bin_edges[1])
    
    return bin_centers, bin_values


def flatten_2x1s(image):
    """
    Takes a non-2D image : either (32, 185, 388) or (4, 8, 185, 388) and returns
    a 2d version of that image (for visualization).
    """
    
    flat_image = np.zeros((185*8, 388*4))
    
    if image.shape == (32, 185, 388):
        for i in range(4):
            for j in range(8):
                xs = 185 * j
                xe = 185 * (j + 1)
                ys = 388 * i
                ye = 388 * (i + 1)
                flat_image[xs:xe,ys:ye] = image[i*4+j,:,:]
        
    elif image.shape == (4, 8, 185, 388):
        for i in range(4):
            for j in range(8):
                xs = 185 * j
                xe = 185 * (j + 1)
                ys = 388 * i
                ye = 388 * (i + 1)
                flat_image[xs:xe,ys:ye] = image[i,j,:,:] 
        
    else:
        raise ValueError('Invalid shape for arg `image`')
    
    return flat_image
    
    
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
    
    if filename.endswith('.h5'):
        f = tables.File(filename)
        try:
            # psana format
            raw_image = f.getNode('/data%d/raw' % image_in_file).read()
        except:
            # cheetah format
            raw_image = f.root.data.data.read()
        finally:
            f.close()
        
    elif filename.endswith('.npz'):
        raw_image = np.load(filename)['arr_%d' % image_in_file]
        
    elif filename.endswith('txt'):
        raw_image = np.loadtxt(filename)
        
    else:
        raise ValueError('Cannot understand format of file: %s' % fn)
    
    return raw_image

    
def cheetah_to_3Dpsana(cheetah_image):
    """
    Takes a raw cheetah image (2D) and returns it in psana format (3D)
    """

    psana_image = np.zeros((32, 185, 388))
    assert cheetah_image.shape == (1480, 1552)

    for i in range(8):
        for j in range(4):
            x_start = 185 * i
            x_stop  = 185 * (i+1)
            y_start = 388 * j
            y_stop  = 388 * (j+1)
            psind = i + j * 8 # confirmed visually
            psana_image[psind,:,:] = cheetah_image[x_start:x_stop,y_start:y_stop]

    return psana_image
