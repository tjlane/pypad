
"""
utils.py

Various utility functions. Filtering, image processing, etc.
"""

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
