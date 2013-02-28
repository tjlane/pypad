
"""
utils.py

Various utility functions. Filtering, image processing, etc.
"""

import tables

import numpy as np
from scipy.ndimage import filters
from scipy import interpolate


def find_rings(raw_image, threshold=0.0025, minf_size=1, medf_size=8, sobel=True):
    """
    Applies an edge filter followed by a noise reduction filter. Very good
    at locating powder rings and filtering everything else out.
    
    Parameters
    ----------
    raw_image : ndarray
        An image to find the edges of
        
    Returns
    -------
    binary_image : ndarray, np.bool
        A binary image, with "1" where there are powder rings/strong edges
    """
    
    # flatten the image into a two-D array and later re-process it
    # convert to cheetah-like format
    if raw_image.shape == (4, 8, 185, 388):
        non_flat_img = True
        image = np.zeros((1480, 1552), dtype=np.float) # flat image
        for i in range(8):
            for j in range(4):
                x_start = 185 * i
                x_stop  = 185 * (i+1)
                y_start = 388 * j
                y_stop  = 388 * (j+1)
                psind = i + j * 8
                image[x_start:x_stop,y_start:y_stop] = raw_image[j,i,:,:].astype(np.float)
    elif len(raw_image.shape) == 2:
        non_flat_img = False
        image = raw_image.astype(np.float)
    else:
        raise ValueError('`raw_image` should be 2d or shape-(4,8,185,388), got'
                         ': %s' % str(raw_image.shape))
    
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
    
    if non_flat_img:
        image = enforce_raw_img_shape( image.astype(np.bool) )
    else:
        image = image.astype(np.bool)
    
    return image
    
    
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
    
    raw_image = enforce_raw_img_shape(raw_image)
    
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
    

def enforce_raw_img_shape(raw_image):
    """
    Make sure that the `raw_image` has shape: (4,8,185,388).
    
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
    if raw_image.shape == (4,8,185,388):
        new_image = raw_image # we're good
        
    # PyCSPad format
    elif raw_image.shape == (32,185,388):
        new_image = np.zeros((4, 8, 185, 388), dtype=raw_image.dtype)
        for i in range(8):
            for j in range(4):
                psind = i + j * 8
                new_image[j,i,:,:] = raw_image[psind,:,:]
        
    # Cheetah format
    elif raw_image.shape == (1480, 1552):
        new_image = np.zeros((4, 8, 185, 388), dtype=raw_image.dtype)
        for i in range(8):
            for j in range(4):
                x_start = 185 * i
                x_stop  = 185 * (i+1)
                y_start = 388 * j
                y_stop  = 388 * (j+1)
                psind = i + j * 8
                new_image[j,i,:,:] = raw_image[x_start:x_stop,y_start:y_stop]
    
    else:
        raise ValueError("Cannot understand `raw_image`: does not have any"
                         " known dimension structure")
    
    return new_image


def _assemble_implicit(xyz, raw_image, num_x=2000, num_y=2000):
    """
    For testing purposes ONLY. This is just a slower, shittier version of the
    assembly algorithm provided by cspad.CSPad that was used to ensure
    internal consistency.
    """
    
    assert xyz.shape[1] == 3

    points = xyz[:,:2] # ignore z-comp. of detector
    
    x = np.linspace(points[:,0].min(), points[:,0].max(), num_x)
    y = np.linspace(points[:,1].min(), points[:,1].max(), num_y)
    grid_x, grid_y = np.meshgrid(x,y)
    
    # flatten out the raw image
    # assert raw_image.shape == (4,8,185,388)
    # flat_image = np.zeros( np.product(raw_image.shape) )
    # interval = 185 * 388
    # 
    # for i in range(4):
    #     for j in range(8):
    #         ind = i*8 + j
    #         flat_image[interval*ind:interval*(ind+1)] = raw_image[i,j,:,:].flatten()
    
    flat_image = raw_image.flatten()
    
    grid_z = interpolate.griddata(points, flat_image, (grid_x,grid_y), 
                                  method='linear', fill_value=0.0)
    
    return grid_z


def sketch_2x1s(pixel_positions, mpl_axes=None):
    """
    Draw a rough sketch of the layout of the CSPAD

    Parameters
    ----------
    pixel_positions : np.ndarray
        The x,y,z coordinates of the pixels on the CSPAD
    """
    
    quad_color = ['k', 'g', 'purple', 'b']

    if not mpl_axes:
        plt.figure()
        ax = plt.subplot(111)
    else:
        ax = mpl_axes

    for i in range(4):
        for j in range(8):
            x = pixel_positions[0,i,j,:,:]
            y = pixel_positions[1,i,j,:,:]
            corners = np.zeros((5,2))

            corners[0,:] = np.array([ x[0,0],   y[0,0] ])     # bottom left
            corners[1,:] = np.array([ x[0,-1],  y[0,-1] ])    # bottom right
            corners[3,:] = np.array([ x[-1,0],  y[-1,0] ])    # top left
            corners[2,:] = np.array([ x[-1,-1], y[-1,-1] ])   # top right
            corners[4,:] = np.array([ x[0,0],   y[0,0] ])     # make rectangle

            ax.plot(corners[:,0], corners[:,1], lw=2, color=quad_color[i])
            ax.scatter(x[0,0], y[0,0])

    if mpl_axes:
        return ax
    else:
        plt.show()

