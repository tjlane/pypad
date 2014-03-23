
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
utils.py

Various utility functions. Filtering, image processing, etc.
"""

import numpy as np
from scipy.ndimage import filters
from scipy import interpolate

import matplotlib.pyplot as plt

from pypad import read


def update_average(n, A, B):
    '''
    updates a numpy matrix A that represents an average over the previous n-1 shots
    by including B into the average, B being the nth shot
    '''
    if n == 0:
        A += B
        return
    else:
        A *= (n-1)/float(n)
        A += (1.0/float(n))*B
        return


def normalize(q_values, intensities, q_min=0.5, q_max=3.5):
    assert q_values.shape == intensities.shape
    inds = (q_values > q_min) * (q_values < q_max)
    factor = float(np.sum(intensities[inds])) / float(np.sum(inds))
    return intensities / factor


class RadialAverager(object):
    
    def __init__(self, q_values, mask, n_bins=101):
        """
        Parameters
        ----------
        q_values : np.ndarray (float)
            For each pixel, this is the momentum transfer value of that pixel
        mask : np.ndarray (int)
            A boolean (int) saying if each pixel is masked or not
        n_bins : int
            The number of bins to employ. If `None` guesses a good value.
        """
        
        self.q_values = read.enforce_raw_img_shape(q_values)
        self.mask = read.enforce_raw_img_shape(mask).astype(np.int)
        self.n_bins = n_bins
        
        # figure out the number of bins to use
        if n_bins != None:
            self.n_bins = n_bins
            self._bin_factor = float(self.n_bins-0.5) / self.q_values.max()
        else:
            self._bin_factor = 25.0
            self.n_bins = (self.q_values.max() * self._bin_factor) + 1
        
        self._bin_assignments = np.floor( q_values * self._bin_factor ).astype(np.int32)
        self._normalization_array = (np.bincount( self._bin_assignments.flatten(),
                                     weights=self.mask.flatten() ) \
                                     + 1e-100).astype(np.float)

        #print self.n_bins, self._bin_assignments.max() + 1 
        assert self.n_bins == self._bin_assignments.max() + 1, 'bin mismatch in init'
        self._normalization_array = self._normalization_array[:self.n_bins]
        
        return
    
    def __call__(self, image):
        """
        Bin pixel intensities by their momentum transfer.
        
        Parameters
        ----------            
        image : np.ndarray
            The intensity at each pixel, same shape as pixel_pos


        Returns
        -------
        bin_centers : ndarray, float
            The q center of each bin.

        bin_values : ndarray, int
            The average intensity in the bin.
        """

        image = read.enforce_raw_img_shape(image)
        
        if not (image.shape == self.q_values.shape):
            raise ValueError('`image` and `q_values` must have the same shape')
        if not (image.shape == self.mask.shape):
            raise ValueError('`image` and `mask` must have the same shape')

        weights = image.flatten() * self.mask.flatten()
        bin_values = np.bincount(self._bin_assignments.flatten(), weights=weights)
        bin_values /= self._normalization_array
   
        assert bin_values.shape[0] == self.n_bins, 'bin number mismatch (%d, %d)' \
                                                   % (bin_values.shape[0], self.n_bins)
    
        return bin_values

    @property
    def bin_centers(self):
        return np.arange(self.n_bins) / self._bin_factor


def arctan3(y, x):
    """
    Compute the inverse tangent. Like arctan2, but returns a value in [0,2pi].
    """
    theta = np.arctan2(y,x)
    if type(theta) == np.ndarray:
        theta[theta < 0.0] += 2 * np.pi
    else:
        if theta < 0.0: theta += 2 * np.pi
    return theta


def preprocess_image(raw_image, threshold=0.0025, sigma=1.0, minf_size=1,
                     rank_size=1, sobel=True):
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
    if raw_image.shape == (4,16,185,194):
        non_flat_img = True
        image = np.zeros((1480, 1552), dtype=np.float) # flat image
        for i in range(8):
            for j in range(4):
                x_start = 185 * i
                x_stop  = 185 * (i+1)
                y_start = 388 * j
                y_stop  = 388 * (j+1)
                
                two_by_one = np.hstack(( raw_image[j,i*2,:,:], raw_image[j,i*2+1,:,:])).astype(np.float)
                image[x_start:x_stop,y_start:y_stop] = two_by_one
                
    elif len(raw_image.shape) == 2:
        non_flat_img = False
        image = raw_image.astype(np.float)
        
    else:
        raise ValueError('`raw_image` should be 2d or shape-(4,16,185,194), got'
                         ': %s' % str(raw_image.shape))
    
    # apply rank filter & gaussian filter
    if rank_size > 2:
        image = filters.rank_filter(image, -1, size=rank_size)
    if sigma > 0.1:
        image = filters.gaussian_filter(image, sigma=sigma)
    
    image -= image.min()
    assert image.min() == 0
    assert image.max() > 0
    
    # threshold
    image = (image > (image.max() * threshold))
    
    if minf_size > 2:
        image = filters.minimum_filter(image, size=minf_size)
    if sobel:
        image = np.abs(filters.sobel(image, 0)) + np.abs(filters.sobel(image, 1))    
    
    if non_flat_img:
        image = read.enforce_raw_img_shape( image.astype(np.bool) )
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


def multi_for(iterables):
    """
    Some python magic for multi-dimensional for loop
    
    Parameters
    ----------
    iterables: a list of iterables
        e.g. [xrange(2), xrange(2), xrange(1), xrange(2)]
    
    Returns
    -------
    The multi-dimensional iterable
        e.g. (0, 0, 0, 0) ... (1, 1, 0, 1)
    """
    if not iterables:
        yield ()
    else:
        for item in iterables[0]:
            for rest_tuple in multi_for(iterables[1:]):
                yield (item,) + rest_tuple


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
        The average intensities corresponding to `bin_centers`.
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
    
    # normalize by number of pixels in each bin to obtain average intensity
    bin_normalizations = np.histogram( r, bins=bin_edges )
    
    bin_values = bin_values/bin_normalizations[0]
    bin_centers = np.array([(bin_edges[i] + bin_edges[i+1])/2 for i in range(len(bin_values))])
    
    return bin_centers, bin_values


def flatten_2x1s(image):
    """
    Takes a non-2D image : (32, 185, 388), (4, 8, 185, 388), (4, 16, 185, 194) 
    and returns a 2d version of that image (for visualization).
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
                
                
    elif image.shape == (4, 16, 185, 194):
        for i in range(4):
            for j in range(8):
                xs = 185 * j
                xe = 185 * (j + 1)
                ys = 388 * i
                ye = 388 * (i + 1)
                flat_image[xs:xe,ys:ye] = np.hstack((image[i,j*2,:,:], image[i,j*2+1,:,:]))
        
    else:
        raise ValueError('Invalid shape for arg `image`')
    
    return flat_image
    

def _assemble_implicit(xyz, raw_image, num_x=2000, num_y=2000):
    """
    For testing purposes ONLY. This is just a slower, shittier version of the
    assembly algorithm provided by cspad.CSPad that was used to ensure
    internal consistency.
    """
    
    if xyz.shape == (3, 4, 8, 185, 388):
        xyz = np.array(( xyz[0].flatten(),
                         xyz[1].flatten(),
                         xyz[2].flatten(), )).T
    
    assert xyz.shape[1] == 3
    assert raw_image.shape == (4, 8, 185, 388)
    
    points = xyz[:,:2] # ignore z-comp. of detector
    
    x = np.linspace(points[:,0].min(), points[:,0].max(), num_x)
    y = np.linspace(points[:,1].min(), points[:,1].max(), num_y)
    grid_x, grid_y = np.meshgrid(x,y)
    
    flat_image = raw_image.flatten()
    
    print "interpolating"
    grid_z = interpolate.griddata(points, flat_image, (grid_x,grid_y), 
                                  method='nearest', fill_value=0.0)
    
    return grid_z

    
def ER_rotation_matrix(axis, theta):
    """
    Compute the rotation matrix defining a rotation in 3D of angle `theta` 
    around `axis` using the Euler-Rodrigues formula.

    Parameters
    ----------
    axis : np.ndarray, float
        A 3-vector defining the axis of rotation.

    theta : float
        The rotation angle, in radians.

    Returns
    -------
    R : np.ndarray, float
        A 3x3 array defining a rotation matrix.

    Citation
    --------
    ..[1] https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_parameters
    """

    # note: the a,b,c,d parameters here don't correspond with their equivalents
    # elsewhere in the code

    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    R = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])
    return R
