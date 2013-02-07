
import os

import numpy as np
from scipy import optimize, interpolate
from scipy.ndimage import filters
from scipy.stats.mstats import gmean

import matplotlib.pyplot as plt
import matplotlib.patches as plt_patches



def load_example_data(fname='sibeh_image.npz'):
    path = os.path.join('..', 'test_data', fname)
    image = np.load(path)['arr_0']
    return image
    

def compute_radii(image, center):
    """
    Center is (x,y) tuple in pixel units.
    """
    
    x = np.arange(image.shape[0])
    y = np.arange(image.shape[1])
    
    XX, YY = np.meshgrid(y, x)
    
    dx = np.power( XX - center[0], 2 )
    dy = np.power( YY - center[1], 2 )
    r = np.sqrt( dx + dy )
    
    assert r.shape == image.shape
    
    return r
    
    
def bin_intensities_by_radius(image, center, n_bins=None):

    r = compute_radii(image, center)
    
    if n_bins == None:
        n_bins = max(image.shape) # / 2
    
    # assume we've got a binary filter applied for now (!)
    bin_values, bin_edges = np.histogram( r * image, bins=n_bins )
    bin_values = bin_values[1:]
    bin_centers = bin_edges[1:-1] + np.abs(bin_edges[2] - bin_edges[1])
        
    bin_values = smooth(bin_values)
    
    return bin_centers, bin_values
    
    
def find_edges(image, threshold=0.01):
    """
    applies an edge filter
    """
    
    image = np.abs(filters.sobel(image, 0)) + np.abs(filters.sobel(image, 1))
    image -= image.min()
    
    assert image.min() == 0
    assert image.max() > 0
    
    print 'threshold value: %d' % (image.max() * threshold)
    image = (image > (image.max() * threshold)).astype(np.int8)
            
    return image

    
def smooth(x, beta=10.0, window_size=None):
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
    
    if window_size == None:
        window_size = 20 #len(x) / 20
    
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
    
    
def fwhm(vector):
    """
    Compute the full-width half-max of "vector", where "vector" is a single
    peak.
    """

    vector -= vector.min()
    x = np.arange(len(vector))
    spline = interpolate.UnivariateSpline(x, vector-np.max(vector)/2, s=3)

    roots = spline.roots() # find the roots
    if not len(roots) == 2:
        raise RuntimeError('Could not find FWHM of `vector`')
    r1, r2 = roots
    
    return r1, r2
    
    
def extract_peak(vector, m_ind=None):
    """
    Gets the biggest peak.
    """
    
    if m_ind == None:
        m_ind = np.argmax(vector)
    
    if m_ind == 0:
        left_index = 0
    else:
        left_index = None
        i = m_ind
        while not left_index:
            diff = vector[i-1] - vector[i]
            if (diff >= 0) or (i-1 == 0):
                left_index = i
            i -= 1
        
    if m_ind == len(vector):
        right_index = len(vector)
    else:
        right_index = None            
        i = m_ind
        while not right_index:
            diff = vector[i+1] - vector[i]
            if (diff >= 0) or (i+1 == len(vector)):
                right_index = i
            i += 1
        
    assert left_index < right_index
                
    return vector[left_index:right_index+1]
    
    
def maxima_indices(a):
    maxima = np.where(np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True] == True)[0]
    return maxima
    
    
def optimize_center(image, use_edge_mask=True, alpha=1e-2):
    
    if use_edge_mask:
        image = find_edges(image)
    
    center_guess = (image.shape[0] / 2., image.shape[1] / 2.)
    
    def objective(c):
        bin_centers, bin_values = bin_intensities_by_radius(image, c)
        peak = extract_peak(bin_values)
        try:
            r,l = fwhm(peak)
            v = l-r
        except RuntimeError as e:
            v = 1e300
        v += bin_values.max() * alpha
        return v
        
    print "Initial guess:", center_guess
    print "Initial FWHM:", objective(center_guess)
    
    # guess some bounds
    buff = 10.0 # pixel units
    xb = (center_guess[0] - buff, center_guess[0] + buff)
    yb = (center_guess[1] - buff, center_guess[1] + buff)
    
    opt_center = optimize.fmin(objective, center_guess)
    print "opt center:", opt_center
    
    return opt_center
    
    
def plot_center(image, center, use_edge_mask=True):
    
    print "plotting..."
    
    if use_edge_mask:
        image = find_edges(image)
    
    fig = plt.figure(figsize=(15,6))
    
    ax = plt.subplot(121)
    ax.imshow(image.T)
    blob_circ = plt_patches.Circle(center, 15, fill=False, lw=2, ec='red')
    ax.add_patch(blob_circ)
    
    bin_centers, bin_values = bin_intensities_by_radius(image, center)
    ax = plt.subplot(122)
    ax.plot(bin_centers, bin_values, lw=2)
    ax.set_xlabel('Radius / Pixel Units')
    ax.set_ylabel('Radial Average')
    
    #plt.savefig('trial.pdf')
    plt.show()
    
    return
    
    
def main():

    for image_file in ['sibeh_image.npz', 'silver_sim.npz', 'ssrl_silver.npz']:
        image = load_example_data(image_file)
        opt_center = optimize_center(image, use_edge_mask=True)
        plot_center(image, opt_center)
    
    return
    
    
if __name__ == '__main__':
    main()
