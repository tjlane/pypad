
import os, sys

import numpy as np
from scipy import optimize, interpolate, misc
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
        window_size = len(x) / 10
    
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
    
    plt.figure()
    plt.scatter(x, vector)
    plt.plot(x, spline(x))
    plt.show()
    
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
                break
            i -= 1
        
    if m_ind+1 == len(vector):
        right_index = len(vector)
    else:
        right_index = None            
        i = m_ind
        while not right_index:
            diff = vector[i+1] - vector[i]
            if (diff >= 0) or (i+2 == len(vector)):
                right_index = i
                break
            i += 1
        
    assert left_index <= right_index
                
    return vector[left_index:right_index+1]
    
    
def peak_widths(vector):
    
    x = np.arange(len(vector))    
    spline = interpolate.UnivariateSpline(x, vector, s=3)
    maxima = maxima_indices( spline(x) )
    ddx = misc.derivative( spline, x, n=2 )
    inflections = np.where( np.abs(ddx) < 1e-2 )[0]
    
    peak_widths = []
    
    for m_ind in maxima:
        diffs = 1.0 / (m_ind - inflections)
        left  = inflections[np.argmax(diffs)]
        right = inflections[np.argmin(diffs)]
        assert right > left
        width = right - left
        peak_widths.append(width)
        
    peak_widths = np.array(peak_widths)
    
    plt.figure()
    plt.plot(x, vector)
    plt.vlines(maxima, vector.min(), vector.max(), color='red')
    plt.vlines(inflections, vector.min(), vector.max())
    plt.show()
        
    return maxima, peak_widths
    
    
def maxima_indices(a):
    a = smooth(a)
    maxima = np.where(np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True] == True)[0]
    if maxima[-1] == len(a)-1:
        maxima = maxima[:-1]
    return maxima
    
    
def optimize_center(image, use_edge_mask=True, alpha=0.1):
    
    if use_edge_mask:
        image = find_edges(image)
    
    center_guess = (image.shape[0] / 2., image.shape[1] / 2.)
    
    def objective(c):
        bin_centers, bin_values = bin_intensities_by_radius(image, c)
        
        # maxima = maxima_indices(bin_values)
        # obj_values = []
        
        # plt.figure()
        # plt.plot(bin_centers, bin_values, lw=2)
        # plt.vlines(bin_centers[maxima], bin_values.min(), bin_values.max())
        # plt.show()
        
        # for i,m_ind in enumerate(maxima):
        #     peak = extract_peak(bin_values, m_ind)
        #     try:
        #         r,l = fwhm(peak)
        #         v = (l-r) + bin_values[m_ind] * alpha
        #         obj_values.append(v)
        #     except:
        #         # if we fail, just don't include that one
        #         pass
        
        max_inds, widths = peak_widths(bin_values)
        obj = np.mean( (widths - alpha) * bin_values[max_inds] )
        
        print "obj:", obj, len(max_inds)
        
        return obj
        
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

    for image_file in ['silver_sim.npz', 'sibeh_image.npz', 'ssrl_silver.npz']:
        image = load_example_data(image_file)
        opt_center = optimize_center(image, use_edge_mask=True)
        plot_center(image, opt_center)
    
    return
    
    
if __name__ == '__main__':
    main()
