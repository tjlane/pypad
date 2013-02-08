
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
        n_bins = max(image.shape) / 2
    
    # assume we've got a binary filter applied for now (!)
    bin_values, bin_edges = np.histogram( r * image, bins=n_bins )
    bin_values = bin_values[1:]
    bin_centers = bin_edges[1:-1] + np.abs(bin_edges[2] - bin_edges[1])
        
    bin_values = smooth(bin_values)
    
    return bin_centers, bin_values
    
    
def find_edges(image, threshold=0.01, minf_size=3, medf_size=10):
    """
    applies an edge filter
    """
    
    image = np.abs(filters.sobel(image, 0)) + np.abs(filters.sobel(image, 1))
    image -= image.min()
    
    assert image.min() == 0
    assert image.max() > 0
    
    print 'threshold value: %d' % (image.max() * threshold)
    image = (image > (image.max() * threshold)).astype(np.bool)
    
    image = filters.minimum_filter(image, size=minf_size)
    image = filters.median_filter(image, size=medf_size)
            
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
        window_size = len(x) / 20
    
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
    
    
def peak_widths(vector, spacing=1e6):
    
    x = np.arange(len(vector)) * ( spacing / float(len(vector)) )
    spline = interpolate.UnivariateSpline(x, vector, s=3)
    maxima = maxima_indices( spline(x) )
    
    dx  = misc.derivative( spline, x, n=1 )
    ddx = misc.derivative( spline, x, n=2 )
    inflections = np.where( (np.abs(ddx) < 5e-7) * (np.abs(dx) > 1e-3) )[0]
    assert len(inflections) > 1
    
    peak_widths = []
    maxima_location = []
    
    for m_ind in maxima:
        diffs = 1.0 / (m_ind - inflections)
        left  = inflections[np.argmax(diffs)]
        right = inflections[np.argmin(diffs)]
        if right > left:
            width = right - left
            peak_widths.append(width)
            maxima_location.append(m_ind)
        
    peak_widths = np.array(peak_widths)
    maxima_location = np.array(maxima_location, dtype=np.int)
    
    # plt.figure()
    # plt.plot(x, vector)
    # plt.vlines(x[maxima], vector.min(), vector.max(), color='red')
    # plt.vlines(x[inflections], vector.min(), vector.max())
    # plt.show()
    # sys.exit()
    
        
    return maxima_location, peak_widths
    
    
def peak_widths_v2(vector, window=20):
    
    x = np.arange(len(vector))
    x2 = np.arange(window*2) 
    spline = interpolate.UnivariateSpline(x, vector, s=3)
    maxima = maxima_indices( spline(x) )
    maxima = maxima[ vector[maxima] > vector.mean() ]
    
    peak_widths = []
    real_maxima = []
    
    for m_ind in maxima:
        start = max([0, m_ind-window])
        end   = min([len(vector)-1, m_ind+window])
        slc = vector[start:end]
        
        spline = interpolate.UnivariateSpline( np.arange(end-start), slc-slc.mean() )
        
        roots = spline.roots()
        if not len(roots) == 2:
            pass
        else:
            real_maxima.append(m_ind)
            peak_widths.append( roots[1] - roots[0] )

    peak_widths = np.array(peak_widths)
    assert len(peak_widths) == len(real_maxima)
    
    # plt.figure()
    # plt.plot(x, vector)
    # plt.vlines(x[maxima], vector.min(), vector.max(), color='red')
    # plt.show()

    return real_maxima, peak_widths

    
def peak_widths_v3(vector):

    x = np.arange(len(vector))
    spline = interpolate.UnivariateSpline(x, vector, s=3)
    
    maxima = maxima_indices( spline(x) )
    ddx = misc.derivative( spline, x, n=2 )
    
    peak_widths = ddx[maxima]
    
    # plt.figure()
    # plt.plot(x, vector)
    # plt.plot(x, ddx * 30)
    # plt.vlines(maxima, 0.0, peak_widths*10, lw=2, color='red')
    # plt.show()
    # sys.exit()
    
    return maxima, peak_widths
    
    
def all_widths(vector, bar=0.1):
    
    x = np.arange(len(vector))
    spline = interpolate.UnivariateSpline(x, vector - vector.max()*bar, s=3)
    roots = spline.roots()
    
    if len(roots) % 2 == 0:
        width_sum = np.mean(roots[1::2] - roots[::2])
    elif len(roots) == 0:
        raise RuntimeError('Width finder failed.')
    else:
        newbar = bar + (1.-bar) / 8.
        width_sum = all_widths(vector, newbar)
    
    return width_sum
    
    
    
def maxima_indices(a):
    a = smooth(a)
    maxima = np.where(np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True] == True)[0]
    if maxima[-1] == len(a)-1:
        maxima = maxima[:-1]
    return maxima
    
    
def objective(center, image, alpha=0.5, beta=100.0):

    bin_centers, bin_values = bin_intensities_by_radius(image, center)
    obj = all_widths(bin_values)
    
    print "obj:", obj
    
    return obj
    

def objective_help(args):
    return objective(*args)
    
    
def optimize_center(image, use_edge_mask=True):
    
    if use_edge_mask:
        image = find_edges(image)
    
    center_guess = (image.shape[0] / 2., image.shape[1] / 2.)
        
    print "Initial guess:", center_guess
    print "Initial FWHM:", objective(center_guess, image)
    
    # guess some bounds
    buff = 10.0 # pixel units
    xb = (center_guess[0] - buff, center_guess[0] + buff)
    yb = (center_guess[1] - buff, center_guess[1] + buff)
    
    opt_center = optimize.fmin(objective, center_guess, args=(image,))
    
    # opt_center = optimize.fmin_bfgs(objective, center_guess, 
    #                                 args=(image,))
    
    
    print "opt center:", opt_center
    
    return opt_center
    
    
def plot_surface(image):
    
    image = find_edges(image)
    center_guess = (image.shape[0] / 2., image.shape[1] / 2.)
    
    buff = 5 # pixel units
    xb = (center_guess[0] - buff, center_guess[0] + buff)
    yb = (center_guess[1] - buff, center_guess[1] + buff)
    
    Xb = np.arange(*xb)
    Yb = np.arange(*yb)
    
    z = np.zeros(( len(Xb), len(Yb) ))
    for ix,x in enumerate(Xb):
        for iy,y in enumerate(Yb):
            z[ix, iy] = objective((x,y), image)
            
    plt.figure()
    plt.imshow(z.T, interpolation='nearest')
    plt.show()
    
    wopt = np.where( z == z.min() )
    cx = wopt[0] + center_guess[0] - buff
    cy = wopt[1] + center_guess[1] - buff
    
    return (cx, cy)
    
    
    
def plot_center(image, center, use_edge_mask=True):
    
    print "plotting..."
    
    if use_edge_mask:
        image = find_edges(image)
    
    fig = plt.figure(figsize=(15,6))
    
    ax = plt.subplot(121)
    ax.imshow(image.T)
    blob_circ = plt_patches.Circle(center, 15, fill=False, lw=2, ec='orange')
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

    for image_file in ['silver_sim.npz', 'ssrl_silver.npz', 'sibeh_image.npz']:
        image = load_example_data(image_file)
        
        #opt_center = plot_surface(image)
        opt_center = optimize_center(image, use_edge_mask=True)
        
        plot_center(image, opt_center)
    
    return
    
    
if __name__ == '__main__':
    main()
