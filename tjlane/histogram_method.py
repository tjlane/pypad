
import numpy as np
from scipy import optimize
from scipy.ndimage import filters

import matplotlib.pyplot as plt
import matplotlib.patches as plt_patches



def load_example_data():
    image = np.load('sibeh_image.npz')['arr_0']
    return image
    
    
def load_example2():
    fileID = open('p2sec_p3d_11_1_00010.bin', 'r')
    fileID.seek(1024)
    data = np.fromfile(fileID, dtype=np.float32)
    data2D = np.reshape(data,(2527, 2463))
    data2D = data2D * (data2D > 0.)
    return data2D
    

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
    
    
def bin_intensities_by_radius(image, center, n_bins=1000):

    r = compute_radii(image, center)
    
    # assume we've got a binary filter applied for now (!)
    bin_values, bin_edges = np.histogram( r*image, bins=n_bins )
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
    
    
def fwhm(vector):
    
    # create a spline of x and blue-np.max(blue)/2 
    
    spline = UnivariateSpline(x, vector-np.max(vector)/2, s=0)
    r1, r2 = spline.roots() # find the roots

    
def optimize_center(image, use_edge_mask=True):
    
    if use_edge_mask:
        image = find_edges(image)
    
    center_guess = (image.shape[0] / 2., image.shape[1] / 2.)
    print "Initial guess:", center_guess
    
    # guess some bounds
    buff = 10.0 # pixel units
    xb = (center_guess[0] - buff, center_guess[0] + buff)
    yb = (center_guess[1] - buff, center_guess[1] + buff)
    
    def objective(c):
        bin_centers, bin_values = bin_intensities_by_radius(image, c)
        return bin_values.max()
    
    opt_center = optimize.fmin(objective, center_guess)
    
    print "opt center:", opt_center
    
    return opt_center
    
    
def plot_center(image, center):
    
    print "plotting..."
    
    fig = plt.figure()
    plt.imshow(image.T)
    blob_circ = plt_patches.Circle(center, 15, fill=False, lw=2, ec='red')
    plt.gca().add_patch(blob_circ)
    plt.show()
    
    return
    
    
def main():

    image = load_example_data()
    center = (image.shape[0] / 2., image.shape[1] / 2.)
    
    image = find_edges(image)
    print "Smoothed"

    bin_centers, bin_values = bin_intensities_by_radius(image, center)
    print "initial plot"
    plt.figure()
    plt.plot(bin_centers, bin_values, lw=2)
    plt.show()
    
    opt_center = optimize_center(image)
    plot_center(image, opt_center)
    
    bin_centers, bin_values = bin_intensities_by_radius(image, opt_center)
    plt.figure()
    plt.plot(bin_centers, bin_values, lw=2)
    plt.show()
    
    return
    
    
if __name__ == '__main__':
    main()
