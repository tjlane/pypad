
import os, sys

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')

import numpy as np
from scipy import optimize, interpolate
from scipy.ndimage import filters

import matplotlib.pyplot as plt
import matplotlib.patches as plt_patches

from odin.math import smooth

def load_example_data(fname='sibeh_image.npz'):
    path = os.path.join('..', 'test_data', fname)
    image = np.load(path)['arr_0']
    return image
    
    
class AssembledImage(object):
    
    def __init__(self, images):
        """
        init
        """
        
        self.images = images
        
        # parameters
        self.n_bins = None
        self.peak_regulization = 10.0
        self.threshold = 0.01
        self.minf_size = 3
        self.medf_size = 10
        self.horizontal_cut = 0.1
        self.use_edge_mask = True
        self.beta = 10.0
        self.window_size = 10
        
        
        # decide what to do
        if type(self.images) == np.ndarray:
            self._single_image = True
            self._pixel_center = self.find_center(self.images, 
                                               use_edge_mask=self.use_edge_mask)
            
        elif type(images) == list:
            self._single_image = False
            self._pixel_center = []
            
            for image in images:
                center = self.find_center(image, use_edge_mask=self.use_edge_mask)
                self._pixel_center.append( center )
                
        else:
            raise TypeError('`images` must be one of {np.ndarray, list}')
            
            
    @property
    def pixel_center(self):
        if self._single_image:
            return self._pixel_center
        else:
            raise NotImplementedError()
        return
    
        
    @property
    def image_corners(self):
        return
    
        
    @property
    def image_angles(self):
        return
    
        
    @property
    def assembled_image(self):
        if self._single_image:
            assembled_image = self.images
        else:
            pass
        return assembled_image
    
    
    @property
    def basis_vector_representation(self):
        return
    
    @property
    def xyz_representation(self):
        return
        
        
    def _compute_radii(self, image, center):
        """
        Compute the radii of each pixel, in pixel units, with respect to 
        `center`.
        
        Parameters
        ----------
        center : tuple
            (x,y) in pixel units.
            
        image : ndarray
            The image
            
        Returns
        -------
        r : ndarray
            The radius of each pixel.
        """
    
        x = np.arange(image.shape[0])
        y = np.arange(image.shape[1])
    
        XX, YY = np.meshgrid(y, x)
    
        dx = np.power( XX - center[0], 2 )
        dy = np.power( YY - center[1], 2 )
        r = np.sqrt( dx + dy )
    
        assert r.shape == image.shape
    
        return r
    
    
    def _bin_intensities_by_radius(self, center, binary_image):
        """
        Bin binary pixel intensities by their radius.
        
        Parameters
        ----------
        binary_image : np.ndarray, np.bool
            A binary image.
            
        center : tuple
            (x,y) in pixel units.
            
        Returns
        -------
        bin_centers : ndarray, float
            The radial center of each bin.
        
        bin_values : ndarray, int
            The number of pixels of value "1" in the bin.
        """
        
        if not binary_image.dtype == np.bool:
            raise TypeError('`binary_image` must be dtype np.bool')

        r = self._compute_radii(binary_image, center)
    
        if self.n_bins == None:
            n_bins = max(binary_image.shape) / 2
        else:
            n_bins = self.n_bins
    
        # assume we've got a binary filter applied for now (!)
        bin_values, bin_edges = np.histogram( r * binary_image, bins=n_bins )
        bin_values = bin_values[1:]
        bin_centers = bin_edges[1:-1] + np.abs(bin_edges[2] - bin_edges[1])
        
        bin_values = smooth(bin_values, beta=self.beta, 
                            window_size=self.window_size)
    
        return bin_centers, bin_values
    
    
    def _find_edges(self, image):
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
    
        image = np.abs(filters.sobel(image, 0)) + np.abs(filters.sobel(image, 1))
        image -= image.min()
    
        assert image.min() == 0
        assert image.max() > 0
    
        logger.debug('threshold value: %d' % (image.max() * self.threshold))
        image = (image > (image.max() * self.threshold)).astype(np.bool)
    
        image = filters.minimum_filter(image, size=self.minf_size)
        image = filters.median_filter(image, size=self.medf_size)
            
        return image.astype(np.bool)
    
    
    def _all_widths(self, vector, bar=None):
        """
        Compute the sum of all peak widths.
        
        Parameters
        ----------
        vector : ndarray
            A vector describing a plot with peaks.
            
        Returns
        -------
        width_sum : float
            The sum of all the widths of the peaks in `vector`.
        """
        
        if bar == None:
            bar = self.horizontal_cut # rename for conv.
        
        x = np.arange(len(vector))
        spline = interpolate.UnivariateSpline(x, vector - vector.max()*bar, s=3)
        roots = spline.roots()
    
        if len(roots) % 2 == 0:
            width_sum = np.sum(roots[1::2] - roots[::2])
        elif len(roots) == 0:
            raise RuntimeError('Width finder failed -- odd number of roots at '
                               'all horizontal cuts.')
        else:
            newbar = bar + (1.-bar) / 8.
            width_sum = self._all_widths(vector, newbar)
    
        return width_sum
    
    
    
    def _maxima_indices(self, a):
        """
        Return the indicies of all local maxima in `a`.
        """
        a = smooth(a)
        maxima = np.where(np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True] == True)[0]
        if maxima[-1] == len(a)-1:
            maxima = maxima[:-1]
        return maxima
    
    
    def _objective(self, center, image):
        """
        The objective function for finding a good center. Minimize this.
        """

        bin_centers, bin_values = self._bin_intensities_by_radius(center, image)
        n_maxima = len(self._maxima_indices(bin_values))
        obj = self._all_widths(bin_values) + self.peak_regulization * n_maxima
    
        logger.debug("objective value: %f" % obj)
    
        return obj
    
    
    def find_center(self, image, use_edge_mask=True):
        """
        Find the center of `image`, which contains one or more concentric rings.
        
        Parameters
        ----------
        image : ndarray
            The image.
            
        use_edge_mask : bool
            Whether or not to apply an edge mask noise filter before finding
            the center. Highly recommended.
            
        Returns
        -------
        center : tuple
            The (x,y) position of the center, in pixel units
        """
    
        if use_edge_mask:
            image = self._find_edges(image)
    
        center_guess = (image.shape[0] / 2., image.shape[1] / 2.)
        
        logger.debug("Initial guess for center: %s" % str(center_guess))
        logger.debug("Initial FWHM: %f" % self._objective(center_guess, image))
        
        # run simplex minimization
        opt_center = optimize.fmin(self._objective, center_guess, args=(image,))
    
        logger.debug("optimzed center: %s" % str(opt_center))
    
        return opt_center
    

    def plot_assembled(self, use_edge_mask=True, fig_name=None):
    
        image = self.assembled_image
        
        logger.info("Plotting assembled image...")
    
        if use_edge_mask:
            image = self._find_edges(image)
    
        fig = plt.figure(figsize=(15,6))
    
        ax = plt.subplot(121)
        ax.imshow(image.T)
        blob_circ = plt_patches.Circle(self.pixel_center, 15, fill=False, lw=2, ec='orange')
        ax.add_patch(blob_circ)
    
        bin_centers, bin_values = self._bin_intensities_by_radius(self.pixel_center, image)
        ax = plt.subplot(122)
        ax.plot(bin_centers, bin_values, lw=2)
        ax.set_xlabel('Radius / Pixel Units')
        ax.set_ylabel('Radial Average')
    
        if fig_name:
            plt.savefig(fig_name)
            logger.info("Saved: %s" % fig_name)
        plt.show()
    
        return
        

def plot_objective_surface(image):

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
    
def main():

    for image_file in ['silver_sim.npz', 'cxi0112_image_r096_ev1.npz', 'ssrl_silver.npz', 'sibeh_image.npz']:
        image = load_example_data(image_file)
        ai = AssembledImage(image)
        ai.plot_assembled()
    
    return
    
    
if __name__ == '__main__':
    main()
