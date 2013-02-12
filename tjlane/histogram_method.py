
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
    
    
class AssembledImage(object):
    
    def __init__(self, images, center_guesses=None, **kwargs):
        """
        init
        """
        
        self.images = images
        
        # parameters -- default values
        self.n_bins             = None
        self.peak_regulization  = 10.0
        self.threshold          = 0.01
        self.minf_size          = 3
        self.medf_size          = 10
        self.horizontal_cut     = 0.1
        self.use_edge_mask      = True
        self.beta               = 10.0
        self.window_size        = 10
        
        
        # parse kwargs into self
        for key in kwargs:
            if hasattr(self, key):
                self.key = kwargs[key]
            else:
                raise ValueError('Invalid Parameter: %s' % key)
        
        
        # decide what to do
        if type(self.images) == np.ndarray:
            self._single_image = True
            
            if center_guesses:
                if not type(center_guesses) == tuple:
                    raise TypeError('`center_guesses` should be tuple if passing'
                                    ' only one image.')
            
            self._fragment_center = self.find_center(self.images, 
                                              center_guess=center_guesses,
                                              use_edge_mask=self.use_edge_mask)
            
        elif type(images) == list:
            self._single_image = False
            self.reconstruct_fragmented(images, center_guesses=center_guesses)
            
            

                
        else:
            raise TypeError('`images` must be one of {np.ndarray, list}')
            
            
    @property
    def fragment_centers(self):
        if self._single_image:
            return self._fragment_center
        else:
            raise NotImplementedError()
        return
    
        
    @property
    def fragment_corners(self):
        return
    
        
    @property
    def assembled_image(self):
        if self._single_image:
            assembled_image = self.images
        else:
            assembled_image = self._assemble_image(self.images)
        return assembled_image
    
    
    @property
    def basis_vector_representation(self):
        return
    
    @property
    def xyz_representation(self):
        return
        
        
    def _compute_radii(self, center, image):
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
    
    
    def _bin_intensities_by_radius(self, center, binary_image, radii=None):
        """
        Bin binary pixel intensities by their radius.
        
        Parameters
        ----------
        binary_image : np.ndarray, np.bool
            A binary image.
            
        center : tuple
            (x,y) in pixel units.
            
        Optional Parameters
        -------------------
        radii : ndarray
            
        Returns
        -------
        bin_centers : ndarray, float
            The radial center of each bin.
        
        bin_values : ndarray, int
            The number of pixels of value "1" in the bin.
        """
        
        if not binary_image.dtype == np.bool:
            raise TypeError('`binary_image` must be dtype np.bool')

        if radii == None:
            radii = self._compute_radii(center, binary_image)
        else:
            if not radii.shape == binary_image.shape:
                raise ValueError('`radii` and `binary_image` must have same shape')
    
        if self.n_bins == None:
            n_bins = max(binary_image.shape) / 2
        else:
            n_bins = self.n_bins
    
        # assume we've got a binary filter applied for now (!)
        bin_values, bin_edges = np.histogram( radii * binary_image, bins=n_bins )
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
        
        
    def _assemble_image(self, images, centers):
        """
        Take `images` and 
        """
        
        raise NotImplementedError()
        
        points = self.detector.xyz[:,:2]

        x = np.linspace(points[:,0].min(), points[:,0].max(), num_x)
        y = np.linspace(points[:,1].min(), points[:,1].max(), num_y)
        grid_x, grid_y = np.meshgrid(x,y)

        grid_z = interpolate.griddata(points, self.intensities, (grid_x,grid_y), 
                                      method='nearest')
        
        return grid_z
    
    
    def _objective(self, center, image):
        """
        The objective function for finding a good center. Minimize this.
        """

        bin_centers, bin_values = self._bin_intensities_by_radius(center, image)
        n_maxima = len(self._maxima_indices(bin_values))
        obj = self._all_widths(bin_values) + self.peak_regulization * n_maxima
        
        logger.debug("objective value: %f" % obj)
        
        return obj
    
        
    def _multi_image_objective(self, centers, images):
        """
        Assemble `images` using `centers`, then evaluate the usual objective
        function (self._objective) on the assembled image.
        
        Parameters
        ----------
        centers : ndarray, float
            An array with the coordinates of each centers.
        """
        
        assert len(centers) == len(images) * 2
        n_images = len(images)
        
        # compute the radii of each image wrt to its center
        flat_radii  = []
        flat_images = []
        for i in range(n_images):
            flat_images.append( images[i].flatten() )
            flat_radii.append( self._compute_radii(centers[i*2:i*2+2], images[i]).flatten() )
            
        # transform lists into arrays
        flat_radii  = np.concatenate(flat_radii)
        flat_images = np.concatenate(flat_images)
        assert flat_radii.shape == flat_images.shape
        assert flat_images.dtype == np.bool
            
        # bin & evaluate the objective function
        bin_centers, bin_values = self._bin_intensities_by_radius((0.0, 0.0), 
                                                flat_images, radii=flat_radii)
                                                    
        n_maxima = len(self._maxima_indices(bin_values))
        obj = self._all_widths(bin_values) + self.peak_regulization * n_maxima
        
        logger.debug("multi image objective value: %f" % obj)
        
        return obj
        
    
    def find_center(self, image, center_guess=None, use_edge_mask=True):
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
    
        if not center_guess:
            center_guess = (image.shape[0] / 2., image.shape[1] / 2.)
        
        logger.debug("Initial guess for center: %s" % str(center_guess))
        logger.debug("Initial FWHM: %f" % self._objective(center_guess, image))
        
        # run simplex minimization
        opt_center = optimize.fmin(self._objective, center_guess, args=(image,))
    
        logger.debug("optimzed center: %s" % str(opt_center))
    
        return opt_center
    

    def reconstruct_fragmented(self, images, center_guesses=None):
        """
        Reconstruct a fragmented image of rings (e.g. CSPAD quads) into a 
        unified image.
        """
        
        if not center_guesses:
            raise NotImplementedError()
            
        binary_images = []
        for image in images:
            binary_images.append( self._find_edges(image) )

        # first, get a guess by finding each center individually
        initial_centers = []
        for i,bimage in enumerate(binary_images):
            center = self.find_center(bimage, center_guess=center_guesses[i],
                                      use_edge_mask=self.use_edge_mask)
            initial_centers.extend( center )
        initial_centers = np.array(initial_centers)
            
        # next, do a global optimization of the geometry as a function of
        # all the centers
        opt_center = optimize.fmin(self._multi_image_objective, 
                                   initial_centers, args=(binary_images,))
        

        self._fragment_centers = []

        return
        

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
        
        
def load_example_data(fname='sibeh_image.npz', subdir='assembled'):
    path = os.path.join('..', 'test_data', subdir, fname)
    image = np.load(path)['arr_0']
    return image
        

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
    

def test_assembled():

    for image_file in ['silver_sim.npz', 'cxi0112_image_r096_ev1.npz', 'ssrl_silver.npz', 'sibeh_image.npz']:
        image = load_example_data(image_file, subdir='assembled')
        ai = AssembledImage(image)
        ai.plot_assembled()
    
    return
    
    
def test_quads():
    
    silver_files = ['silver_sim_bottomleft.npz',  'silver_sim_topleft.npz',
                    'silver_sim_bottomright.npz', 'silver_sim_topright.npz']
                    
                    
    sibeh_files = ['sibeh_bottomleft.npz', 'sibeh_topleft.npz', 
                   'sibeh_bottomright.npz', 'sibeh_topright.npz']
                   
    center_guesses = [ (1,99), (98,99), (0,3), (99,1) ]
        
        
    images = []
    for i,f in enumerate(silver_files):
        images.append(load_example_data(f, subdir='quads'))
        
    ai = AssembledImage(images, center_guesses=center_guesses)
    ai.plot_assembled()
                   
    
    return
    
    
if __name__ == '__main__':
    test_quads()
