
import os, sys
import tables

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

import assemble
    
    
class AssembledImage(object):
    
    def __init__(self, raw_image, initial_param_dict=None, params_to_optimize=None,
                 **kwargs):
        """
        init
        """
        
        self.raw_image = raw_image
        
        if initial_param_dict:
            assert type(initial_param_dict) == dict
            self.initial_param_dict = initial_param_dict
        else:
            # todo should be psana defaults
            raise NotImplementedError()
            
        # add an absolute center to the optimization no matter what
        self.initial_param_dict['abs_center'] = np.array([1200., 1200.])
        params_to_optimize.append('abs_center')
        
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
        self.params_to_optimize = assemble.essential_params.keys()
        
        # parse kwargs into self
        for key in kwargs:
            if hasattr(self, key):
                self.key = kwargs[key]
            else:
                raise ValueError('Invalid Parameter: %s' % key)
                
                
        # optimize the geometry
        self.opt_param_dict = self.optimize_geometry(raw_image)
        
        return

    
    @property
    def basis_vector_representation(self):
        raise NotImplementedError()
        return
    
    
    @property
    def xyz_representation(self):
        raise NotImplementedError()
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
    
    
    def _inject_params_into_dict(self, param_dict, list_of_params, param_values):
        """
        """
        
        start = 0
        for p in list_of_params:
            param_arr_shape = assemble.geometry_params[p]
            num_params_expected = np.product( param_arr_shape )
            end = start + num_params_expected
            param_dict[p] = param_values[start:end].reshape(param_arr_shape)
            start = end
        
        return param_dict
    
    
    def _objective(self, param_vals, list_of_params, raw_image):
        """
        The objective function for finding a good center. Minimize this.
        
        Parameters
        ----------
        param_vals : ndarray
            A flat array of the parameters we're optimizing. Which parameters
            these are depends on the entires in `list_of_params`.
            
        list_of_params : list of str
            A list of which parameters are contained in `params_to_opt`.
        
        Returns
        -------
        obj : float
            The value of the function. Lower is better.
            
        See Also
        --------
        self._inject_params_into_dict : function
            Makes sense of `params_to_opt` and `list_of_params`
        """
        
        param_dict_u = self._inject_params_into_dict(self.initial_param_dict, 
                                                     list_of_params, 
                                                     param_vals)
        
        assembled_image = assemble.assemble_image_from_params(raw_image, param_dict_u)
        
        plt.imshow(assembled_image.T)
        plt.show()
        assembled_image = assembled_image.astype(np.bool)
        plt.imshow(assembled_image.T)
        plt.show()
        
        abs_center = self.initial_param_dict['abs_center']
        bin_centers, bin_values = self._bin_intensities_by_radius(abs_center, assembled_image)
        n_maxima = len(self._maxima_indices(bin_values))
        
        # HERE IS THE OBJECTIVE FUNCTION -- MODIFY TO PLAY
        obj = self._all_widths(bin_values) + self.peak_regulization * n_maxima
        
        logger.debug("objective value: %f" % obj)
        
        return obj
    
        
    def optimize_geometry(self, raw_image):
        """
        Find the center of `image`, which contains one or more concentric rings.

        Parameters
        ----------
        raw_image : ndarray
            The image, in pyana/psana's raw format.

        use_edge_mask : bool
            Whether or not to apply an edge mask noise filter before finding
            the center. Highly recommended.

        Returns
        -------
        params_dict : dict
            Dict of the pyana parameters used to optimize the geometry.
        """

        if self.use_edge_mask:
            raw_image = self._find_edges(raw_image)

        print self.params_to_optimize    
        print self.initial_param_dict

        initial_guesses = np.concatenate([ self.initial_param_dict[p].flatten() \
                                           for p in self.params_to_optimize ])

        # run simplex minimization
        opt_params = optimize.fmin(self._objective, initial_guesses, 
                                   args=(self.params_to_optimize, raw_image))
        opt_param_dict = self._inject_params_into_dict(self.initial_param_dict, 
                                                       self.params_to_optimize, 
                                                       opt_params)

        return opt_param_dict
        
        
    def plot_assembled(self):
        """
        """
        
        i = assemble_image_from_params(self.raw_image, self.opt_param_dict)
        assemble.plot_assembled_image(i)
        
        return
    
        
def cheetah_to_psana(cheetah_image):
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
    
    
def load_AgBe():
    
    f = tables.File('hdf5_examples/AgBe/r0003-RawSum.h5')
    cheetah_agbe = f.root.data.data.read()
    psana_agbe = cheetah_to_psana(cheetah_agbe)
    
    return psana_agbe
    
    
    
def test_cheetah_conv():
    raw_image = load_AgBe()
    ai = assemble.assemble_image_from_dir(raw_image, 'example_calibration_dir')
    assemble.plot_assembled_image(ai)
    return
    
    
def test_agbe_assembly():
    
    params_to_opt = ['offset_corr']
    
    raw_image = load_AgBe()
    initial_param_dict = assemble.load_params_from_dir('example_calibration_dir')
    ai = AssembledImage(raw_image, initial_param_dict, params_to_opt)
    ai.plot_assembled()
    
    return
    
    
    
if __name__ == '__main__':
    test_agbe_assembly()
