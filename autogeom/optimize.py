
import numpy as np
from scipy import optimize, interpolate

from autogeom import cspad
from autogeom import utils

import matplotlib.pyplot as plt
import matplotlib.patches as plt_patches
    
class Optimizer(object):
    """
    Modify CSPad geometry parameters to optimize the geometry by some criterion.
    """
    
    def __init__(self, initial_cspad=None, params_to_optimize=None,
                 **kwargs):
        """
        Initialize an optimizer class.
        
        Parameters
        ----------
        initial_cspad : autogeom.cspad.CSpad
            A CSPad object containing the CSPad parameters to start from. If
            `None`, then use the default parameter values.
            
        params_to_optimize : list of str
            A list of the parameters to optimize. Can be any of
            
            ['pedestals', 'offset_corr', 'offset', 'common_mode', 
             'marg_gap_shift', 'rotation', 'center_corr', 'center', 'tilt',
             'filter', 'quad_rotation', 'pixel_status', 'quad_tilt']
             
            if `None`, defaults to:
                 
             ['offset_corr', 'marg_gap_shift']
             
            the only one you really might consider adding to this is:
            
              'center_corr'
            
            many of the options you don't want to touch -- they are measured
            optically to better precision than you can hope for!
            
            Note that the absolute center of the image (beam position) is also
            always found by optimization.
        """
                
        if initial_cspad:
            if type(initial_cspad) == dict:
                self.cspad = cspad.CSPad(initial_cspad)
            elif isinstance(initial_cspad, cspad.CSPad):
                self.cspad = initial_cspad
            else:
                raise TypeError('`initial_cspad` must be one of {dict, CSPad}')
        else:
            # todo should be psana defaults
            raise NotImplementedError()
        
        # parameters -- default values
        self.n_bins              = None
        self.peak_regulization   = 10.0
        self.threshold           = 4.585e-04
        self.minf_size           = 3
        self.medf_size           = 8
        self.horizontal_cut      = 0.1
        self.use_edge_filter     = True
        self.beta                = 10.0
        self.window_size         = 5
        self.pixel_size          = 0.109 # mm
        self.radius_range        = []
        self.beam_loc            = (900.0, 870.0)
        self.plot_each_iteration = True
        
        if params_to_optimize:
            self.params_to_optimize = params_to_optimize
        else:
            self.params_to_optimize = ['offset_corr'] #, 'marg_gap_shift']
        
        # parse kwargs into self -- this will replace the defaults above
        print ""
        for key in kwargs:
            if hasattr(self, key):
                self.__dict__[key] = kwargs[key]
                print "Set parameter : %s --> %s" % (key, str(kwargs[key]))
            else:
                raise ValueError('Invalid Parameter: %s' % key)
                
        # check radius_range is sane
        if not len(self.radius_range) % 2 == 0:
            raise ValueError('`radius_range`, which defines which regions of the'
                             ' radial projection to optimize, must contain an '
                             'even number of entries.')
        else:
            self.radius_range.sort()
        
        return

    
    def __call__(self, raw_image, return_maxima_locations=False):
        """
        Takes a raw_image and produces an optimized CSPad geometry.
        """
                                        
        self.optimize_geometry(raw_image)
        
        if return_maxima_locations:
            if self.use_edge_filter:
                raw_image = utils.find_rings(raw_image)
            else:
                raw_image = ( raw_image > self.threshold ).astype(np.bool)
            assembled_image = self.cspad(raw_image)
            bc, bv = self._bin_intensities_by_radius(self.beam_loc,
                                                     assembled_image)
            maxima_locations = self._maxima_indices(bv)
            return self.cspad, maxima_locations
            
        else:
            return self.cspad
    
        
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
    
    
    def _bin_intensities_by_radius(self, center, image, radii=None):
        """
        Bin pixel intensities by their radius.
        
        Parameters
        ----------
        image : np.ndarray, np.bool
            An image.
            
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
            The total intensity in the bin.
        """
        
        # if not image.dtype == np.bool:
        #     raise TypeError('`image` must be dtype np.bool')

        if radii == None:
            radii = self._compute_radii(center, image)
        else:
            if not radii.shape == image.shape:
                raise ValueError('`radii` and `image` must have same shape')
        
        if self.n_bins == None:
            n_bins = max(image.shape) / 2
        else:
            n_bins = self.n_bins
        
        if image.dtype == np.bool:
            bin_values, bin_edges = np.histogram( radii * image, bins=n_bins )
        else:
            bin_values, bin_edges = np.histogram( radii, weights=image, bins=n_bins )
            
        bin_values = bin_values[1:]
        bin_centers = bin_edges[1:-1] + np.abs(bin_edges[2] - bin_edges[1])
        
        bin_values = utils.smooth(bin_values, beta=self.beta, 
                                  window_size=self.window_size)
                                  
        # slice out only the requested parts of the radial profile
        if len(self.radius_range) > 0:
            include = np.zeros( len(bin_values), dtype=np.bool)
            for i in range( len(self.radius_range)/2 ):
                include += (bin_centers > self.radius_range[i]) * \
                           (bin_centers < self.radius_range[i+1])
            
            if np.sum(include) == 0:
                raise RuntimeError('`radius_range` slices were not big enough to '
                                   'inlcude any data!')
                
            bin_centers = bin_centers[include]
            bin_values  = bin_values[include]
        
        return bin_centers, bin_values
    
        
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
        a = utils.smooth(a)
        maxima = np.where(np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True] == True)[0]
        if maxima[-1] == len(a)-1:
            maxima = maxima[:-1]
        return maxima
    
    
    def _unravel_params(self, flat_param_array):
        """
        We have to pass the parameters to scipys optimization routines as a
        flat array. This function takes that flat array and returns a dict of
        arrays of the correct shape -- the shapes that cspad.CSPad expects.
        
        Returns
        -------
        param_dict : dict
            Dictionary mapping param_name --> param_value,
            with param_value an np.ndarray of the shape expected by CSPad
        """
        
        param_dict = {}

        start = 2 # the first two parameter values are reserved for beam_loc
        for p in self.params_to_optimize:
            param_arr_shape = cspad._array_sizes[p]
            num_params_expected = np.product( param_arr_shape )
            end = start + num_params_expected
            param_dict[p] = flat_param_array[start:end].reshape(param_arr_shape)
            start = end

        return param_dict
    
    
    def _objective(self, param_vals, raw_image):
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
        
        # un-ravel & inject the param values in the CSPad object
        param_dict = self._unravel_params(param_vals)
        self.cspad.set_many_params(param_dict.keys(), param_dict.values())
        
        assembled_image = self.cspad(raw_image)
        
        # TJL : since the edge filtering *should* be the same for all images, 
        # it would make sense to do it once early on and not every time we 
        # evaluate the obj fxn. However, for unknown reasons, the filters seem 
        # to work much better on the assembled image, so here it is...
        if self.use_edge_filter:
            assembled_image = utils.find_rings(assembled_image, 
                                               threshold=self.threshold, 
                                               minf_size=self.minf_size, 
                                               medf_size=self.medf_size)
        
        # the absolute center will always be the first two elements by convention
        self.beam_loc = param_vals[:2]
        
        bin_centers, bin_values = self._bin_intensities_by_radius(self.beam_loc, assembled_image)
        max_inds = self._maxima_indices(bin_values)
        
        # only count big maxima for regularization purposes
        n_maxima = 0
        for ind in max_inds:
            if bin_values[ind] > bin_values.mean() / 5.:
                n_maxima += 1
        
        #n_maxima = len(max_inds)
        
        # if plotting is requested, plot away!
        if self.plot_each_iteration:
            
            self._axL.cla()
            self._axR.cla()
            
            self._axL.imshow(assembled_image.T)
            self._axR.plot(bin_centers, bin_values, lw=2, color='k')
            
            blob_circ = plt_patches.Circle(self.beam_loc, 15, fill=False, lw=2, 
                                           ec='orange')
            self._axL.add_patch(blob_circ)
            
            if len(self.radius_range) > 0:
                for i in range( len(self.radius_range)/2 ):
                    self._axR.fill_between(np.arange(self.radius_range[i], 
                                                     self.radius_range[i+1]),
                                           self._axR.get_ylim()[0], 
                                           self._axR.get_ylim()[1], 
                                           facecolor='blue', alpha=0.5)
                
                self._axR.vlines(self.radius_range, self._axR.get_ylim()[0], 
                                 self._axR.get_ylim()[1], lw=2, color='k')
            plt.draw()
                  
        # --------- HERE IS THE OBJECTIVE FUNCTION -- MODIFY TO PLAY -----------
        # TJL: I will think about smart ways to inject other functions in here,
        # but it may also be good to keep this static, once we have something
        # nice and robust.
        
        obj = self._all_widths(bin_values) - 0.1 * np.mean(bin_values[max_inds]) + \
              self.peak_regulization * n_maxima
        
        # ----------------------------------------------------------------------
        
        print "objective value: %f, number of peaks: %d" % (obj, n_maxima)
        
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

        print ""
        print "Beginning optimization..."
        print "Optimizing:", self.params_to_optimize
        print "Objective function:  Sum[peak_widths] + %.2f * num_peaks" % self.peak_regulization
        print ""

        initial_guesses = np.concatenate([ self.cspad.get_param(p).flatten() \
                                           for p in self.params_to_optimize ])
                                           
        # add in the absolute center -- we need to optimize this as well!
        initial_guesses = np.concatenate([ self.beam_loc, initial_guesses ])

        # turn on interactive plotting -- this is the only way I've gotten it to work
        if self.plot_each_iteration:
            plt.ion()
            self._fig = plt.figure(figsize=(18,9))
            self._axL = self._fig.add_subplot(121)
            self._axR = self._fig.add_subplot(122)
            self._axR.set_xlabel('Radius')
            self._axR.set_ylabel('Intensity')

        # run simplex minimization
        opt_params = optimize.fmin_powell(self._objective, initial_guesses, 
                                   args=(raw_image,), xtol=1e-2, ftol=1e-2)
                                   
        # turn off interactive plotting
        if self.plot_each_iteration: plt.ioff()
                                   
        # un-ravel & inject the param values in the CSPad object
        param_dict = self._unravel_params(opt_params)
        self.cspad.set_many_params(param_dict.keys(), param_dict.values())

        return
    