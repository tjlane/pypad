
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
    
    def __init__(self, geometry=None, params_to_optimize=['offset_corr_xy'],
                 **kwargs):
        """
        Initialize an optimizer class.
        
        Parameters
        ----------
        geometry : autogeom.cspad.CSpad OR autogeom.cspad.Metrology
            An initial specification of the CSPad geometry, via either a set
            of psana alignment parameters (CSPad object) or an optical metrology
            (Metrology object). If `None`, defaults to a (probably bad) geometry.
            
        params_to_optimize : list of str
            A list of the parameters to optimize. If `geometry` is a CSPad 
            object, can be any of
            
            ['pedestals', 'offset_corr', 'offset', 'common_mode', 
             'marg_gap_shift', 'rotation', 'center_corr', 'center', 'tilt',
             'filter', 'quad_rotation', 'pixel_status', 'quad_tilt']
             
            defaults to:
                 
             ['offset_corr_xy'], a custom one
             
            the only one you really might consider adding to this is:
            
              'center_corr'
            
            many of the options you don't want to touch -- they are measured
            optically to better precision than you can hope for!
            
            If `geometry` is a Metrology object, then the only things updated
            will be the quad positions and beam_location.
            
            Note that the absolute center of the image (beam position) is also
            always found by optimization.
        """
        
        # handle all possible `geometry` cases
        if isinstance(geometry, cspad.CSPad):
            self.cspad = geometry
            
        elif isinstance(geometry, cspad.Metrology):
            if params_to_optimize != ['offset_corr_xy']:
                raise ValueError('Can only optimize `offset_corr_xy` with '
                                 'Metrology object as geometry')
            self.cspad = geometry
            
        elif type(initial_cspad) == dict:
            self.cspad = cspad.CSPad(geometry)
            
        elif geometry == None:
            self.cspad = cspad.CSPad.default()
            
        else:
            raise TypeError('`geometry` must be one of {None, dict, CSPad, Metrology}')
        
        # parameters -- default values
        self.objective_type      = 'overlap'
        self.n_bins              = 1000
        self.peak_weight         = 0.0
        self.width_weight        = 10.0
        self.threshold           = 4.5e-04
        self.sigma               = 1.0
        self.minf_size           = 3
        self.rank_size           = 8
        self.sobel               = False
        self.horizontal_cut      = 0.25
        self.use_edge_filter     = True
        self.beta                = 10.0
        self.window_size         = 3
        self.pixel_size          = 0.10992 # mm
        self.radius_range        = []
        self.beam_location       = np.array((900.0, 870.0))
        self.plot_each_iteration = True
        
        self.params_to_optimize = params_to_optimize
        
        
        # parse kwargs into self -- this will replace the defaults above
        print ""
        for key in kwargs:
            if hasattr(self, key):
                self.__dict__[key] = kwargs[key]
                print "Set parameter : %s \t--> %s" % (key, str(kwargs[key]))
            else:
                raise ValueError('Invalid Parameter: %s' % key)
                
        # check radius_range is sane
        if not len(self.radius_range) % 2 == 0:
            raise ValueError('`radius_range`, which defines which regions of the'
                             ' radial projection to optimize, must contain an '
                             'even number of entries.')
        else:
            self.radius_range = np.sort(np.array(self.radius_range, dtype=np.float))
        
        if not len(self.beam_location) == 2:
            raise ValueError('`beam_location` must be len 2')
        self.beam_location = np.array(self.beam_location)
        
        return

    
    def __call__(self, raw_image, return_maxima_locations=False):
        """
        Takes a raw_image and produces an optimized CSPad geometry.
        """
                                        
        self.optimize_geometry(raw_image)
        
        if return_maxima_locations:
            if self.use_edge_filter:
                raw_image = utils.find_rings(raw_image,
                                             threshold=self.threshold,
                                             sigma=self.sigma,
                                             minf_size=self.minf_size,
                                             rank_size=self.rank_size,
                                             sobel=self.sobel)
            assembled_image = self.cspad(raw_image)
            bc, bv = self._bin_intensities_by_radius(self.beam_location,
                                                     assembled_image)
            maxima_locations = self._maxima_indices(bv)
            return self.cspad, maxima_locations
            
        else:
            return self.cspad
    
        
    def _slice(self, bin_centers, bin_values):
        """
        slice out only the requested parts of the radial profile
        """
        
        if len(self.radius_range) > 0:
            include = np.zeros( len(bin_values), dtype=np.bool)
            for i in range( len(self.radius_range)/2 ):
                include += (bin_centers > self.radius_range[i]) * \
                           (bin_centers < self.radius_range[i+1])

            if np.sum(include) == 0:
                raise RuntimeError('`radius_range` slices were not big enough to '
                                   'include any data!')

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
        
        # measure the distance between the roots
        if len(roots) % 2 == 0:
            width_sum = np.sum(roots[1::2] - roots[::2])
            
        # if we missed all the peaks return a big number to let optimizer know 
        # this is a bad parameter set
        elif len(roots) < 2:
            print ('Warning: width finder failed -- odd number of roots at '
                   'all horizontal cuts.')
            width_sum = 1e10 
            
        # other wise move the bar up and try again
        else:
            newbar = bar + (1.-bar) / 4.
            width_sum = self._all_widths(vector, newbar)
        
        return width_sum
    
    
    def _simple_width(self, vector, bar=0.25):
        m = bar * vector.max()
        width = np.sum((vector > m))
        return width
    
        
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

        start = 2 # the first two positions are reserved for beam_location
        #start = 0 # got rid of center opt. -- now it defines origin
        for p in self.params_to_optimize:
            param_arr_shape = cspad._array_sizes[p]
            num_params_expected = np.product( param_arr_shape )
            end = start + num_params_expected
            if end > len(flat_param_array):
                raise RuntimeError('array overrun on `flat_param_array`')
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
        
        # un-ravel & inject the param values in the geometry object
        param_dict = self._unravel_params(param_vals)
        self.cspad.set_many_params(param_dict.keys(), param_dict.values())
        
                
        # the absolute center will always be the first two elements by convention
        # note this means we're actually overspecifying the number of free
        # parameters (x/y for each quad + center), but this seems to optimize
        # to a better geometry, since the center acts as a kind of global move
        
        self.beam_location = param_vals[:2]
        self.cspad.set_param('beam_location', param_vals[:2])
        
        pp = self.cspad.pixel_positions / self.pixel_size
        bc, bv = self.cspad.intensity_profile(raw_image, n_bins=self.n_bins,
                                   beta=self.beta, window_size=self.window_size)
        bin_centers, bin_values = self._slice(bc / self.cspad.pixel_size, bv)
                                                                  
        max_inds = self._maxima_indices(bin_values)
        
        # only count big maxima for regularization purposes
        n_maxima = 0
        for ind in max_inds:
            if bin_values[ind] > bin_values.mean() / 2.:
                n_maxima += 1
        
        # if plotting is requested, plot away!
        if self.plot_each_iteration:
            self._axL.cla()
            self._axR.cla()
            utils.sketch_2x1s(pp, self._axL)
            self._axR.plot(bin_centers, bin_values, lw=2, color='k')
            blob_circ = plt_patches.Circle(self.beam_location, 15, fill=False, lw=2, 
                                           ec='orange')
            self._axL.add_patch(blob_circ)
            self._axR.set_xlabel('Radius')
            self._axR.set_ylabel('Intensity')
            plt.draw()
                  
        # --------- HERE IS THE OBJECTIVE FUNCTION -- MODIFY TO PLAY -----------

        # new objective function : overlap integral
        if self.objective_type == 'overlap':
            
            bins = bin_centers.copy() * self.cspad.pixel_size         
            quad_profiles = np.zeros((4, len(bins)-2))

            for i in range(4):
                bc, bv = self.cspad.intensity_profile(raw_image,
                                                      n_bins=bins,
                                                      beta=self.beta, 
                                                      window_size=self.window_size,
                                                      quad=i)
                quad_profiles[i,:] = bv + 1e-100
            
            obj = - np.log10( np.sum(np.product(quad_profiles, axis=0)) )
                                       
        # old objective function : peak height
        elif self.objective_type == 'peak_height':
            obj = - bin_values.max()
        
            if self.width_weight != 0.0:
                obj += self.width_weight * self._simple_width(bin_values, 
                                                            bar=self.horizontal_cut)
            
            if self.peak_weight != 0.0:
                obj += self.peak_weight * n_maxima
        
        else:
            raise ValueError('No implemented objective_type: %s' % objective_type)
        
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
        print "Objective function:  %s" % self.objective_type
        print ""

        initial_guesses = np.concatenate([ self.cspad.get_param(p).flatten() \
                                           for p in self.params_to_optimize ])
                                           
        # add in the absolute center -- we need to optimize this as well!
        initial_guesses = np.concatenate([ self.beam_location, initial_guesses ])

        # turn on interactive plotting -- this is the only way I've gotten it to work
        if self.plot_each_iteration:            
            plt.ion()
            self._fig = plt.figure(figsize=(12,6))
            self._axL = self._fig.add_subplot(121)
            self._axR = self._fig.add_subplot(122)

        if self.use_edge_filter:
            image = utils.find_rings(raw_image, threshold=self.threshold,
                                     sigma=self.sigma,
                                     minf_size=self.minf_size,
                                     rank_size=self.rank_size,
                                     sobel=self.sobel)

        # run minimization
        opt_params = optimize.fmin_powell(self._objective, initial_guesses, 
                                          args=(image,), xtol=1e-2, ftol=1e-3)
                                   
        # turn off interactive plotting
        if self.plot_each_iteration: plt.ioff()
                                   
        # un-ravel & inject the param values in the CSPad object
        param_dict = self._unravel_params(opt_params)
        self.cspad.set_many_params(param_dict.keys(), param_dict.values())

        return
    