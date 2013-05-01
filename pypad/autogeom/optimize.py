
import numpy as np
from scipy import optimize, interpolate

from pypad import cspad
from pypad import utils
from pypad import plot

import matplotlib.pyplot as plt
import matplotlib.patches as plt_patches
    
class Optimizer(object):
    """
    Modify CSPad geometry parameters to optimize the geometry by some criterion.
    """
    
    def __init__(self, cspad_instance, **kwargs):
        """
        Initialize an optimizer class.
        
        Parameters
        ----------
        
            
        params_to_optimize : list of str
            
        """
        
        if not isinstance(cspad_instance, cspad.CSPad):
            raise TypeError('`cspad_instance` must be a pypad.cspad.CSPad object')
        self.cspad = cspad_instance
        
        # -------- optimization algorithm parameters -- default values ------- #
        # below are the defaults for each parameter used in the optimization   #
        # kwargs is parsed for values that match these parameters, and if      #
        # provided those parameter values are substituted here                 #
        
        # optimization / objective function
        self.objective_type      = 'overlap'
        self.params_to_optimize  = ['quad_offset','quad_rotation']
        
        # for objective_type : peak_height
        self.peak_weight         = 0.0
        self.width_weight        = 10.0
        
        # filtering / image processing
        self.use_edge_filter     = True
        self.threshold           = 4.5e-04
        self.sigma               = 1.0
        self.minf_size           = 3
        self.rank_size           = 8
        self.sobel               = False
        
        # smoothing / peak counting
        self.horizontal_cut      = 0.25 # defines the start % for finding peaks
        self.beta                = 10.0
        self.window_size         = 3

        # misc but important
        self.radius_range        = []
        self.plot_each_iteration = True
        self.dilation            = 5.0
        
        # probably deprecated
        self.n_bins              = 200
        
        # -------------------------------------------------------------------- #
        
        # parse kwargs into self -- this will replace the defaults above
        print ""
        for key in kwargs:
            if hasattr(self, key):
                self.__dict__[key] = kwargs[key]
                print "Set parameter : %s \t--> %s" % (key, str(kwargs[key]))
            else:
                raise ValueError('Invalid Parameter: %s' % key)
                
        # --- now, do some checking to make sure the more complicated values
        #     that got passed are ok
                
        # check radius_range is sane
        if self.radius_range in [None, 'None', False]:
            self.radius_range = None
        elif not len(self.radius_range) % 2 == 0:
            raise ValueError('`radius_range`, which defines which regions of the'
                             ' radial projection to optimize, must contain an '
                             'even number of entries.')
        else:
            self.radius_range = np.sort(np.array(self.radius_range, dtype=np.float))
        
        # apply the dilation
        self.cspad.dilate(self.dilation)
        
        return

    
    def __call__(self, raw_image, return_maxima_locations=False):
        """
        Takes a raw_image and produces an optimized CSPad geometry.
        """
                                        
        self.optimize_geometry(raw_image)
        
        if return_maxima_locations:
            if self.use_edge_filter:
                raw_image = utils.preprocess_image(raw_image,
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

        rmin = np.min( self.radius_range )
        rmax = np.max( self.radius_range )
        if (not np.any(bin_centers > rmax)) or (not np.any(bin_centers < rmin)):
            raise ValueError('Invalid radius range -- out of bounds of measured radii.')
        
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

        start = 0 # got rid of center opt. -- now it defines origin
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
        
        # compute the radial profile
        bc, bv = self.cspad.intensity_profile(raw_image, n_bins=self.n_bins)
        if self.radius_range == None:
            bin_centers, bin_values = bc, bv
        else:
            bin_centers, bin_values = self._slice(bc, bv)
        
        # if plotting is requested, plot away!
        if self.plot_each_iteration:
            self._axL.cla()
            self._axR.cla()
            plot.sketch_2x1s(self.cspad.pixel_positions, self._axL)
            self._axR.plot(bin_centers, bin_values, lw=2, color='k')
            self._axR.set_xlabel('Radius')
            self._axR.set_ylabel('Intensity')
            plt.draw()
            
        # --------- HERE IS THE OBJECTIVE FUNCTION -- MODIFY TO PLAY -----------

        # new objective function : overlap integral
        if self.objective_type == 'overlap':
            
            bins = bin_centers
            quad_profiles = []

            for i in range(4):
                bc, bv = self.cspad.intensity_profile(raw_image,
                                                      n_bins=bins,
                                                      quad=i)
                quad_profiles.append(bv)
                
            m = np.min([ a.shape for a in quad_profiles ])
            quad_profiles = np.array([ a[:m] for a in quad_profiles ])
            
            obj = - np.sum(np.product(quad_profiles, axis=0))
                                       
        # old objective function : peak height
        elif self.objective_type == 'peak_height':
            obj = - bin_values.max()
            
            # # only count big maxima for regularization purposes
            # max_inds = self._maxima_indices(bin_values)
            # n_maxima = 0
            # for ind in max_inds:
            #     if bin_values[ind] > bin_values.mean() / 2.:
            #         n_maxima += 1
        
            if self.width_weight != 0.0:
                obj += self.width_weight * self._simple_width(bin_values, bar=self.horizontal_cut)
            
            if self.peak_weight != 0.0:
                obj += self.peak_weight
        
        else:
            raise ValueError('No implemented objective_type: %s' % objective_type)
        
        # ----------------------------------------------------------------------
        
        print "objective value: %.4e" % obj
        
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

        # turn on interactive plotting -- this is the only way I've gotten it 
        # to work
        if self.plot_each_iteration:            
            plt.ion()
            self._fig = plt.figure(figsize=(12,6))
            self._axL = self._fig.add_subplot(121)
            self._axL.set_aspect('equal')
            self._axR = self._fig.add_subplot(122)

        if self.use_edge_filter:
            image = utils.preprocess_image(raw_image, threshold=self.threshold,
                                           sigma=self.sigma,
                                           minf_size=self.minf_size,
                                           rank_size=self.rank_size,
                                           sobel=self.sobel)
        else:
            image = raw_image

        # run minimization -- downhill simplex
        opt_params = optimize.fmin(self._objective, initial_guesses, 
                                   args=(image,), xtol=1e-3, ftol=1e-3)
                                          
        # turn off interactive plotting
        if self.plot_each_iteration: plt.ioff()
                                   
        # un-ravel & inject the param values in the CSPad object
        param_dict = self._unravel_params(opt_params)
        self.cspad.set_many_params(param_dict.keys(), param_dict.values())

        return
    
