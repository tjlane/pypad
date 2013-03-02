
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
            will be the quad positions and beam_loc.
            
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
        self.n_bins              = None
        self.peak_weight         = 10.0
        self.max_weight          = 0.1
        self.threshold           = 4.5e-04
        self.minf_size           = 3
        self.medf_size           = 8
        self.horizontal_cut      = 0.25
        self.use_edge_filter     = True
        self.beta                = 10.0
        self.window_size         = 10
        self.pixel_size          = 0.10992 # mm
        self.radius_range        = []
        self.beam_loc            = np.array((900.0, 870.0))
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
            #self.radius_range *= self.pixel_size # convert to mm
        
        if not len(self.beam_loc) == 2:
            raise ValueError('`beam_loc` must be len 2')
        self.beam_loc = np.array(self.beam_loc)
        #self.beam_loc *= self.pixel_size
        
        return

    
    def __call__(self, raw_image, return_maxima_locations=False):
        """
        Takes a raw_image and produces an optimized CSPad geometry.
        """
                                        
        self.optimize_geometry(raw_image)
        
        if return_maxima_locations:
            if self.use_edge_filter:
                raw_image = utils.find_rings(raw_image)
            assembled_image = self.cspad(raw_image)
            bc, bv = self._bin_intensities_by_radius(self.beam_loc,
                                                     assembled_image)
            maxima_locations = self._maxima_indices(bv)
            return self.cspad, maxima_locations
            
        else:
            return self.cspad
    
        
    def _bin_intensities_by_radius(self, center, pixel_pos, intensities):
        """
        Bin pixel intensities by their radius.
        
        Parameters
        ----------            
        center : tuple
            (x,y) in pixel units.
            
        pixel_pos : np.ndarray
            The x,y,z positions of each pixel
            
        intensities : np.ndarray
            The intensity at each pixel, same shape as pixel_pos
            
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

        if not pixel_pos.shape[1:] == intensities.shape:
            raise ValueError('`pixel_pos` and `intensities` must have same'
                             ' shape format. Current shapes: %s and %s respectively.' %\
                              (str(pixel_pos.shape), str(intensities.shape)))
            
        if not ((pixel_pos.shape[0] == 2) or (pixel_pos.shape[0] == 3)):
            raise ValueError('`pixel_pos` must have first dimension be 2 or 3.'
                             'Current shape: %s' % str(pixel_pos.shape))
            
        # use only x,y for now
        if (pixel_pos.shape[0] == 3):
            pixel_pos = pixel_pos[:2]
            
        # compute radii
        radii = np.sqrt( np.power(pixel_pos[0]-center[0], 2) + \
                         np.power(pixel_pos[1]-center[1], 2) )
        
        # generate the histogram
        if self.n_bins == None:
            n_bins = np.sqrt( np.product(intensities.shape) ) / 2.
        else:
            n_bins = self.n_bins
        
        if intensities.dtype == np.bool:
            bin_values, bin_edges = np.histogram( radii * intensities, bins=n_bins )
        else:
            bin_values, bin_edges = np.histogram( radii, weights=intensities, bins=n_bins )
            
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

        start = 2 # the first two positions are reserved for beam_loc
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
        self.beam_loc = param_vals[:2]
        
        pp = self.cspad.pixel_positions / self.pixel_size
        bin_centers, bin_values = self._bin_intensities_by_radius(self.beam_loc, 
                                                                  pp, raw_image)
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
            blob_circ = plt_patches.Circle(self.beam_loc, 15, fill=False, lw=2, 
                                           ec='orange')
            self._axL.add_patch(blob_circ)
            plt.draw()
                  
        # --------- HERE IS THE OBJECTIVE FUNCTION -- MODIFY TO PLAY -----------
        # TJL: I will think about smart ways to inject other functions in here,
        # but it may also be good to keep this static, once we have something
        # nice and robust.
                
        obj = self._simple_width(bin_values, bar=self.horizontal_cut) - \
              self.max_weight * bin_values.max() + \
              self.peak_weight * n_maxima
        
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
        print "Objective function:  Sum[peak_widths] + %.2f * num_peaks" % self.peak_weight
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

        if self.use_edge_filter:
            image = utils.find_rings(raw_image, threshold=self.threshold, 
                                     minf_size=self.minf_size,
                                     medf_size=self.medf_size)

        # run simplex minimization
        opt_params = optimize.fmin_powell(self._objective, initial_guesses, 
                                          args=(image,), xtol=1e-2, ftol=1e-2)
                                   
        # turn off interactive plotting
        if self.plot_each_iteration: plt.ioff()
                                   
        # un-ravel & inject the param values in the CSPad object
        param_dict = self._unravel_params(opt_params)
        self.cspad.set_many_params(param_dict.keys(), param_dict.values())

        return
    