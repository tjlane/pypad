# THIS FILE IS PART OF PyPad, AND IS GOVERENED BY A PERMISSIBILITY LICENSE 
# GOVERNING ITS USE AND DISTRIBUTION. YOU SHOULD HAVE RECIEVED A COPY OF THIS
# LICENSE WITH THE SOFTWARE; IF NOT PROVIDED, WRITE TO <tjlane@stanford.edu>.
#
# AUTHORS:
# TJ Lane <tjlane@stanford.edu>
# Jonas Sellberg <sellberg@slac.stanford.edu>
#
# Apr 30, 2013

"""
score.py

Methods for evaluating optimized geometries.
"""

import yaml
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as plt_patches
from scipy import optimize

from pypad import cspad
from pypad import read
from pypad import utils
from pypad import plot


#  -- FUNDAMENTAL CONSTANTS --
h = 4.135677516e-15            # Planks constant | eV s
c = 299792458 * 1.0e10         # speed of light  | Ang / s
# ----------------------------


class PowderReference(object):
    """
    A class that provides for scoring a detector geometry based on a known
    powder reference sample.
    """
    
    def __init__(self, lattice_spacing, millers, calibration_samples, geometry,
                 energy_guess=9600.0, distance_offset_guess=0.0, 
                 opt_energy=False, unit_cell='cubic'):
        """
        Generate a powder reference for assessing an x-ray scattering detector
        geometry.
        
        Parameters
        ----------
        lattice_spacing : float
            The latice spacing of the material, in Angstroms
            
        millers : list of 3-tuples
            A list of (a,b,c) Miller indices to include in the fit.
            
        calibration_samples : list of dicts
            A list of dictionaries of the form
            
                { filename : x,
                  distance : y }
                  
            which specify the calibration samples to score. The energy these
            data were taking at will be determined automatically.
            
        geometry : cspad.CSPad
            The geometry to score.
        """
        
        self.unit_cell = unit_cell
        self.lattice_spacing = lattice_spacing
        self.millers = millers
        
        self.energy          = energy_guess
        self.distance_offset = distance_offset_guess
        self.opt_energy      = opt_energy
        
        if not isinstance(geometry, cspad.CSPad):
            raise TypeError('`geometry` must be an instance of cspad.CSPad')
        self.cspad = geometry
        
        
        self.sample_peak_locs    = [] # real space
        self.sample_peak_heights = [] # arb units
        self.sample_distances    = [] # relative dist
        
        
        # the following populates the lists above
        self.calibration_samples = calibration_samples
        for sample in self.calibration_samples:
            self._load_calibration_sample(sample['filename'], sample['distance'])
        assert len(self.sample_peak_locs) == len(self.sample_distances)
        
        
        return
        
        
    def _load_calibration_sample(self, filename, distance):
        """
        Loads a calibration file, measures the peak positions of rings in that
        file, and stores that data in the internal structure self.samples
        """
        
        img = read.load_raw_image(filename)
        bin_centers, a = self.cspad.intensity_profile(img, n_bins=800)
        
        #sa = a
        sa = utils.smooth(a, beta=10.0, window_size=20)
        # finds all local maxima, i.e. all points that has both adjacent points with lower values after smoothing
        max_inds = np.where(np.r_[True, sa[1:] > sa[:-1]] & np.r_[sa[:-1] > sa[1:], True] == True)[0]
        real_peak_locations = bin_centers[max_inds]
        
        self.sample_peak_locs.append( np.array(real_peak_locations) )
        self.sample_peak_heights.append( a[max_inds] )
        self.sample_distances.append(float(distance))
        
        print "    found: %d peaks" % len(real_peak_locations)
        
        return
    
        
    @classmethod
    def load(cls, filename):
        """
        Load a score_params yaml file from disk and create a PowderReference
        instance.
        
        Parameters
        ----------
        filename : str
            The name of the score yaml file to load.
        
        Returns
        -------
        cls : score.PowderReference
            An instance containing the data referenced in the input file. See
            the documentation for more information.
        """
        
        f = open(filename, 'r')
        params = yaml.load(f)
        f.close()
        
        print "Loaded: %s" % filename
        
        geom = cspad.CSPad.load( params['geometry'] )
        
        for x in ['lattice', 'millers', 'samples']:
            if not x in params.keys():
                raise IOError('Could not find required input field: '
                              '`%s` in %s' % (x, filename))
  
        kwargs = {}
        if 'energy' in params.keys():
            kwargs['energy_guess'] = float(params['energy'])
        if 'opt_energy' in params.keys():
            kwargs['opt_energy'] = bool(params['opt_energy'])
        if 'initial_offset' in params.keys():
            kwargs['distance_offset_guess'] = float(params['initial_offset'])
        if 'unit_cell' in params.keys():
            kwargs['unit_cell'] = params['unit_cell']
        
        return cls(params['lattice'], params['millers'], params['samples'], 
                   geom, **kwargs)
    
        
    @property
    def k(self):
        return ( 2.0 * np.pi / (h * c) ) * self.energy
        
        
    @property
    def num_samples(self):
        assert len(self.sample_peak_locs) == len(self.sample_distances)
        return len(self.sample_distances)
        
        
    @property
    def num_millers(self):
        return len(self.millers)
    
        
    @property
    def observed(self):
        """
        The observed peak locations, trimmed to match those close to expected
        peak locations (self.expected). Reciprocal space (inv. Ang.).
        
        Returns
        --------
        observed : np.ndarray
            A two-d array, with the first dimension the sample number, the 
            second dimension the peak position in q-space, such that it lines
            up w/self.expected.
        """
        
        observed = np.zeros(( self.num_samples, self.num_millers ))
        expected = self.expected
        
        for i in range(self.num_samples):
            observed[i] = self._match_peaks(i)
                
        return observed
                 
                 
    @property
    def expected(self):
        """
        The expected peak locations. Reciprocal space (inv. Ang.).
        """
        return self._compute_expected_ring_locations()
    
        
    def reciprocal(self, real_space, path_length):
        """
        Convert the real-space peaks into reciprocal space.
        """
        two_thetas = np.arctan(real_space / (path_length))
        reciprocal_space = 2.0 * self.k * np.sin( 0.5 * two_thetas )
        return reciprocal_space
        
        
    def real_space(self, reciprocal_space, path_length):
        """
        Convert from momentum to real-space
        """
        q = reciprocal_space
        real = (path_length) * np.tan(2.0 * np.arcsin( q / (2.0*self.k) ))
        return real
    
        
    def _compute_expected_ring_locations(self):
        """
        Returns the q-space location (in inv. Angstoms) of each ring
        corresponding to the Miller indices asked for.
        """
        
        expected = np.zeros(len(self.millers))
        
        # calculate reciprocal unit vectors b1, b2, b3 for various lattices
        if (self.unit_cell == 'cubic' or self.unit_cell == 'sc' or self.unit_cell == 'bcc' or self.unit_cell == 'fcc'):
            # lattice spacing
            try:
                if (len(self.lattice_spacing) > 1):
                    print "WARNING: edge lengths a = b = c in cubic unit cells, ignoring", self.lattice_spacing[1:]
                a = self.lattice_spacing[0]
            except TypeError:
                a = self.lattice_spacing
            
            # reciprocal unit vectors
            b1 = np.array([2.0*np.pi/a, 0, 0])
            b2 = np.array([0, 2.0*np.pi/a, 0])
            b3 = np.array([0, 0, 2.0*np.pi/a])
            
        elif (self.unit_cell == 'hexagonal'):
            # lattice spacing
            try:
                if (len(self.lattice_spacing) > 2):
                    print "WARNING: edge lengths a = b in cubic unit cells, ignoring", self.lattice_spacing[2:]
                a = self.lattice_spacing[0]
                c = self.lattice_spacing[1]
            except TypeError:
                print "ERROR: need to specify two edge lengths in hexagonal unit cells, aborting."
                sys.exit(1)
            
            # reciprocal unit vectors
            b1 = np.array([2.0*np.pi/(a*np.sqrt(3)), 2.0*np.pi/a, 0])
            b2 = np.array([-2.0*np.pi/(a*np.sqrt(3)), 2.0*np.pi/a, 0])
            b3 = np.array([0, 0, 2.0*np.pi/c])
            
        else:
            print "ERROR: unknown unit cell `%s`, aborting." % self.unit_cell
            sys.exit(1)
        
        # calculate reciprocal lattice vector for each set of Miller indices,
        # its length determines the momentum transfer
        for i,miller_index in enumerate(self.millers):
            # reciprocal lattice vector
            G = miller_index[0]*b1 + miller_index[1]*b2 + miller_index[2]*b3
            
            # momentum transfer q = |G|
            expected[i] = np.sqrt( np.sum( np.power( G, 2 ) ) )
        
        return expected
    
    
    def _match_peaks(self, sample_index):
        """
        Automatically match observed powder rings to miller indices. Currently
        will match observed peaks to the closest expected peak weighted by its
        peak intensity.
        
        Returns
        --------
        m_observed: np.ndarray
            A one-d array containing the q-positions of the observed peaks that
            matches the expected peaks the best
        """
        
        observed = self.reciprocal(self.sample_peak_locs[sample_index], \
                               self.sample_distances[sample_index] + \
                                   self.distance_offset )
        expected = self.expected
        
        ### New algorithm by Jonas ###
        # Will optimize the observed peak indices that corresponds to the expected peaks
        # based on the proximity of the peak position in q weighted by its peak intensity.
        # This algorithm assumes:
        # (1) all expected peaks exist in observed peaks
        #     (i.e. has to be found by the peak finder in _load_calibration_sample()
        #      and be present on the detector at all detector distances)
        # (2) the peaks in observed and expected are ordered after increasing q
        
        observed_matches = np.zeros( self.num_millers, dtype=np.int32 )
        assert len(expected) == self.num_millers
        
        opt = 1e30
        for i in utils.multi_for(map( xrange, np.arange(self.num_millers), np.ones(self.num_millers, dtype=np.int32 )*len(observed) - np.arange(self.num_millers)[::-1] )):
            # only try peak index combinations where all indices are different
            if len(i) == len(set(i)):
                indices = [x for x in i] # can't access indices with tuple, need to make it into a list
                new_opt = sum(((expected-observed[indices])/self.sample_peak_heights[sample_index][indices])**2)
                if new_opt < opt:
                    opt = new_opt
                    observed_matches[:] = indices
        
        m_observed = observed[observed_matches]
        
        ### OPTION #2 by Jonas ###
        # currently doesn't work since the optimization algorithm fails for discrete variables,
        # could be worth trying if a good linear optimizer for discrete variables is found. This could be a candidate:
        # http://stackoverflow.com/questions/5179889/optimization-problem-in-python
        
#        def peak_objective(args):
#            """
#            The peak objective function -- the difference between the measured/
#            expected peak locations weighted by their intensity.
#            """
#            
#            # make sure that indices are integer numbers
#            for i in range(len(args)):
#                if args[i] >= (len(observed) - 1):
#                    args[i] = len(observed) - 1
#                elif args[i] < 0:
#                    args[i] = 0
#                else:
#                    args[i] = round(i)
#            args = np.array( args, dtype=np.int32 ).flatten()
#            
#            obj = (expected - observed[args])/self.sample_peak_heights[sample_index][args]
#            
#            return obj.flatten()
#        
#        i0 = np.arange(self.num_millers)
#        #opt = optimize.leastsq(peak_objective, i0, args=(expected, observed), full_output=1)
#        opt = optimize.leastsq(peak_objective, i0, full_output=1)
#        
#        observed_matches = np.array( opt[0], dtype=np.int32 )
#        m_observed = observed[observed_matches]
        
        ### OPTION #3 by Jonas ###
        # A similar algorithm to the one currently used, except that it reduces the
        # number of loops by adjusting the lower range of the inner for loop depending
        # on the outer loop index. Unfortunately it is hard to generalize for arbitrary
        # number of Bragg peaks (i.e. self.num_millers)
        
#        observed_matches = np.zeros( self.num_millers, dtype=np.int32 )
#        assert len(expected) == self.num_millers
#        
#        opt = 1e30
#        for i in range(len(observed) - 1):
#            for j in range(i + 1, len(observed)):
#                indices = [i, j]
#                new_opt = sum(((expected-observed[indices])/self.sample_peak_heights[sample_index][indices])**2)
#                if new_opt < opt:
#                    opt = new_opt
#                    observed_matches[:] = indices
        
        return m_observed
    
    
    def evaluate(self):
        """
        The central function -- evaluates the geometry by fitting the energy,
        sample-to-detector distance by optimizing the expected and observed
        positions of powder rings on the detector.
        """
        
        # TJL note to self : we probably want to *also* evaluate each quad
        # separately and look to see if one is not-so-good.
        
        def objective(args):
            """
            The objective function -- the difference between the measured/
            expected peak locations.
            """
            if self.opt_energy:
                self.distance_offset = args[0]
                self.energy          = args[1]
            else:
                self.distance_offset = args[0]
            obj = np.abs(self.observed - self.expected)
            return obj.flatten()
        
        
        x0 = (self.distance_offset, self.energy)
        opt = optimize.leastsq(objective, x0, full_output=1)
        
        print
        print " --- Optimized Energy & Detector Distance --- "
        print " Detector offset: %.3f mm " % opt[0][0]
        print " Energy:          %.4f keV" % (opt[0][1] / 1000.0,)
        print
        print " Total Residuals: %f inv. Angstroms" % float( np.sum(np.abs( opt[2]['fvec'] )) )
        print " Residuals for each peak:"
        print opt[2]['fvec']
        print 
        
        return
        
        
    def display(self):
        """
        Show an image of the current fit.
        """
        
        # make a plot!
        self._fig = plt.figure(figsize=(12,6))
        self._axL = plt.subplot(121)
        self._axR = plt.subplot(122)
        
        # provide functionality so that the user can cycle through images by
        # pressing n/l
        
        self._current_image = 0
        connection = plt.connect('key_press_event', self._on_keypress)
        self._plot_cal_sample(self._current_image)
        
        return
    
        
    def _plot_cal_sample(self, index):
        """
        Plots a single calibration sample.
        """
        
        self._axL.cla()
        self._axR.cla()
        
        if (index < 0) or (index >= self.num_samples):
            print "Cannot access sample: %d" % index
            print "Total %d samples available" %  self.num_samples
            return
        
        print
        print "Plotting calibration sample: %d" % index
        print "  (may take a moment)"
        
        # load up the calibration sample requested
        fn = self.calibration_samples[index]['filename']
        d  = self.calibration_samples[index]['distance'] + self.distance_offset
        print "  distance: %.2f mm" % d
        
        # -- plot left panel, the assemled image with ring predictions overlaid
        img = read.load_raw_image(fn)
        plot.imshow_cspad( self.cspad(img), vmin=0, ax=self._axL )
        
        # plot circles on the image, over where the powder rings should be
        # note that (1000,1000) is where the center is in pixel units
        # for our fxn imshow_cspads
        
        real_expected = self.real_space(self.expected, d) / 0.10992
        for r in real_expected:
            blob_circ = plt_patches.Circle((1000,1000), r, fill=False, lw=1, ec='white')
            self._axL.add_patch(blob_circ)
        
        # plot beam center
        beam_center = plt_patches.Circle((1000,1000), 2, fill=True, lw=1, color='r')
        self._axL.add_patch(beam_center)
        
        
        # --- plot the right image
        n_bins = 800
        
        bin_centers, a = self.cspad.intensity_profile(img, n_bins=n_bins)
        q_bin_centers = self.reciprocal(bin_centers, d)
        self._axR.plot(q_bin_centers, a / a.max(), color='orange', lw=4)
        
        for i in range(4):
            bin_centers, a = self.cspad.intensity_profile(img, n_bins=n_bins, quad=i)
            a /= a.max()
            a += 1.0 * i + 1.0
            q_bin_centers = self.reciprocal(bin_centers, d)
            self._axR.plot(q_bin_centers, a, color=plot.quad_colors[i], lw=2)

        self._axR.text(0.5, -0.2, 'All Quads')
        self._axR.text(0.5,  0.8, 'Quad 0')
        self._axR.text(0.5,  1.8, 'Quad 1')
        self._axR.text(0.5,  2.8, 'Quad 2')
        self._axR.text(0.5,  3.8, 'Quad 3')
        
        self._axR.set_ylim([-0.3, 5.3])

        self._axR.vlines(self.expected, 0, a.max()*1.2, color='k', linestyles='dashed')
        # self._axR.vlines(self.observed[index,:], 0, a.max(), color='r', linestyles='dashed')

        self._axR.set_xlabel(r'q ($\AA^{-1}$)')
        self._axR.set_ylabel('Intensity')
        
        plt.show()
        print
        
        return
        
    
    def _on_keypress(self, event):
        if event.key == 'l':
            if (self._current_image % self.num_samples) == 0:
                self._current_image += self.num_samples
            self._current_image -= 1
            self._plot_cal_sample(self._current_image)
        elif event.key == 'n':
            self._current_image += 1
            if (self._current_image % self.num_samples) == 0:
                self._current_image -= self.num_samples
            self._plot_cal_sample(self._current_image)
        else:
            return
    

def test_agbe_score():

    # peak locations in mm, from Jonas
    real_peak_locations = [  2.80296,  5.6334,  8.60124, 11.62404, 14.64684,
                            17.66964, 20.7474, 23.79768, 26.95788 ]

    path_length = 129.0148
    energy      = 9394.363725

    sref = PowderReference.agbe(real_peak_locations, energy, path_length)

    sref.score()

    return
