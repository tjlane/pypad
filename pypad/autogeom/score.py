
"""
This file contains methods to score a geometry optimization method against
known standards: AgBe, Ag.
"""

import yaml

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
                 opt_energy=False):
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
    def obsd(self):
        """
        The observed peak locations, trimmed to match those close to expected
        peak locations (self.expt). Reciprocal space (inv. Ang.).
        
        Returns
        --------
        obsd : np.ndarray
            A two-d array, with the first dimension the sample number, the 
            second dimension the peak position in q-space, such that it lines
            up w/self.expt.
        """
        
        obsd = np.zeros(( self.num_samples, self.num_millers ))
        expt = self.expt
        
        for i in range(self.num_samples):
            obsd[i] = self._match_peaks(i)
                
        return obsd
                 
                 
    @property
    def expt(self):
        """
        The expected peak locations. Reciprocal space (inv. Ang.).
        """
        return self._compute_expected_ring_locations()
    
        
    def reciprocal(self, real_space, path_length):
        """
        Convert the real-space peaks into reciprocal space.
        """
        two_thetas = np.arctan(real_space / (path_length + self.distance_offset))
        reciprocal_space = 2.0 * self.k * np.sin( 0.5 * two_thetas )
        return reciprocal_space
        
        
    def real_space(self, reciprocal_space, path_length):
        """
        Convert from momentum to real-space
        """
        q = reciprocal_space
        real = (path_length + self.distance_offset) * np.tan(2.0 * np.arcsin( q / (2.0*self.k) ))
        return real
    
        
    def _compute_expected_ring_locations(self):
        """
        Returns the q-space location (in inv. Angstoms) of each ring
        corresponding to the Miller indicies asked for.
        """
        
        expected = np.zeros(len(self.millers))
        
        for i,miller_index in enumerate(self.millers):
            zf = np.sqrt( np.sum( np.power( np.array(miller_index), 2 ) ) )
            expected[i] = (2.0 * np.pi * zf) / self.lattice_spacing
        
        return expected
    
        
    def _match_peaks(self, sample_index):
        """
        Automatically match observed powder rings to miller indices. Currently
        will match observed peaks to the closest expected peak.
        """
        
        # (1) match observed peaks to their closest expt peak
        # (2) for expt peaks with more than one match, keep only the tallest
        #     peak
        
        obsd = self.reciprocal(self.sample_peak_locs[sample_index], \
                               self.sample_distances[sample_index])
        expt = self.expt
                               
        obsd_matches = np.zeros( len(obsd), dtype=np.int32 )
        m_obsd = np.zeros_like(expt)
        
        # match each obsd value with its closest expt value
        for i in range(len(obsd)):
            match_index     = np.argmin( np.abs(expt - obsd[i]) )
            obsd_matches[i] = match_index
            
        # for each expt value, if has multi matches, choose best
        for i in range(self.num_millers):
            w = (obsd_matches == i)
            if np.sum(w) > 1:
                m_obsd[i] = obsd[ np.argmax( self.sample_peak_heights[sample_index][w] ) ]
        
        return m_obsd
    
        
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
            # if self.opt_energy:
            #     self.distance_offset = args[0]
            #     self.energy          = args[1]
            # else:
            #     self.distance_offset = args[0]
            obj = np.abs(self.obsd - self.expt)
            return obj.flatten()
        
        
        x0 = (self.distance_offset, self.energy)
        opt = optimize.leastsq(objective, x0, full_output=1)
        
        print
        print " --- Optimized Energy & Detector Distance --- "
        print " Detector offest: %.2f mm " % opt[0][0]
        print " Energy:          %.3f keV" % (opt[0][1] / 1000.0,)
        print
        print " Total Residuals: %f inv. Angstroms" % float( np.sum(opt[2]['fvec']) )
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
        
        if (index < 0) or (index > self.num_samples):
            print "Cannot access sample: %d" % index
            print "Total %d samples available" %  self.num_samples
            return
        
        print
        print "Plotting calibration sample: %d..." % index
        print "  (may take a moment)"
        
        # load up the calibration sample requested
        fn = self.calibration_samples[index]['filename']
        d  = self.calibration_samples[index]['distance'] + self.distance_offset
        
        # print "Plotting sample: %s" % fn
        print "  distance: %.2f mm" % d
        
        # -- plot left panel, the assemled image with ring predictions overlaid
        img = read.load_raw_image(fn)
        plot.imshow_cspad( self.cspad(img), vmin=0, ax=self._axL)
        
        # plot circles on the image, over where the powder rings should be
        # note that (1000,900) is where the center is in pixel units
        # for our fxn imshow_cspads
        
        real_expt = self.real_space(self.expt, d) / 0.10992
        for r in real_expt:
            blob_circ = plt_patches.Circle((1000,900), r, fill=False, lw=1, ec='white')
            self._axL.add_patch(blob_circ)
        
        
        # --- plot the right image
        n_bins = 800
        
        for i in range(4):
            bin_centers, a = self.cspad.intensity_profile(img, n_bins=n_bins, quad=i)
            a /= a.max()
            a += 0.7 * i
            q_bin_centers = self.reciprocal(bin_centers, d)
            self._axR.plot(q_bin_centers, a, color=plot.quad_colors[i], lw=2)

        self._axR.vlines(self.expt, 0, a.max(), color='k', linestyles='dashed')
        self._axR.vlines(self.obsd[index,:], 0, a.max(), color='r', linestyles='dashed')

        self._axR.set_xlabel(r'q ($\AA^{-1}$)')
        self._axR.set_ylabel('Intensity')
        
        plt.show()
        print
        
        return
        
    
    def _on_keypress(self, event):
        if event.key == 'l':
            self._current_image -= 1
            self._plot_cal_sample(self._current_image)
        elif event.key == 'n':
            self._current_image += 1
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