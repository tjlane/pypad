
"""
This file contains methods to score a geometry optimization method against
known standards: AgBe, Ag.
"""

import numpy as np
from scipy import optimize


#  -- FUNDAMENTAL CONSTANTS --
h = 4.135677516e-15            # Planks constant | eV s
c = 299792458 * 1.0e10         # speed of light  | Ang / s
# ----------------------------


class PowderReference(object):
    """
    A class that provides for scoring a detector geometry based on a known
    powder reference sample.
    """
    
    def __init__(self, lattice_spacing, real_peak_locations, energy, path_length,
                 millers_included=False):
        """
        Generate a powder reference for assessing an x-ray scattering detector
        geometry.
        
        Parameters
        ----------
        lattice_spacing : float
            The latice spacing of the material, in Angstroms
            
        real_peak_locations : ndarray, float
            The locations of peak maxima in real space.
            
        energy : float
            Beam energy in eV
        
        path_length : float
            The distance from the sample to the detector.
            
        millers_included : list of 3-tuples
            A list of (a,b,c) Miller indices to include in the fit. If not
            included, the algorithm will automagically try and match 
        """
        
        self.lattice_spacing = lattice_spacing
        self.real_peak_locations = real_peak_locations
        self.path_length = path_length
        self.energy = energy
        self.k = ( 2.0 * np.pi / (h * c) ) * energy # inverse angstroms
        self._compute_reciprocal_peaks()
        
        if millers_included == False:
            self.millers_included = self._determine_miller_indicies_to_use()
        else:
            self.millers_included = millers_included
        
        return
        
        
    @classmethod
    def agbe(cls, real_peak_locations, energy, path_length):
        """
        Factory method that sets the lattice spacing to that for SilverBehenate.
        """
        
        real_peak_locations = np.sort(real_peak_locations)
        lattice_spacing = 58.3803 # ang
        
        n_peaks = len(real_peak_locations)
        millers_included = [ (0,0,x) for x in range(1, n_peaks+1) ]
        
        return cls(lattice_spacing, real_peak_locations, energy, path_length,
                   millers_included=millers_included)
                   

    @classmethod               
    def silver(cls, real_peak_locations, energy, path_length):
        """
        Factory method that sets the lattice spacing to that for silver.
        """
        
        real_peak_locations = np.sort(real_peak_locations)
        lattice_spacing = 4.09 # ang
        
        n_peaks = len(real_peak_locations)
        millers_included = [(1,1,1), (2,0,0), (2,2,0), (3,1,1)]
        
        return cls(lattice_spacing, real_peak_locations, energy, path_length,
                   millers_included=millers_included)
    
                   
    @classmethod               
    def gold(cls, real_peak_locations, energy, path_length):
        """
        Factory method that sets the lattice spacing to that for gold.
        """

        real_peak_locations = np.sort(real_peak_locations)
        lattice_spacing = 4.076

        n_peaks = len(real_peak_locations)
        millers_included = [(1,1,1), (2,0,0), (2,2,0), (3,1,1)]

        return cls(lattice_spacing, real_peak_locations, energy, path_length,
                   millers_included=millers_included)
    

    def _compute_reciprocal_peaks(self):
        """
        Convert the real-space peaks into reciprocal space.
        """
        self.reciprocal_peak_locations = self.reciprocal(self.real_peak_locations)
    
        
    def reciprocal(self, real_space):
        """
        Convert the real-space peaks into reciprocal space.
        """
        two_thetas = np.arctan(real_space / self.path_length)
        reciprocal_space = 2.0 * self.k * np.sin( 0.5 * two_thetas )
        return reciprocal_space
        
        
    def real_space(self, reciprocal_space):
        """
        Convert from momentum to real-space
        """
        q = reciprocal_space
        real = self.path_length * np.tan(2.0 * np.arcsin( q / (2.0*self.k) ))
        return real
    
        
    def _compute_expected_ring_locations(self):
        """
        Returns the q-space location (in inv. Angstoms) of each ring
        corresponding to the Miller indicies asked for.
        """
        
        expected = np.zeros(len(self.millers_included))
        
        for i,miller_index in enumerate(self.millers_included):
            zf = np.sqrt( np.sum( np.power( miller_index, 2 ) ) )
            expected[i] = (2.0 * np.pi * zf) / self.lattice_spacing
        
        return expected
    
        
    def _match_peaks(self, obsd, expt):
        """
        Automatically match observed powder rings to miller indices. Currently
        will match observed peaks to the closest expected peak.
        """
        
        m_obsd = obsd.copy()
        m_expt = np.zeros_like(m_obsd)
        m_millers = []
        
        assert len(expt) == len(self.millers_included)
        
        for i in range(len(obsd)):
            match_index  = np.argmin( np.abs(expt - obsd[i]) )
            m_expt[i]    = expt[match_index]
            m_millers.append( self.millers_included[match_index] )
        
        return m_obsd, m_expt, m_millers
    
        
    def score(self, verbose=True):
        """
        Compute a score based on the differences between the observed and
        expected locations of the powder rings. A lower score is better,
        with zero being perfect.
        
        Returns
        -------
        total_score : float
            The final total score.
        """
        
        total_score = 0.0
                
        self.obsd = self.reciprocal_peak_locations
        self.expt = self._compute_expected_ring_locations()
        m_obsd, m_expt, m_millers = self._match_peaks(self.obsd, self.expt)
        
        total_score = np.sqrt( np.sum( np.power(m_obsd - m_expt, 2) ) ) / float(len(m_obsd))
        
        if verbose:
                        
            print ""
            print "\t\t\tREFERENCE SAMPLE SCORING"
            print "\t\t\t========================"
            print ""
            print "Peak\t  Miller \tObsd (1/A)\tExpd (1/A)\t  Diff.  "
            print "----\t---------\t----------\t----------\t---------"
            
            for i,miller_index in enumerate(m_millers):
                diff = m_obsd[i] - m_expt[i]
                print "   %d\t%s\t%.4e\t%.4e\t%.2e" % (i, str(miller_index),
                                                     m_obsd[i], m_expt[i], diff)
            print ""
            print "FINAL SCORE:\t%.4e" % total_score
            print ""
        
        return total_score
    
    
    def _update_pl_and_score(self, path_length):
        """
        Helper function for `determine distance`
        """
        
        self.path_length = path_length
        self._compute_reciprocal_peaks()
        
        total_score = self.score(verbose=False)
        
        return total_score
        
    
    def determine_distance(self, verbose=True):
        """
        Given the geometry provided, guess as to the detector distance. 
        
        Does this by minimizing the `score` (squared residuals between the 
        predicted and observed peak positions) as a function of distance.
        """
        
        if verbose:
            print "\t\t\tOPTIMIZING PATH LENGTH"
            print "\t\t\t======================"
        
        p0 = self.path_length
        opt_path_length = optimize.fmin(self._update_pl_and_score, p0)
        
        if verbose:
            print "\nOptimal path length:   %f" % opt_path_length
            print   "Original path length:  %f" % p0
            print "Correction: %f" % (p0 - opt_path_length)
            print ""
            self.score()
        
        return opt_path_length
    
    
    def _update_energy_and_score(self, energy):
        """
        Helper function for `determine distance`
        """
        
        self.energy = energy
        self.k = ( 2.0 * np.pi / (h * c) ) * energy
        self._compute_reciprocal_peaks()
        
        total_score = self.score(verbose=False)
        
        return total_score
    
        
    def determine_energy(self, verbose=True):
        """
        Given the geometry provided, guess as to the beam energy.
        
        Does this by minimizing the `score` (squared residuals between the 
        predicted and observed peak positions) as a function of distance.
        """
        
        if verbose:
            print "\t\t\tOPTIMIZING BEAM ENERGY"
            print "\t\t\t======================"
        
        e0 = self.energy
        opt_energy = optimize.fmin(self._update_energy_and_score, e0)
        
        if verbose:
            print "\nOptimal energy (eV):   %f" % opt_energy
            print   "Original energy (eV):  %f" % e0
            print "Correction: %f" % (e0 - opt_energy)
            print ""
            self.score()
        
        return opt_energy


def test_agbe_score():

    # peak locations in mm, from Jonas
    real_peak_locations = [  2.80296,  5.6334,  8.60124, 11.62404, 14.64684,
                            17.66964, 20.7474, 23.79768, 26.95788 ]

    path_length = 129.0148
    energy      = 9394.363725

    sref = PowderReference.agbe(real_peak_locations, energy, path_length)

    sref.score()

    return