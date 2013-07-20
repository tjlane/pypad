
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
cspad.py

An interface to the CSPad geometry. Converts LCLS optical metrologies into
explicit pixel positions.
"""

import sys
import os
import cPickle
import h5py
import re
from glob import glob
from os.path import join as pjoin

import numpy as np
import scipy.ndimage.interpolation as interp
import matplotlib.pyplot as plt
from matplotlib.nxutils import points_inside_poly

from pypad import utils
from pypad import read
from pypad import export
from pypad import default


class BasisGrid(object):
    """
    A class representing a set of rectangular grids in space -- specifically,
    x-ray scattering detectors. Does not contain all the metadata associated
    with a full-fledged Detector class (e.g. the wavelength, etc).

    Note that the geometry below is definied in "slow" and "fast" scan
    dimensions. These are simply the two dimensions that define the plane
    a single rectangular pixel grid lives in. They may also be called the y and
    x dimensions without any loss of generality.

    The convention here -- and in all of ODIN -- is one of Row-Major ordering,
    which is consistent with C/python. This means that y is the slow dim, x is
    the fast dim, and when ordering these two they will appear as (slow, fast).

    Note on units: units are arbitrary -- all the units must be the same for
    this to work. We don't keep track of units here.

    The class is a set of rectangular grids, with each grid defined by four
    quantities:

        -- p vector : DEFINES A GRIDS POSITION IN SPACE.
                      The vector between a chosen origin (possibly interaction
                      site) and the corner of the grid that is smallest in both
                      slow and fast (x/y) dimensions of the coordinate system.
                      Usually this will be the "bottom left" corner, but due to
                      the generality of the coordinates used, this is not
                      necessarily true.

        -- s/f vect : DEFINES A GRIDS ORIENTATION IN SPACE
                      Vectors pointing along the slow/fast-scan direction,
                      respectively. These define the plane containing the pixels.
                      The magnitudes of these vectors defines the size of the
                      pixel in that dimension.

        -- shape    : DEFINES GRID DIMENSIONS
                      The number of pixels in the fast/slow direction. Ints.
    """


    def __init__(self, list_of_grids=[]):
        """
        Initialize a BasisGrid object.

        Parameters
        ----------
        list_of_grids : list
            A list of tuples of the form  (p, s, f, shape). See the doc
            for the `add_grid` method on this class for more information. May
            be an empty list (default) in which case a GridList with no pixels
            is created.

        See Also
        --------
        add_grid
        add_grid_using_center
        """

        if not type(list_of_grids) == list:
            raise TypeError('`list_of_grids` must be a list')

        self._num_grids = 0
        self._ps        = [] # p-vectors
        self._ss        = [] # slow-scan vectors
        self._fs        = [] # fast-scan vectors
        self._shapes    = [] # shapes

        if len(list_of_grids) > 0:
            for grid in list_of_grids:
                self.add_grid(*grid)

        return
        

    def _check_valid_basis(self, p, s, f, shape):
        """
        Check to make sure that all the inputs look good.
        """

        if not (p.shape == (3,)) and (s.shape == (3,)) and (f.shape == (3,)):
            raise ValueError('`p`, `s`, `f` must be 3-vectors')

        if not (len(shape) == 2):
            raise ValueError('`shape` must be len 2')

        return


    def _assert_list_sizes(self):
        """
        A simple sanity check
        """
        assert len(self._ps)     == self.num_grids
        assert len(self._ss)     == self.num_grids
        assert len(self._fs)     == self.num_grids
        assert len(self._shapes) == self.num_grids
        return


    @property
    def num_pixels(self):
        """
        Return the total number of pixels in the BasisGrid.
        """
        n = np.sum([np.product(self._shapes[i]) for i in range(self.num_grids)])
        return int(n)


    @property
    def num_grids(self):
        return self._num_grids


    def add_grid(self, p, s, f, shape):
        """
        Add a grid (detector array) to the basis representation.

        Parameters
        ----------
        p : np.ndarray, float
            3-vector from the origin to the pixel on the grid with
            smallest coordinate in all dimensions.

        s : np.ndarray, float
            3-vector pointing in the slow scan direction

        f : np.ndarray, float
            3-vector pointing in the slow scan direction

        shape : tuple or list of float
            The number of pixels in the (slow, fast) directions. Len 2.

        See Also
        --------
        add_grid_using_center
        """
        self._check_valid_basis(p, s, f, shape)
        self._ps.append(p)
        self._ss.append(s)
        self._fs.append(f)
        self._shapes.append(shape)
        self._num_grids += 1
        self._assert_list_sizes()
        return


    def add_grid_using_center(self, p_center, s, f, shape):
        """
        Add a grid (detector array) to the basis representation. Here, though,
        the p-vector points to the center of the array instead of the slow/fast
        smallest corner.

        Parameters
        ----------
        p_center : np.ndarray, float
            3-vector from the origin to the center of the grid.

        s : np.ndarray, float
            3-vector pointing in the slow scan direction

        f : np.ndarray, float
            3-vector pointing in the slow scan direction

        shape : tuple or list of float
            The number of pixels in the (slow, fast) directions. Len 2.
        """

        p_center = np.array(p_center)
        if not p_center.shape == (3,):
            raise ValueError('`p_center` must have shape (3,)')

        # just compute where `p` is then add the grid as usual
        x = (np.array(shape) - 1)
        center_correction =  ((x[0] * s) + (x[1] * f)) / 2.
        p  = p_center.copy()
        p -= center_correction

        self.add_grid(p, s, f, shape)

        return


    def get_grid(self, grid_number):
        """
        Return a grid for grid `grid_number`.

        Parameters
        ----------
        grid_number : int
            The index of the grid to get.

        Returns
        -------
        p_center : np.ndarray, float
            3-vector from the origin to the center of the grid.

        s : np.ndarray, float
            3-vector pointing in the slow scan direction

        f : np.ndarray, float
            3-vector pointing in the slow scan direction

        shape : tuple or list of float
            The number of pixels in the (slow, fast) directions. Len 2.
        """

        if grid_number >= self.num_grids:
            raise ValueError('Only %d grids in object, you asked for the %d-th'
                             ' (zero indexed)' % (self.num_grids, grid_number))

        grid_tuple = (self._ps[grid_number], self._ss[grid_number],
                      self._fs[grid_number], self._shapes[grid_number])

        return grid_tuple


    def get_grid_corners(self, grid_number):
        """
        Return the positions of the four corners of a grid.

        Parameters
        ----------
        grid_number : int
            The index of the grid to get the corners of.

        Returns
        -------
        corners : np.ndarray, float
            A 4 x 3 array, where the first dim represents the four corners, and
            the second is x/y/z. Note one corner is always just the `p` vector.
        """

        if grid_number >= self.num_grids:
            raise ValueError('Only %d grids in object, you asked for the %d-th'
                             ' (zero indexed)' % (self.num_grids, grid_number))

        # compute the lengths of the parallelogram sides
        s_side = self._fs[grid_number] * float(self._shapes[grid_number][0])
        f_side = self._ss[grid_number] * float(self._shapes[grid_number][1])
        pc = self._ps[grid_number].copy()

        corners = np.zeros((4,3))

        corners[0,:] = pc
        corners[1,:] = pc + s_side
        corners[2,:] = pc + f_side
        corners[3,:] = pc + s_side + f_side

        return corners


    def to_explicit(self):
        """
        Return the entire grid as an n x 3 array, defining the x,y,z positions
        of each pixel.

        Returns
        -------
        xyz : np.ndarray, float
            An N x 3 array of the x,y,z positions of each pixel. Note that this
            is a flattened version of what you get for each grid individually
            using `grid_as_explicit`.

        See Also
        --------
        grid_as_explicit
        """
        ex_grids = [ self.grid_as_explicit(i) for i in range(self.num_grids) ]
        xyz = np.concatenate([ g.reshape((g.shape[0]* g.shape[1], 3)) for g in ex_grids ])
        return xyz


    def grid_as_explicit(self, grid_number):
        """
        Get the x,y,z coordiantes for a single grid.

        Parameters
        ----------
        grid_number : int
            The index of the grid to get.

        Returns
        -------
        xyz : np.ndarray, float
            An (shape) x 3 array of the x,y,z positions of each pixel

        See Also
        --------
        to_explicit
        """

        p, s, f, shape = self.get_grid(grid_number)

        # xyz = i * s + j * f, where i,j are ints running over range `shape`
        mg = np.mgrid[0:shape[0]-1:1j*shape[0], 0:shape[1]-1:1j*shape[1]]
        xyz = np.outer(mg[0].flatten(), s) + np.outer(mg[1].flatten(), f)
        xyz += p # translate
        xyz = xyz.reshape( (shape[0], shape[1], 3) )

        return xyz


pixel_size = 0.10992 # sort of a universal constant, in mm
_array_sizes = {'quad_rotation'       : (4,),
                'quad_offset'         : (4,2),
                'quad_offset_bydiag'  : (4,2)}

class CSPad(object):
    """
    This is a container class for saving, loading, and interacting with the 
    CSPad geometry. It contains an interface to read the optical metrologies
    performed at LCLS by the CSPAD group, and can write many file formats.
    
    Further, this class lets you visualize data on a specific geometry.
    Imagine you have a 'raw_image', an np array from psana with the measured
    intensities for each pixel. You also have a directory 'my_params', which
    contains a directory structure like 
    
    my_params/
      rotation/
        0-end.data
      tilt/
        0-end.data
      ...    
    
    Then, you should be able to turn those mysterious parameters into a CSPad
    geometry by doing something like this:
    
    >>> from pypad import cspad, plot
    >>> geom = cspad.CSPad('path/to/metrology.txt')
    >>> assembled_image = geom(raw_image)
    >>> plot.imshow_cspad(assembled_image)
    
    Of course, the 'geom' object can be re-used to assemble many images.
    """
    
    def __init__(self, metrology, quad_offset=np.zeros((4,2)), 
                 quad_rotation=np.zeros(4), verbose=True):
        """
        Initialize an instance of CSPad, corresponding to a single CSPad
        geometry.             
        """
        
        self.verbose = verbose
        self.metrology = metrology
        
        self._param_list = _array_sizes.keys()
        self.pixel_size  = pixel_size
        
        if not quad_offset.shape == _array_sizes['quad_offset']:
            raise ValueError('quad_offset must have shape (4,2), got: %s' % str(quad_offset.shape))
        if not quad_rotation.shape == _array_sizes['quad_rotation']:
            raise ValueError('quad_rotation must have shape (4,), got: %s' % str(quad_rotation.shape))
        
        
        # internalize params
        self.quad_offset   = quad_offset
        self.quad_rotation = quad_rotation
        
        
        # read the optical metrology
        if type(metrology) == str:
            print "Loading metrology from: %s" % metrology
            qms = self._read_metrology(metrology)

                    
        elif type(metrology) == np.ndarray:
            if not metrology.shape == (4,32,3):
                raise ValueError('`metrology` array is wrong shape. Must be (4,32,3),'
                                 ' got %s' % str(metrology.shape))
            qms = metrology
            
        elif type(metrology) == BasisGrid:
            
            # in this case we assume that the overall positions and rotations
            # of the geometry are more or less set, so we set the base rotations
            # to zero and don't reflect over x
            
            self._metrology_basis    = metrology
            self._reflect_xaxis      = False
            self._base_quad_rotation = [0.0, 0.0, 0.0, 0.0]  # deg ccw from upstream
            
            # skip the rest, which is about converting an optical metrology
            return 
                                 
        else:
            raise TypeError('`metrology` must be {str, np.ndarray, BasisGrid}')


        # generate the _metrology_basis -- this is the first of two BasisGrid
        # objects that form the back end of a CSPad object. This BasisGrid has
        # no information encoded about the relative positions of the quads.
        # That information is incorporated in _generate_positional_basis()
        
        self._metrology_basis = BasisGrid()
        self._reflect_xaxis = True
        self._base_quad_rotation = [90.0, 0.0, 270.0, 180.0] # deg ccw from upstream
        
        for q in range(4):
            if self.verbose: print "\nParsing: quad %d" % q
            for two_by_one in range(8):
                asic_geoms = self._twobyone_to_bg(qms[q], two_by_one)
                for asic in range(2):
                    self._metrology_basis.add_grid( *asic_geoms[asic] )
        
        return
    
        
    def __call__(self, raw_image):
        """
        Takes a raw image and assembles it into a two-dimensional
        view.
        
        Returns
        -------
        assembled_image : ndarray, float
            The assembled image.
        """
        if not raw_image.shape == (4,16,185,194):
            raise ValueError('CSPad cannot assemble image: wrong shape. Raw '
                             'image must be shape (4,16,185,194), which is '
                             '(quads, ASICs, slow, fast). Got an image w/shape:'
                             '%s' % str(raw_image.shape))
        return self._assemble_image(raw_image)
    
        
    def dilate(self, dilation):
        """
        Apply an isotropic dilation to the CSPad, in the same way the motors
        move to make the central beam pass bigger. Can also shrink the dilation
        by passing a negative value for `dilation`.
        
        Parameters
        ----------
        dilation : float
            The change in size of the beam pass (along one edge), in mm.
        """
        
        if np.abs(dilation) > 10.0:
            print("Warning: the maximum possible dialtion at CXI as of Apr 2013"
                  " was 10.0mm. You asked for a dilation of %f mm" % dilation)
        
        dilation_offset = np.zeros((4, 2))
        step_size = dilation / 2.0 # this gets applied to both sizes
        
        # move each quad outward along a diagonal
        dilation_offset[0,:] = np.array([ step_size, -step_size])
        dilation_offset[1,:] = np.array([ step_size,  step_size])
        dilation_offset[2,:] = np.array([-step_size,  step_size])
        dilation_offset[3,:] = np.array([-step_size, -step_size])
        
        self.quad_offset[:,:2] += dilation_offset
        return
    
    
    @staticmethod
    def _unit(v):
        return v / np.linalg.norm(v)
        
        
    def _read_metrology(self, metrology_file):
        """
        Parse a flat text file containing a metrology. We assume the metrology
        is of the form:
        
        # quad 0
        1 x1 y1 z1
        2 x2 y2 z2
        3 ...
        
        # quad 1
        1 x1 y1 z1
        2 x2 y2 z2
        3 ...
        
        ...
        
        """
        
        # read the metrology and discard the first col
        met = np.genfromtxt(metrology_file)[:,1:]
        
        if not met.shape == (32*4, 3):
            raise IOError('metrology file appears to be in the wrong format... '
                          'could not understand file format of: %s' % metrology_file)
        
        quad_metrologies = np.array( np.vsplit(met, 4) )
        return quad_metrologies
    
        
    def _qc_angle(self, v, w, two_by_one_index, tol=0.148):
        """
        Perform a quality control check on the angle between two basis vectors.
        
        A note on the `tol` value:
        
        The longest 2x1 side is 388 pixels long. Therefore to achive single
        pixel accuracy in the metrology, the angle between the two 2x1 sides
        should be less than
        
            theta = | arctan( 1 pixel / 388 pixels ) - pi / 2 |
        
        which is ~0.0026 rad = 0.148 and this is what we've
        set the tolerance to for now.
        """
        
        # compute the angle between the vectors
        value = np.degrees( np.arcsin(np.dot(v, w) / ( np.linalg.norm(v) * np.linalg.norm(w) ) ))
        
        if not np.abs(value) <= tol:
            if self.verbose:
                print "WARNING: Metrology quality control failed for 2x1: %d" % two_by_one_index
                print '--> s/f vectors are not orthogonal :: enforcing orthogonality!'
                print "    Angle: %f // tol: %f" % (value, tol)
            passed = False
        else:
            passed = True
            
        return passed
    

    def _twobyone_to_bg(self, quad_metrology, two_by_one_index):
        """
        Convert a 2x1 in the optical metrology into two basis grid elements.

        Parameters
        ----------
        quad_metrology : np.ndarray, float
            A 32 x 3 array of the corner positions of each 2x1 on the quad.

        two_by_one_index : int
            The index of the 2x1 to check out.
            
        Optional Parameters
        -------------------
        tol : float
            The tolerance for a quality control assessment that is performed on
            the optical measurement, including checking things that should be
            orthogonal are, distances are correct, etc. The float `tol` can be
            interpreted as a tolerane for error for each of these in units of
            mm.
            
        Returns
        -------
        bgs : tuple
            A 2-tuple of (p,s,f,shape) for each ASIC in the two by one.
        """

        if not quad_metrology.shape == (32,3):
            raise ValueError('Invalid quad_metrology, must be shape (32,3), got: %s' % str(quad_metrology.shape))
        

        # below is the sequence in which each 2x1 is optically measured, and
        # thus this sequence is read out of the metrology
        
        scan_sequence = np.array([1,0,3,2,4,5,7,6])
        shape = (185, 194) # always the same, for each ASIC

        # rip out the correct corner positions for the 2x1 we want
        i = int(np.where( scan_sequence == two_by_one_index )[0]) * 4
        xyz_2x1 = quad_metrology[i:i+4,:] # the four corners of the 2x1

        
        # The metrology is for sensor points that are actually *outside* the
        # physical ASIC chip. To account for this, we apply a `p_offset`,
        # which effectively centers the ASIC between these measured points.
        # Practically, we take the average center of each side of the rectangle
        # measured in the metrology
        
        # 0.10992 mm : pixel size
        # 0.27480 mm : gap between ASICS in a 2x1
        
        fl = 2 * 194 * 0.10992 + 0.27480  # length of the long side of a 2x1
        sl = 185 * 0.10992                # length of the short side of a 2x1
        h  = np.sqrt( fl*fl + sl*sl )     # length of hypothenuse of an ASIC

        # Next, we build a basis grid representation of each 2x1, which consists
        # of two individual ASICS. Each ASIC has an "origin" which denotes the
        # position of the first pixel in memory -- this is denoted p0 for the 
        # first ASIC and p1 for the second. The values "s" and "f" define the
        # grid of pixels along the slow and fast scan directions (here, short
        # and long sizes) of the 2x1 respectively.
        
        # The s/f vectors are computed by averaging the vectors that define each
        # side in the optical measurement.
        
        # The p vector is measured based on placing a point the correct distance
        # along the diagonal between two measured points. This is checked for
        # correctness by ensuring that the diagonal used is orthogonal to the
        # other 2x1 diagonal measured
        
        # NOTE : later, we swap x to reach CXI convention -- this is done after
        #        we position the quads in _generate_positional_basis()
        
        if two_by_one_index in [0,1]:

            s = 0.10992 * self._unit( (xyz_2x1[0,:] - xyz_2x1[1,:]) +
                                      (xyz_2x1[3,:] - xyz_2x1[2,:]) )
            f = 0.10992 * self._unit( (xyz_2x1[2,:] - xyz_2x1[1,:]) +
                                      (xyz_2x1[3,:] - xyz_2x1[0,:]) )
            
            diagonal = (xyz_2x1[3,:] - xyz_2x1[1,:]) / 1000.
            center = np.sum(xyz_2x1, axis=0) / 1000.
            offset = (np.linalg.norm(diagonal) - h) / 2.
            
            p0 = xyz_2x1[1,:] / 1000. + offset * self._unit(center)
            p1 = p0 + shape[1] * f + 0.27480 * self._unit(f) # for 3px gap


        elif two_by_one_index in [2,3,6,7]:

            s = 0.10992 * self._unit( (xyz_2x1[1,:] - xyz_2x1[2,:]) +
                                      (xyz_2x1[0,:] - xyz_2x1[3,:]) )
            f = 0.10992 * self._unit( (xyz_2x1[3,:] - xyz_2x1[2,:]) + 
                                      (xyz_2x1[0,:] - xyz_2x1[1,:]) )
                                      
            diagonal = (xyz_2x1[0,:] - xyz_2x1[2,:]) / 1000.
            center = np.sum(xyz_2x1, axis=0) / 1000.
            offset = (np.linalg.norm(diagonal) - h) / 2.

            p0 = xyz_2x1[2,:] / 1000. +  offset * self._unit(center)
            p1 = p0 + shape[1] * f + 0.27480 * self._unit(f) # for 3px gap


        elif two_by_one_index in [4,5]:

            s = 0.10992 * self._unit( (xyz_2x1[2,:] - xyz_2x1[3,:]) +
                                      (xyz_2x1[1,:] - xyz_2x1[0,:]) )
            f = 0.10992 * self._unit( (xyz_2x1[0,:] - xyz_2x1[3,:]) +
                                      (xyz_2x1[1,:] - xyz_2x1[2,:]) )

            diagonal = (xyz_2x1[1,:] - xyz_2x1[3,:]) / 1000.
            center = np.sum(xyz_2x1, axis=0) / 1000.
            offset = (np.linalg.norm(diagonal) - h) / 2.

            p0 = xyz_2x1[3,:] / 1000. +  offset * self._unit(center)
            p1 = p0 + shape[1] * f + 0.27480 * self._unit(f) # for 3px gap


        else:
            raise ValueError('two_by_one_index must be in 0...7')
            
            
        # --- perform some quality control checks ---
                
        # (1) ensure s/f orthogonal
        
        self._qc_angle(s, f, two_by_one_index)
        
        # --- end QC ---------------------------------
        
        
        # no matter what, correct the s/f vectors so they are orthogonal
        axis = np.cross(s, f)

        # the angle between s/f
        theta = np.arccos(np.dot(s, f) / ( np.linalg.norm(s) * np.linalg.norm(f) ))
        rot = ((np.pi / 2.0) - theta) / 2.0
        
        # rotate each vector -- below `rot` is deg ccw wrt x-axis in the ref
        # frame of '+axis'
        Rs = utils.ER_rotation_matrix(axis,  rot)
        Rf = utils.ER_rotation_matrix(axis, -rot)
        
        s = np.dot(Rs, s)
        f = np.dot(Rf, f)
        
        assert np.abs( np.dot(s, f) ) < 1e-10
        
        # return a tuple of basis grid objects
        bgs = ( (p0.copy(), s.copy(), f.copy(), shape), 
                (p1.copy(), s.copy(), f.copy(), shape) )
                
        
        return bgs
    
        
    def get_param(self, param_name):
        """
        Parameter getter function.
        """
        if param_name in self._param_list:
            if hasattr(self, param_name):
                return getattr(self, param_name)
            else:
                raise AttributeError('Current CSPad object does not have '
                                     'parameter: %s instantiated' % param_name)
        else:
            raise ValueError('No known parameter with name: %s, (check input file)' % param_name)
    
            
    def set_param(self, param_name, value):
        """
        Parameter setter.
        """
        if hasattr(self, param_name):
            if value.shape == _array_sizes[param_name]:
                if param_name == 'quad_offset_bydiag':
                    self._set_offset_bydiag(value)
                else:
                    self.__dict__[param_name] = value
            else:
                raise ValueError('`value` has wrong shape for: %s' % param_name)
        else:
            raise AttributeError('No parameter with name: %s' % param_name)
    
    
    def set_many_params(self, param_names, param_values):
        """
        Set many params at once.
        
        Here, param_names, param_values are list.
        """
        if (type(param_names) == list) and (type(param_values) == list):
            if len(param_names) == len(param_values):
                for i,pn in enumerate(param_names):
                    self.set_param(pn, param_values[i])
            else:
                raise ValueError('`param_names` & `param_values` must be same len')
        else:
            raise TypeError('`param_names` & `param_values` must be type list')
    
        
    @property
    def quad_offset_bydiag(self):
        """
        This is a parameterization of `self.quad_offset` in terms of x+y, x-y
        instead of x/y. Useful for the optimization routine.
        """
        quad_offset_bydiag = np.zeros_like(self.quad_offset)
        quad_offset_bydiag[0,:] = self.quad_offset[0,:] + self.quad_offset[1,:]
        quad_offset_bydiag[1,:] = self.quad_offset[0,:] - self.quad_offset[1,:]
        return quad_offset_bydiag
        
        
    def _set_offset_bydiag(self, quad_offset_bydiag):
        self.quad_offset[0,:] = (quad_offset_bydiag[0,:] + quad_offset_bydiag[1,:]) / 2.0
        self.quad_offset[1,:] = (quad_offset_bydiag[0,:] - quad_offset_bydiag[1,:]) / 2.0
        return
    

    @property
    def basis_repr(self):
        return self._generate_positional_basis()
        
    
    @property
    def pixel_positions(self):
        """
        Compute and return the x,y,z positions of each pixel. The format is
        shape-(3, 4, 16, 185, 194), indexed as: (x/y/z, quad, ASIC, slow, fast)
        """
        
        bg = self.basis_repr
        pix_pos = np.zeros((3, 4, 16, 185, 194))
        for i in range(4):
            for j in range(16):
                grid = bg.grid_as_explicit(i*16 + j)
                
                if grid.shape == (194,185,3):
                    for k in range(3):
                        pix_pos[k,i,j,:,:] = grid[:,:,k].T
                else:
                    for k in range(3):
                        pix_pos[k,i,j,:,:] = grid[:,:,k]
        
        return pix_pos
    
        
    def _asic_index(self, quad, two_by_one, asic):
        """
        Given a quad index, two_by_one index, asic index, all ints, return the 
        absolute index of the basis grid.
        """
        if quad > 4:
            raise ValueError('There are only 4 quads! Asked for quad: %d' % quad)
        if two_by_one > 8:
            raise ValueError('There are only 8 2x1s per quad! Asked for 2x1: %d' % two_by_one)
        if asic > 2:
            raise ValueError('There are only 2 ASICs per 2x1! Asked for ASIC: %d' % asic)
        return int(quad * 16 + two_by_one * 2 + asic)
    
        
    @property
    def do_asics_overlap(self):
        """
        Returns `True` if the CSPAD has two+ asics that overlap when projected
        into the x-y plane. Else returns `False`.
        """
        
        bg = self.basis_repr
        
        # compute an array of the corners of each ASIC
        corners = np.zeros(( 64, 4, 2 ))
        
        for i in range(64):
            corners[i,:,:] = bg.get_grid_corners(i)[:,:2]
        
        # loop over each ASIC and make sure that the corners of the others
        # are not inside the area it takes up in the xy plane
        
        asics_overlap = False
        
        for i in range(64):
            asic   = corners[i,:,:]
            others = np.concatenate((corners[:i], corners[i+1:])).reshape(63*4, 2)
            any_inside = np.sum( points_inside_poly(others, asic) ).astype(np.bool)
            if any_inside:
                asics_overlap = True                
        
        return asics_overlap
    
        
    def intensity_profile(self, raw_image, n_bins=None, quad='all'):
        """
        Bin pixel intensities by their radius.

        Parameters
        ----------            
        raw_image : np.ndarray
            The intensity at each pixel, same shape as pixel_pos
        n_bins : int
            The number of bins to employ. If `None` guesses a good value.
        quad : int
            Bin only for a single quad. "all" means all quads.

        Returns
        -------
        bin_centers : ndarray, float
            The radial center of each bin.

        bin_values : ndarray, int
            The average intensity in the bin.
        """
        
        # compute radii
        pp = self.pixel_positions
        
        if quad == 'all':
            radii = np.sqrt( np.power(pp[0],2) + np.power(pp[1],2) )
            intensities = raw_image
        elif type(quad) == int:
            radii = np.sqrt( np.power(pp[0,quad], 2) + np.power(pp[1,quad], 2) )
            intensities = raw_image[quad,:,:,:]
        else:
            raise ValueError('`quad` must be {0,1,2,3} or "all", got %s' % str(quad))

        # histogram -- we have two methods : one for binary images
        if intensities.dtype == np.bool:
            
            bin_factor = 100.0 # HARDCODED -- seems to work well
            bin_assignments = np.floor( radii * bin_factor ).astype(np.int32)
        
            assert bin_assignments.shape == radii.shape
            assert intensities.shape     == radii.shape
            
            bin_values  = np.bincount( bin_assignments[intensities].flatten() ).astype(np.float32)
            bin_values /= (np.bincount( bin_assignments.flatten() ) + 1e-100).astype(np.float)[:bin_values.shape[0]]
            bin_centers = np.arange(bin_values.shape[0]) / bin_factor
            
        else:
            if n_bins == None : n_bins = int( np.sqrt(np.product(raw_image.shape)) )
            
            # New algorithm by Jonas to calculate angular average instead of angular sum
            
            bin_values, bin_edges = np.histogram( radii, weights=intensities, bins=n_bins )
            bin_normalizations = np.histogram( radii, bins=n_bins )
            
            bin_values = bin_values / (bin_normalizations[0] + 1e-300)
            bin_centers = np.array([(bin_edges[i] + bin_edges[i+1])/2 for i in range(len(bin_values))])
        
        assert bin_centers.shape == bin_values.shape
        return bin_centers, bin_values
    
        
    def _rotate_xy(self, vector, degrees_ccw):
        """
        Perform a rotation in the x-y plane of a 3-vector
        """
        ct = np.cos( np.deg2rad(degrees_ccw) )
        st = np.sin( np.deg2rad(degrees_ccw) )
        R = np.array([[ct, -st], [st, ct]])
        vector[:2] = np.dot(R, vector[:2])
        return vector
    

    def _generate_positional_basis(self):
        """
        Generate a BasisGrid representation of the CSPad, with all offsets,
        tilts, that dictate the relative positions of the quads applied.
        """

        bg = BasisGrid()
        
        # assemble each 2x1 on the quad
        for quad_index in range(4): # quad in CSPAD
            
            for i in range(8): # 2x1 in quad
                
                for j in range(2): # ASIC in 2x1
                
                    p, s, f, shape = self._metrology_basis.get_grid( self._asic_index(quad_index, i, j) )
                    
                    # replace references w/objects so we dont modify the metrology
                    p = p.copy()
                    s = s.copy()
                    f = f.copy()
                    
                    # in the metrology, each quad is not oriented wrt to one
                    # another -- here we rotate each quad by the appropriate 
                    # amount to get them into the correct relative orientation
                    
                    s = self._rotate_xy( s, self._base_quad_rotation[quad_index] + self.quad_rotation[quad_index])
                    f = self._rotate_xy( f, self._base_quad_rotation[quad_index] + self.quad_rotation[quad_index])
                    p = self._rotate_xy( p, self._base_quad_rotation[quad_index] + self.quad_rotation[quad_index])
                    
                    
                    # we must mirror the x-coordinates of each vector to be consistent with
                    # the CXI coordinate convention, where the x-axis is positive towards
                    # the hutch door (right handed system)
                    if self._reflect_xaxis:
                        p[0] = -p[0]
                        s[0] = -s[0]
                        f[0] = -f[0]
                    
                    
                    # add the quad offset, which defines the relative spatial
                    # orientations of each quad
                    p[:2] += self.quad_offset[quad_index]
                
                    # finally, add these to our basis grid
                    bg.add_grid(p, s, f, shape)
                
        return bg
    
        
    def _asic_rotation(self, quad_index, asic_index):
        """
        Return the xy-rotation of the asic *within the frame of the quad*. Does
        this by comparing the fast-scan vector to [1,0].
        
        Parameters
        ----------
        quad_index : int
            The quad index: 0,1,2,3
            
        asic_index : int
            The index of the asic: 0,1,...,15
            
        Returns
        -------
        theta : float
            The rotation, in degrees.
        """
        
        # determine the rotation
        i = asic_index / 2
        j = asic_index % 2
        p, s, f, shape = self._metrology_basis.get_grid( self._asic_index(quad_index, i, j) )
        theta = utils.arctan3(f[1], f[0]) * (360. / (np.pi * 2.0))
        
        # remove what the default is, due to the CSPad geometry,
        # after our manipulations
        if i in [0,1]:
            base = 0.0
        elif i in [2,3,6,7]:
            base = 270.0
        elif i in [4,5]:
            base = 180.0
            
        theta -= base
        
        return theta
    

    def _twobyone_location(self, quad_index, two_by_one_index):
        """
        Return the lower-left corner of the asic *within the frame of the quad*.
        
        Parameters
        ----------
        quad_index : int
            The quad index: 0,1,2,3
            
        two_by_one_index : int
            The index of the 2x1: 0,1,...,7
            
        Returns
        -------
        c : np.ndarray
            A 3-vector indicating where the lower-left corner of the ASIC is.
            Note that here ***x has not been mirrored*** since we are working
            in the frame of the quad metrology.
        """
        
        # always take the first ASIC
        p, s, f, shape = self._metrology_basis.get_grid( self._asic_index(quad_index, two_by_one_index, 0) )
        
        # CXI is the wild west -- no rules, so roll it manually
        if two_by_one_index in [0,1]:
            c = p.copy() + 185 * s.copy()
        elif two_by_one_index in [2,3,6,7]:
            c = p.copy() + 185 * s.copy() + (194*2 + 3) * f.copy()
        elif two_by_one_index in [4,5]:
            c = p.copy() + (194*2 + 3) * f.copy()
        else:
            raise ValueError('`two_by_one_index` must be in [0,7]')
        
        return c
    
        
    def _assemble_quad(self, quad_image, quad_index):
        """
        Assemble the geometry of an individual quad
        """

        assert quad_image.shape == (16,185,194)

        # make the array for this quadrant
        quad_px_size = 850
        quadrant = np.zeros( (quad_px_size, quad_px_size), dtype=quad_image.dtype )

        for i in range(8):

            # assemble the 2x1 -- insert a 3 px gap
            gap = np.zeros( (185,3), dtype=quad_image.dtype )
            two_by_one = np.hstack( (quad_image[i*2,:,:], gap, 
                                    quad_image[i*2+1,:,:]) )

            # re-orient data
            if i in [0,1]:
                two_by_one = two_by_one[::-1,:]
            elif i in [2,3,6,7]:
                two_by_one = two_by_one[::-1,::-1].T
            elif i in [4,5]:
                two_by_one = two_by_one[:,::-1]


            # rotate the 2x1 to be in the correct orientation
            two_by_one = interp.rotate(two_by_one, -self._asic_rotation(quad_index, i),
                                       output=two_by_one.dtype)

            # find the corner
            c = self._twobyone_location(quad_index, i)
            cs = int( c[0] / self.pixel_size )
            rs = int( c[1] / self.pixel_size )
      
            if (rs < 0) or (rs+two_by_one.shape[0] > quad_px_size):
                raise ValueError('rs: out of bounds in rows')
            if (cs < 0) or (cs+two_by_one.shape[1] > quad_px_size):
                raise ValueError('cs: out of bounds in cols')

            quadrant[rs:rs+two_by_one.shape[0],cs:cs+two_by_one.shape[1]] = two_by_one.copy()
       
        return quadrant
           
           
    def _assemble_image(self, raw_image):
        """
        Build each of the four quads, and put them together.
        """
    
        assert raw_image.shape == (4,16,185,194)
        
        # set up the raw image and the assembled template
        raw_image = read.enforce_raw_img_shape(raw_image)
        
        # for some reason, bool types don't work. Make them ints
        if raw_image.dtype == np.bool:
            raw_image = raw_image.astype(np.int32)
        
        bounds = 2*850 + 300 # JAS: total image range is 2000, ensures beam center is at (1000,1000)
        assembled_image = np.zeros((bounds, bounds), dtype=raw_image.dtype)

        # iterate over quads
        for quad_index in range(4):

            quad_index_image = self._assemble_quad( raw_image[quad_index], quad_index )

            qr = [90, 0 , 270, 180]
            
            # this interp method goes CW, so we have to take the negative...
            quad_index_image = interp.rotate( quad_index_image, 
                                              -(self._base_quad_rotation[quad_index] + self.quad_rotation[quad_index]),
                                              reshape=False,
                                              output=quad_index_image.dtype )
                                  
            # shift the quads into their respective places
            base_row = [850,   850,   0,   0]
            base_col = [  0,   850, 850,   0]
            qoff_row = int(  self.quad_offset[quad_index,1] / self.pixel_size) + \
                             base_row[quad_index] + 150
            qoff_col = int( -self.quad_offset[quad_index,0] / self.pixel_size) + \
                             base_col[quad_index] + 150
                        
            if (qoff_row < 0) or (qoff_row >= bounds):
                raise ValueError('qoff_row: %d out of bounds [0,%d)' % (qoff_row, bounds))
            if (qoff_col < 0) or (qoff_col >= bounds):
                raise ValueError('qoff_col: %d out of bounds [0,%d)' % (qoff_col, bounds))
            
            assembled_image[qoff_row:qoff_row+850, qoff_col:qoff_col+850] = quad_index_image[:,:]
        
        # swap x-axis to conform to CXI convention
        assembled_image = assembled_image[:,::-1]

        return assembled_image
        
        
    # Here we have a series of import/export functions that facilitate
    # conversion between different formats
        
    @classmethod
    def default(cls):
        metrology = np.array([default.q0, default.q1, default.q2, default.q3])
        return cls(metrology, verbose=False)
    
        
    @classmethod
    def load_cspad(cls, filename):
        """
        Loads the a CSPad from disk.

        Parameters
        ----------
        filename : str
            The path to the shotset file.

        Returns
        -------
        cspad : pypad.CSPad
            A CSPad object
        """

        if not filename.endswith('.cspad'):
            raise ValueError('Must load a cspad file (.cspad extension)')

        if not os.path.exists(filename):
            raise IOError('File: %s not found.' % filename)

        hdf = h5py.File(filename, 'r')
        c = cls._from_serial(hdf['/cspad'])
        hdf.close()
        
        return c
    
        
    @classmethod
    def load_cheetah(cls, filename):
        """
        Convert a cheetah pixelmap to a CSPad object.

        Parameters
        ----------
        filename : str
            The path to the cheetah pixelmap.

        Returns
        -------
        cspad : CSPad
            An CSPad object.
        """

        f = h5py.File(filename)

        if not f.keys() == ['x', 'y', 'z']:
            raise IOError('File: %s is not a valid pixel map, should contain fields'
                          ' ["x", "y", "z"] exlusively' % filename)

        # convert m --> mm
        x = read.enforce_raw_img_shape( np.array(f['x']) * 1000.0 )
        y = read.enforce_raw_img_shape( np.array(f['y']) * 1000.0 )

        # for some reason z is in microns, so um --> mm
        z = read.enforce_raw_img_shape( np.array(f['z']) / 1000.0 )

        f.close()

        bg = BasisGrid()
        shape = (185, 194) # will always be this for each ASIC

        # loop over each ASIC, and convert it into a basis grid
        for i in range(4):
            for j in range(16):

                # extract all the corner positions (code ineligant but explicit)
                # corners are numbered 0 -> 4, starting top left and continuing cw
                corners = np.zeros(( 4, 3 ))
                corners[0,:] = ( x[i,j,0,0],   y[i,j,0,0],   z[i,j,0,0]   )
                corners[1,:] = ( x[i,j,0,-1],  y[i,j,0,-1],  z[i,j,0,-1]  )
                corners[2,:] = ( x[i,j,-1,-1], y[i,j,-1,-1], z[i,j,-1,-1] )
                corners[3,:] = ( x[i,j,-1,0],  y[i,j,-1,0],  z[i,j,-1,0]  )

                # dont use this -- takes the vectors from 2 px inside the ASIC
                # if your starting pixel map is good, this should give the same
                # result as the above
                # corners[0,:] = ( x[i,j,2,2],   y[i,j,2,2],   z[i,j,2,2]   )
                # corners[1,:] = ( x[i,j,2,-3],  y[i,j,2,-3],  z[i,j,2,-3]  )
                # corners[2,:] = ( x[i,j,-3,-3], y[i,j,-3,-3], z[i,j,-3,-3] )
                # corners[3,:] = ( x[i,j,-3,2],  y[i,j,-3,2],  z[i,j,-3,2]  )


                # average the vectors formed by the corners to find f/s vects
                # the fast scan direction is the last index, s is next
                # f points left -> right, s points bottom -> top
                f = (( corners[1,:] - corners[0,:] ) + ( corners[2,:] - corners[3,:] ))
                s = (( corners[3,:] - corners[0,:] ) + ( corners[2,:] - corners[1,:] ))

                # make them pixel-size magnitude
                f = f * (pixel_size / np.linalg.norm(f))
                s = s * (pixel_size / np.linalg.norm(s))

                # center is just the average of the 4 corners
                c = np.mean( corners, axis=0 )
                assert c.shape == (3,)
                # c[2] += distance_offset

                bg.add_grid_using_center(c, s, f, shape)

        return cls(bg)
    
        
    @classmethod
    def load_dtc(cls, filename):
        """
        Load an Odin detector into a CSPad object

        Parameters
        ----------
        filename : str
            The path to the Detector object on disk.

        Returns
        -------
        cspad : CSPad
            An CSPad object.
        """
        
        try:
            from odin import xray
        except ImportError as e:
            raise ImportError('Cannot find Odin. You must have Odin installed to '
                              'export to Odin. Download and install Odin from '
                              'https://github.com/tjlane/odin')
        
        dtc = xray.Detector.load(filename)
        
        # we have to convert from an odin-type basis grid
        pypad_bg = BasisGrid()
        odin_bg  = dtc._basis_grid
        
        for i in range(odin_bg.num_grids):
            p, s, f, shp = odin_bg.get_grid(i)
            pypad_bg.add_grid(p, s, f, shp)
            
        return cls(pypad_bg)
    
        
    @classmethod
    def load_crystfel(cls, filename, verbose=False):
        """
        Convert a CrystFEL geom file to a CSPad object.

        Parameters
        ----------
        filename : str
            The path to the geom text file.

        Returns
        -------
        cspad : CSPad
            An CSPad object.
        """
        
        # NOTE ON UNITS: all CrystFEL units are pixel units, except the
        # sample-to-detector offset
        
        if not filename.endswith('.geom'):
            raise IOError('Can only read flat text files with extension `.geom`.'
                          ' Got: %s' % filename)

        if verbose: print "Converting CSPAD geometry in: %s ..." % filename
        f = open(filename, 'r')
        geom_txt = f.read()
        f.close()

        # initialize an odin BasisGrid object -- will add ASICs to this
        bg = BasisGrid()
        shp = (185, 194) # this never changes for the CSPAD

        # now we'll read out the geometry -- do this all in a try block so we can
        # throw an error if something is not there
        try:

            # measure the absolute detector offset
            # right now this appears to be the only z-information in the geometry...
            re_pz = re.search('coffset = (\d+.\d+..\d+)', geom_txt)
            p_z = float(re_pz.group(1)) / 1000.0 # m --> mm

            # iterate over each quad / ASIC
            for q in range(4):
                for a in range(16):

                    if verbose:
                        print "Reading geometry for: QUAD %d / ASIC %d" % (q, a)

                    # match f/s vectors
                    re_fs = re.search('q%da%d/fs =\s+(.\d+.\d+)x\s+(.\d+.\d+)y' % (q, a), geom_txt)
                    f_x = float( re_fs.group(1) )
                    f_y = float( re_fs.group(2) )
                    f = np.array([f_x, f_y, 0.0])
                    f = f * (pixel_size / np.linalg.norm(f))

                    re_ss = re.search('q%da%d/ss =\s+(.\d+.\d+)x\s+(.\d+.\d+)y' % (q, a), geom_txt)
                    s_x = float( re_ss.group(1) )
                    s_y = float( re_ss.group(2) )
                    s = np.array([s_x, s_y, 0.0])
                    s = s * (pixel_size / np.linalg.norm(s))

                    # match corner postions, that become the p vector
                    # note we have to convert from pixel units to mm
                    # and also that CrystFEL measures the corner from the actual
                    # *corner*, and not the center of the corner pixel!
                    re_cx = re.search('q%da%d/corner_x = (.\d+.\d+)' % (q, a), geom_txt)
                    p_x = (float( re_cx.group(1) ) + 0.5) * pixel_size

                    re_cy = re.search('q%da%d/corner_y = (.\d+.\d+)' % (q, a), geom_txt)
                    p_y = (float( re_cy.group(1) ) + 0.5) * pixel_size

                    p = np.array([p_x, p_y, p_z])

                    # finally, add the ASIC to the odin basis grid
                    bg.add_grid(p, s, f, shp)


        # if an Attr error gets thrown, it's because a re failed to match
        except AttributeError as e:
            print e
            raise IOError('Geometry file incomplete -- could not find one or more '
                          'required fields.')

        if verbose: print " ... successfully converted geometry."
        
        return cls(bg)
    

    def to_cheetah(self, filename):
        export.to_cheetah(self, filename)
        return
    
        
    def to_odin(self, energy, distance_offset, filename):
        export.to_odin(self, energy, distance_offset, filename)
        return
    
        
    def to_text(self, filename):
        export.to_text(self, filename)
        return
    
        
    # Below: save/load methods. The saving here is quite lazy on my part, and
    # will be totally incomprehensable to anyone who tries to read the file...
        
    def _to_serial(self):
        """ serialize the object to an array """
        s = np.array( cPickle.dumps(self) )
        s.shape=(1,) # a bit nasty...
        return s
    

    @classmethod
    def _from_serial(self, serialized):
        """ recover a CSPad object from a serialized array """
        if serialized.shape == (1,):
            serialized = serialized[0]
        return cPickle.loads( str(serialized) )
    

    def save(self, filename):
        """
        Writes the current CSPad to disk.

        Parameters
        ----------
        filename : str
            The path to the CSPad file to save.
        """

        if not filename.endswith('.cspad'):
            filename += '.cspad'

        hdf = h5py.File(filename, 'w')
        hdf['/cspad'] = self._to_serial()
        hdf.close()
        
        print 'Wrote %s to disk.' % filename

        return
    

    @classmethod
    def load(cls, filename):
        """
        Load a metrology. Currently supports:
        
            File Type               Extension
            ---------               ---------
            Optical Metrology       .txt
            Serialized CSPad        .cspad
            Odin Detector           .dtc
            Cheetah PixMap          .h5
            CrystFEL Geom           .geom
        
        Parameters
        ----------
        filename : str
            The path to the file to load
            
        Returns
        -------
        cspad : CSPad
            A CSPad object.
        """
        
        if filename.endswith('.txt'):
            geom = cls(filename)
            
        elif filename.endswith('.cspad'):
            geom = cls.load_cspad(filename)
            
        elif filename.endswith('.dtc'):
            geom = cls.load_dtc(filename)
            
        elif filename.endswith('.h5'):
            geom = cls.load_cheetah(filename)
            
        elif filename.endswith('.geom'):
            geom = cls.load_crystfel(filename)
            
        else:
            raise IOError('Cannot understand filetype/extension of: %s' % filename)
            
        return geom


