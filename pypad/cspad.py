#!/usr/bin/env python

"""
An interface to the CSPad geometry. Converts the pyana/psana parameters that
definte the CSPad geometry into actual images or a pixel map of the positions
of each pixel.
"""

import sys
import os
from glob import glob
from os.path import join as pjoin

import numpy as np
import scipy.ndimage.interpolation as interp
import matplotlib.pyplot as plt

from pypad import utils
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



_array_sizes = {'quad_rotation' : (4,),
                'quad_offset'   : (4,2)}

class CSPad(object):
    """
    This is a container class for saving, loading, and interacting with the 
    parameter set that defines the CSPad geometry.
    
    The main idea here is that you can instantiate a CSPad instance that is
    defined by a set of parameters, and use that instance to display images
    collected on the CSPad:
    
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
    
    >>> geom = CSPad.from_dir('my_params', run=0)
    >>> assembled_image = geom(raw_image)
    >>> imshow(assembled_image)
    
    Of course, the 'geom' object can be re-used to assemble many images.
    
    This class is largely based on XtcExplorer's cspad.py
    """
    
    def __init__(self, metrology_dir, quad_offset=np.zeros((4,2)), quad_rotation=np.zeros(4)):
        """
        Initialize an instance of CSPad, corresponding to a single CSPad
        geometry.             
        """
        
        self._param_list = _array_sizes.keys()
        self.pixel_size  = 0.10992
        
        if not quad_offset.shape == _array_sizes['quad_offset']:
            raise ValueError('quad_offset must have shape (4,2), got: %s' % str(quad_offset.shape))
        if not quad_rotation.shape == _array_sizes['quad_rotation']:
            raise ValueError('quad_rotation must have shape (4,), got: %s' % str(quad_rotation.shape))
        
        self.quad_offset   = quad_offset
        self.quad_rotation = quad_rotation
        
        # read the optical metrology
        self.metrology_basis = [[],[],[],[]]
        
        qms = self._read_metrology(metrology_dir)
        for q in range(4):
            for two_by_one in range(8):
                self.metrology_basis[q].append( self._twobyone_to_bg(qms[q], q, two_by_one) )
        
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
        return self._assemble_image(raw_image)
    
        
    @staticmethod
    def _unit(v):
        return v / np.linalg.norm(v)
        
        
    def _read_metrology(self, metrology_dir):
        """
        Tmp file parser.
        """
        
        quad_metrologies = []
        for i in range(4):
            print "loading " + metrology_dir + '/q%d.txt' % i
            qm = np.genfromtxt(metrology_dir+'/q%d.txt' % i)
            quad_metrologies.append( qm.copy() )
            
        return quad_metrologies


    def _twobyone_to_bg(self, quad_metrology, quad_index, two_by_one_index):
        """
        Convert a 2x1 in the optical metrology into two basis grid elements.

        Parameters
        ----------
        quad_metrology : np.ndarray, float
            A 32 x 3 array of the corner positions of each 2x1 on the quad.

        two_by_one_index : int
            The index of the 2x1 to check out.
            
        Returns
        -------
        bgs : tuple
            A 2-tuple of (p,s,f,shape) for each ASIC in the two by one.
        """

        if not quad_metrology.shape == (32,3):
            raise ValueError('Invalid quad_metrology, must be shape (32,3), got: %s' % str(quad_metrology.shape))

        # below is the sequence in which each 2x1 is optically measured. TJL todo :
        # this could be automagically found
        scan_sequence = np.array([1,0,3,2,4,5,7,6])
        shape = (185, 194) # always the same, for each ASIC

        # rip out the correct corner positions for the 2x1 we want
        i = int(np.where( scan_sequence == two_by_one_index )[0]) * 4
        xyz_2x1 = quad_metrology[i:i+4,:] # the four corners of the 2x1


        # find p, s, f
        if two_by_one_index in [0,1]:

            s = 0.10992 * self._unit(xyz_2x1[0,:] - xyz_2x1[1,:])
            f = 0.10992 * self._unit(xyz_2x1[2,:] - xyz_2x1[1,:])

            p0 = xyz_2x1[1,:] / 1000.
            p1 = p0 + shape[1] * f + 0.27480 * self._unit(f) # for 3px gap


        elif two_by_one_index in [2,3,6,7]:

            s = 0.10992 * self._unit(xyz_2x1[1,:] - xyz_2x1[2,:])
            f = 0.10992 * self._unit(xyz_2x1[3,:] - xyz_2x1[2,:])

            p0 = xyz_2x1[2,:] / 1000.
            p1 = p0 + shape[1] * f + 0.27480 * self._unit(f) # for 3px gap


        elif two_by_one_index in [4,5]:

            s = 0.10992 * self._unit(xyz_2x1[2,:] - xyz_2x1[3,:])
            f = 0.10992 * self._unit(xyz_2x1[0,:] - xyz_2x1[3,:])

            p0 = xyz_2x1[3,:] / 1000.
            p1 = p0 + shape[1] * f + 0.27480 * self._unit(f) # for 3px gap


        else:
            raise ValueError('two_by_one_index must be in 0...7')

        bgs = ( (p0.copy(), s.copy(), f.copy(), shape), (p1.copy(), s.copy(), f.copy(), shape) )

        return bgs
    
        
    def get_param(self, param_name):
        """
        Parameter getter function.
        """
        if param_name in self._param_list:
            return self.__dict__[param_name]
        else:
            raise ValueError('No parameter with name: %s, (check input file)' % param_name)
    
            
    def set_param(self, param_name, value):
        """
        Parameter setter.
        """
        if hasattr(self, param_name):
            if value.shape == _array_sizes[param_name]:
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
                    self.set_param(pn, param_values[i], process=False)
            else:
                raise ValueError('`param_names` & `param_values` must be same len')
        else:
            raise TypeError('`param_names` & `param_values` must be type list')
    
        
    @property
    def basis_repr(self):
        return self._generate_basis()
        
    
    @property
    def pixel_positions(self):
        """
        Compute and return the x,y,z positions of each pixel. The format is
        shape-(3, 4, 16, 185, 194), indexed as: (x/y/z, quad, ASIC, slow, fast)
        """
        
        bg = self._generate_basis()
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
    
    
    def intensity_profile(self, raw_image, n_bins=None, quad='all', beta=10.0,
                          window_size=10):
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
            The total intensity in the bin.
        """

        pixel_pos = self.pixel_positions
        center = self.beam_location * self.pixel_size
        
        if not pixel_pos.shape[1:] == raw_image.shape:
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
        if quad == 'all':
            radii = np.sqrt( np.power(pixel_pos[0]-center[0], 2) + \
                             np.power(pixel_pos[1]-center[1], 2) )
            intensities = raw_image
        elif type(quad) == int:
            radii = np.sqrt( np.power(pixel_pos[0,quad]-center[0], 2) + \
                             np.power(pixel_pos[1,quad]-center[1], 2) )
            intensities = raw_image[quad,:,:,:]
        else:
            raise ValueError('`quad` must be {0,1,2,3} or "all", got %s' % str(quad))


        # generate the histogram
        if n_bins == None:
            n_bins = np.sqrt( np.product(intensities.shape) ) / 2.

        assert radii.shape == intensities.shape

        if intensities.dtype == np.bool:
            bin_values, bin_edges = np.histogram( radii * intensities, bins=n_bins )
        else:
            bin_values, bin_edges = np.histogram( radii, weights=intensities, bins=n_bins )

        bin_values = bin_values[1:]
        bin_centers = bin_edges[1:-1] + np.abs(bin_edges[2] - bin_edges[1])

        bin_values = utils.smooth(bin_values, beta=beta, 
                                  window_size=window_size)
        
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
    

    def _generate_basis(self):
        """
        Generate a BasisGrid representation of the CSPad 
        """

        # The loop over quads below is [1,2,3,0] because
        # the quad positions are supposed to be (as viewed from upstream):
        # 
        #     Q0     Q1
        # 
        #         x
        # 
        #     Q3     Q2
        #
        # -- TJL 28.2.13

        bg = BasisGrid()
        
        # assemble each 2x1 on the quad
        for quad_index in range(4): # quad in CSPAD
            
            for i in range(8): # 2x1 in quad
                
                for j in range(2): # ASIC in 2x1
                
                    p, s, f, shape = self.metrology_basis[quad_index][i][j]
                    
                    # rotate each quad by the appropriate amount ( n * 90-deg CW 
                    # rotations for quads n = { 0, 1, 2, 3 } ), then translate
                    # them so that the top-left corners of an 850 x 850 pixel box
                    # overlap. This remains true to psana convention.
                
                    # perform the rotation
                    qr = [90.0, 0.0, 270.0, 180.0]
                    s = self._rotate_xy( s, qr[quad_index] + self.quad_rotation[quad_index])
                    f = self._rotate_xy( f, qr[quad_index] + self.quad_rotation[quad_index])
                    p = self._rotate_xy( p, qr[quad_index] + self.quad_rotation[quad_index])
                    
                    # add the quad offset, which defines the relative spatial
                    # orientations of each quad
                    p[:2] += self.quad_offset[quad_index]
                
                    # finally, add these to our basis grid
                    bg.add_grid(p.copy(), s.copy(), f.copy(), shape)
                
        return bg
    
        
    def _assemble_quad(self, raw_image, quad_index):
           """
           assemble the geometry of an individual quad
           """
           
           pairs = []
           for i in range(8):

               # 1) insert gap between asics in the 2x1
               asics = np.hsplit( raw_image[i], 2)
               gap = np.zeros( (185,3), dtype=raw_image.dtype )
               
               # gap should be 3 pixels wide
               pair = np.hstack( (asics[0], gap, asics[1]) )

               # all sections are originally 185 (rows) x 388 (columns)
               # Re-orient each section in the quad
               if i==0 or i==1:
                   pair = pair[:,::-1].T   # reverse columns, switch columns to rows.
               if i==4 or i==5:
                   pair = pair[::-1,:].T   # reverse rows, switch rows to columns
               pairs.append( pair )

               pair = interp.rotate(pair, self.quad_rotation[quad_index][i],
                                    output=pair.dtype)

           # make the array for this quadrant
           quadrant = np.zeros( (850, 850), dtype=raw_image.dtype )

           # insert the 2x1 sections according to
           for sec in range(8):
               nrows, ncols = pairs[sec].shape

               # old XTCExplorer code -- NOT CXI convention (Mikhail convention)
               # colp,rowp are where the top-left corner of a section should be placed
               # rowp = 850 - self.sec_offset[0] - (self.section_centers[0][quad_index][sec] + nrows/2)
               # colp = 850 - self.sec_offset[1] - (self.section_centers[1][quad_index][sec] + ncols/2)
               # quadrant[rowp:rowp+nrows,colp:colp+ncols] = pairs[sec][0:nrows,0:ncols]
               
               # new code -- TJL and JAS : adopting CXI coordinate convention
               rowp = int( self.sec_offset[0] + (self.section_centers[0][quad_index][sec] + nrows/2) )
               colp = int( self.sec_offset[1] + (self.section_centers[1][quad_index][sec] + ncols/2) )
               quadrant[rowp-nrows:rowp,colp-ncols:colp] = pairs[sec][::-1,::-1]

           return quadrant
           
           
    def _assemble_image(self, raw_image):
        """
        Build each of the four quads, and put them together.
        """
        
        # shift the quads so that the description of the detector begins at
        # (0,0) in the coordinate system, and is not translated far from the
        # origin
        # TJL : could include beam_vector here too

        for i in range(2): # TJL : could cover z as well, right now just x/y
            min_offset = np.min(self.quad_offset[i,:])
            self.quad_offset[i,:] -= min_offset
            self.beam_location[i] -= min_offset
            
                
        # set up the raw image and the assembled template
        raw_image = utils.enforce_raw_img_shape(raw_image)
        bounds = 2*850+100
        assembled_image = np.zeros((bounds, bounds), dtype=raw_image.dtype)
        
        # iterate over quads
        for quad_index in range(4):

            quad_index_image = self._assemble_quad( raw_image[quad_index], quad_index )

            # reorient the quad_index_image as needed
            #quad_index_image = np.rot90( quad_index_image, 4-quad_index )
                
            # we may want a small correction to the quad rotations. psana
            # has a quad_rotation parameter, that is filled with silly
            # values describing the non-relative orientations of the quads
            # (and the value meanings are ambiguous) -- so here we subtract
            # those values and apply just the small correction
            
            quad_rot_corr_defaults = [180.0, 90.0, 0.0, 270.0]
            quad_rot_corr = self.quad_rotation[quad_index] - quad_rot_corr_defaults[quad_index]
            quad_index_image = interp.rotate( quad_index_image, 
                                              90.0*(4-quad_index) + quad_rot_corr,
                                              reshape=False,
                                              output=quad_index_image.dtype )

            qoff_x = int( self.quad_offset[0,quad_index] )
            qoff_y = int( self.quad_offset[1,quad_index] )
            
            # TJL & JAS -- new CXI convention
            qoff_x = bounds - qoff_x - 850
            qoff_y = bounds - qoff_y - 850
            # -------------------------------
            
            if (qoff_x < 0) or (qoff_x >= bounds):
                raise ValueError('qoff_x: %d out of bounds [0,%d)' % (qoff_x, bounds))
            if (qoff_y < 0) or (qoff_y >= bounds):
                raise ValueError('qoff_y: %d out of bounds [0,%d)' % (qoff_y, bounds))
            
            assembled_image[qoff_x:qoff_x+850, qoff_y:qoff_y+850] = quad_index_image

        return assembled_image
            
    
    def to_cheetah(self, filename="pixelmap-cheetah-raw.h5"):
        export.to_cheetah(self, filename)
        return
        
        
    def to_odin(self, filename):
        export.to_odin(self, filename)
        return
    

