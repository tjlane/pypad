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

from autogeom import utils

# the expected size of each parameter array
_array_sizes = { 'center' :         (12, 8),
                 'center_corr' :    (12, 8),
                 'common_mode' :    (3,),
                 'filter' :         (3,),
                 'marg_gap_shift' : (3, 4),
                 'offset' :         (3, 4),
                 'offset_corr' :    (3, 4),
                 'pedestals' :      (5920, 388),
                 'pixel_status' :   (5920, 388),
                 'quad_rotation' :  (4,),
                 'quad_tilt' :      (4,),
                 'rotation' :       (4, 8),
                 'tilt' :           (4, 8),
                 'beam_loc' :       (2,) } # this one is our addition


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
    
    def __init__(self, param_dict):
        """
        Initialize an instance of CSPad, corresponding to a single CSPad
        geometry.
        
        Parameters
        ----------
        param_dict : dict
            A dictionary of the parameters that define the CSPad geometry. Each
            key is a string, and each value is an np.ndarray of size indicated
            below. May include:
            
               Key             Value Size
              ------           ----------
             'center' :         (12, 8),
             'center_corr' :    (12, 8),
             'common_mode' :    (3,),
             'filter' :         (3,),
             'marg_gap_shift' : (3, 4),
             'offset' :         (3, 4),
             'offset_corr' :    (3, 4),
             'pedestals' :      (5920, 388),
             'pixel_status' :   (5920, 388),
             'quad_rotation' :  (4,),
             'quad_tilt' :      (4,),
             'rotation' :       (4, 8),
             'tilt' :           (4, 8)
             
        """

        self._param_list = _array_sizes.keys()
        
        self._check_param_dict(param_dict)
        
        # inject each parameter as an attribute into the class
        for k,v in param_dict.items():
            self.__setattr__(k,v)
        
        self._process_parameters()
        self.small_angle_tilt = True # always true for now
                
        return
    
    
    def __call__(self, raw_image):
        """
        Takes a raw image (8, 188, 3xx) and assembles it into a two-dimensional
        view.
        
        Returns
        -------
        assembled_image : ndarray, float
            The assembled image.
        """
        return self._assemble_image(raw_image)
    
    
    def get_param(self, param_name):
        """
        Parameter getter function.
        """
        if param_name in self._param_list:
            return self.__dict__[param_name]
        else:
            raise ValueError('No parameter with name: %s' % param_name)
    
            
    def set_param(self, param_name, value, process=True):
        """
        Parameter setter.
        """
        if param_name in self._param_list:
            if value.shape == _array_sizes[param_name]:
                self.__dict__[param_name] = value
                if process:
                    self._process_parameters()
            else:
                raise ValueError('`value` has wrong shape for: %s' % param_name)
        else:
            raise ValueError('No parameter with name: %s' % param_name)
    
    
    def set_many_params(self, param_names, param_values):
        """
        Set many params at once.
        
        Here, param_names, param_values are list.
        """
        if (type(param_names) == list) and (type(param_values) == list):
            if len(param_names) == len(param_values):
                for i,pn in enumerate(param_names):
                    self.set_param(pn, param_values[i], process=False)
                self._process_parameters()
            else:
                raise ValueError('`param_names` & `param_values` must be same len')
        else:
            raise TypeError('`param_names` & `param_values` must be type list')
    
        
    def _check_param_dict(self, param_dict):
        """
        Does a sanity check on a parameter dictionary.
        """
        
        if not param_dict.keys().sort() == self._param_list.sort():
            raise ValueError('keys of `param_dict` do not match expected parame'
                             'ter names')

        for key in param_dict.keys():
            if not param_dict[key].shape == _array_sizes[key]:
                raise ValueError('value of %s does not match expected shape: %s'
                                  % (key, str(_array_sizes[key])))
                
        return
        
    
    @property
    def basis_repr(self):
        return self._generate_basis()
        
    
    @property
    def pixel_positions(self):
        """
        Compute and return the x,y,z positions of each pixel.
        """
        bg = self._generate_basis()
        xyz = bg.to_explicit()
        return xyz
    
        
    def _process_parameters(self):
        """
        Alignment calibrations as defined for psana. Injects
        
        -- self.tilt_array
        -- self.sec_offset
        -- self.section_centers
        
        into the local namespace.
        
        This method adapted from XtcExplorer's cspad.py
        """
        
        # angle of each section = rotation (nx90 degrees) + tilt (small angles)
        # ... (4 rows (quads) x 8 columns (sections))
        self.rotation_array = self.rotation
        self.tilt_array     = self.tilt
        self.section_angles = self.rotation_array + self.tilt_array
        
        # position of sections in each quadrant (quadrant coordinates)
        # ... (3*4 rows (4 quads, 3 xyz coordinates) x 8 columns (sections))
        centers              = self.center
        center_corrections   = self.center_corr
        self.section_centers = np.reshape( centers + center_corrections, (3,4,8) )
        
        # quadrant offset parameters (w.r.t. image 0,0 in upper left coner)
        # ... (3 rows (xyz) x 4 columns (quads) )
        quad_pos      = self.offset
        quad_pos_corr = self.offset_corr
        quad_position = quad_pos + quad_pos_corr

        # ... (3 rows (x,y,z) x 4 columns (section offset, quad offset, 
        # quad gap, quad shift)
        marg_gap_shift = self.marg_gap_shift

        # break it down (extract each column, make full arrays to be added to 
        # the above ones)
        self.sec_offset = marg_gap_shift[:,0]

        quad_offset = marg_gap_shift[:,1]
        quad_gap    = marg_gap_shift[:,2]
        quad_shift  = marg_gap_shift[:,3]

        # turn them into 2D arrays, use numpy element-wise multiplication
        quad_offset = np.array( [quad_offset,
                                 quad_offset,
                                 quad_offset,
                                 quad_offset] ).T
                                 
        quad_gap = np.array( [quad_gap*[-1,-1,1],
                              quad_gap*[-1,1,1],
                              quad_gap*[1,1,1],
                              quad_gap*[1,-1,1]] ).T
                              
        quad_shift = np.array( [quad_shift*[1,-1,1],
                                quad_shift*[-1,-1,1],
                                quad_shift*[-1,1,1],
                                quad_shift*[1,1,1]] ).T
                                
        self.quad_offset = quad_position + quad_offset + quad_gap + quad_shift
        
        return
    

    def _generate_basis(self):
        """
        Generate a BasisGrid representation of the CSPad 
        """

        bg = BasisGrid()

        shape = (185, 388)                     # slow, fast
        slow = np.array([ 0.10992, 0.0, 0.0 ]) # basis vector for pixels in x
        fast = np.array([ 0.0, 0.10992, 0.0 ]) # basis vector for pixels in y

        # assemble each 2x1 on the quad
        for quad_index in range(4):
            for i in range(8):

                # we have to re-orient each 2x1
                if i==0 or i==1:
                    # reverse slow dim, switch slow/fast
                    f = -slow.copy()
                    s =  fast.copy()
                    shape = shape[::-1]
                if i==4 or i==5:
                    # reverse fast dim, switch slow/fast
                    f =  slow.copy()
                    s = -fast.copy()
                    shape = shape[::-1]
                else:
                    f = fast.copy()
                    s = slow.copy()

                # now, apply a small rotation in the detector (x/y) plane
                if self.small_angle_tilt:
                    ct = np.cos( np.deg2rad(self.tilt_array[quad_index][i] + (90 * (4-quad_index)) ) )
                    st = np.sin( np.deg2rad(self.tilt_array[quad_index][i] + (90 * (4-quad_index)) ) )
                    R = np.array([[ct, -st], [st, ct]])
                    f[:2] = np.dot(R, f[:2])
                    s[:2] = np.dot(R, s[:2])

                # find the center of the 2x1
                cx = self.sec_offset[0] - self.section_centers[0][quad_index][i] + self.quad_offset[0,quad_index]
                cy = self.sec_offset[1] - self.section_centers[1][quad_index][i] + self.quad_offset[1,quad_index]
                cz = self.sec_offset[2] - self.section_centers[2][quad_index][i] + self.quad_offset[2,quad_index]
                center = np.array([cx, cy, cz])
                center *= 0.10992 # convert to microns

                bg.add_grid_using_center(center, s, f, shape)

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

               if self.small_angle_tilt:
                   pair = interp.rotate(pair, self.tilt_array[quad_index][i],
                                        output=pair.dtype)

           # make the array for this quadrant
           quadrant = np.zeros( (850, 850), dtype=raw_image.dtype )

           # insert the 2x1 sections according to
           for sec in range(8):
               nrows, ncols = pairs[sec].shape

               # colp,rowp are where the top-left corner of a section should be placed
               rowp = 850 - self.sec_offset[0] - (self.section_centers[0][quad_index][sec] + nrows/2)
               colp = 850 - self.sec_offset[1] - (self.section_centers[1][quad_index][sec] + ncols/2)

               quadrant[rowp:rowp+nrows, colp:colp+ncols] = pairs[sec][0:nrows,0:ncols]

           return quadrant
           
           
    def _assemble_image(self, raw_image):
        """
        Build each of the four quads, and put them together.
        """
        
        raw_image = utils.enforce_raw_img_shape(raw_image)
        
        assembled_image = np.zeros((2*850+100, 2*850+100), dtype=raw_image.dtype)
        
        # iterate over quads
        for quad_index in range(4):

            quad_index_image = self._assemble_quad( raw_image[quad_index], quad_index )

            if quad_index>0:
                # reorient the quad_index_image as needed
                quad_index_image = np.rot90( quad_index_image, 4-quad_index )

            qoff_x = self.quad_offset[0,quad_index]
            qoff_y = self.quad_offset[1,quad_index]
            assembled_image[qoff_x:qoff_x+850, qoff_y:qoff_y+850] = quad_index_image

        return assembled_image
        
        
    def to_dir(self, dir_name, run_range=None):
        """
        Write the parameters to `dir_name`, in the standard psana/pyana format.
        
        This method will create a 2-level directory tree under `dir_name` to
        hold these parameters.
        
        Parameters
        ----------
        dir_name : str
            The name of the parent dir (parameter set) to place these parameters
            under.
            
        run_range : tuple
            A tuple of values (X,Y), with X,Y ints, such that these parameters
            will be used for all runs between and including X to Y. If `None`
            (default), then gets set to 0-end.data, meaning all runs.
        """
    
        if os.path.exists(dir_name):
            print "WARNING: rm'ing %s..." % dir_name
            os.system('rm -r %s' % dir_name)
    
        os.mkdir(dir_name)
    
        if run_range == None:
            param_filenames = '0-end.data'
        else:
            param_filenames = '%d-%d.data' % run_range
    
        for key in _array_sizes.keys():
            os.mkdir( pjoin( dir_name, key ))
            fname = pjoin( dir_name, key, param_filenames )
            
            # the file-format for these two parameters is slightly different
            # from all the others... try to maintain that.
            if key == 'center' or key == 'center_corr':
            
                v = self.get_param(key)
                v = v.reshape(3,4,8)
            
                f = open(fname, 'w')
                for i in range(3):
                    for j in range(4):
                        f.write( '\t'.join([ str(x) for x in v[i,j,:] ]) + '\n' )
                    f.write('\n')
            
            else:
                np.savetxt(fname, self.get_param(key), fmt='%.2f')
            print "Wrote: %s" % fname
    
        return
    
    @classmethod
    def from_dir(cls, path, run_number=0):
        """
        Load a parameter set from disk.
        
        Parameters
        ----------
        path : str
            The path to the directory containing the parameters.
            
        run_number : int
            Load the parameters for this run.
        """
        
        param_dict = {}
        
        # if not to_load:
        #     to_load = ['center', 'center_corr', 'filter', 'marg_gap_shift', 
        #                'offset', 'offset_corr', 'pixel_status', 'quad_rotation',
        #                'quad_tilt',  'rotation', 'tilt'] # 'common_mode', 'pedestals'
    
        for p in _array_sizes.keys():
            
            # scan through the possible files and find the right one matching
            # our run number
            files = glob( pjoin(path, p, '*') )
            
            # if theres nothing in the dir, complain
            # we have to deal with beam_loc (aka abs_center) separately, since
            # this has been introduced into the cspad geometry by us...
            if len(files) == 0:
                if p == 'beam_loc':
                    print "Could not locate dir 'beam_loc', using default value"
                    param_dict[p] = np.array([900.0, 870.0]) # default value
                    continue # skip the rest of whats below
                else:
                    raise IOError('No files in parameter dir: %s' % pjoin(path, p))
            
            filename = None
            for f in files:
                
                start, end = os.path.basename(f).split('.')[0].split('-')
                start = int(start)
                if end == 'end':
                    end = 9999
                else:
                    end = int(end)
                    
                if (run_number >= start) and (run_number <= end):
                    filename = f
            
            if filename:
                print "Loading parameters in:", filename
                param_dict[p] = np.loadtxt(filename).reshape(_array_sizes[p])
            else:
                raise IOError('Could not find file for run %d in %s' % (run_number, str(files)))
    
        return cls(param_dict)
        

class Metrology(object):
    """
    Class for interacting with a CSPad optical metrology.
    """
    
    def __init__(self, metrology_file):
        """
        Generate an instance of a CSPad metrology.
        
        Parameters
        ----------
        metrology_file : str
            A file containing the optical metrology readout, formatted in LCLS's
            custom format. See Confluence for more info (good luck...).
        """
        self.load(metrology_file)
        
    
    @property
    def _flat_pixel_positions(self):
        """
        Returns all the pixels as an (x,y,z) / N x 3 array.
        """
        return np.vstack(( self.x_coordinates.flatten(),
                           self.y_coordinates.flatten(),
                           self.z_coordinates.flatten() )).T
    
        
    def load(self, metrology_file, verbose=False):
        """
        Make coordinate maps from meterology file
        """
        
        self.x_coordinates = np.zeros((4,8,185,388), dtype="float")
        self.y_coordinates = np.zeros((4,8,185,388), dtype="float")
        self.z_coordinates = np.zeros((4,8,185,388), dtype="float")


        def get_asics(bigsection):
            """Utility function"""
            asic0 = bigsection[:,0:194]
            asic1 = bigsection[:,(391-194):]
            asics = np.concatenate( (asic0,asic1), axis=1 )
            return asics
        

        # section pixel array / grid
        rr,cc = np.mgrid[0:185:185j, 0:391:391j]

        # now compute the "fractional pixels"
        rrfrac = rr / 185.0
        ccfrac = cc / 391.0

        # remove the 3-pixel gap
        rrfrac = get_asics(rrfrac)
        ccfrac = get_asics(ccfrac)

        sec_coords = np.array([rrfrac,ccfrac])
        sec_coord_order = [ (1,2,0,3), (1,2,0,3), (2,3,1,0), (2,3,1,0),
                            (3,0,2,1), (3,0,2,1), (2,3,1,0), (2,3,1,0) ]

        # load data from metrology file (ignore first column)
        metrology = np.loadtxt(metrology_file)[:,1:]
        metrology = metrology.reshape(4,8,4,3)

        # also, we need to resort the 2x1 sections, they are
        # listed in the file in the order 1,0,3,2,4,5,7,6
        metrology = metrology[:,(1,0,3,2,4,5,7,6),:,:]

        dLong = np.zeros((4,8,2), dtype="float64")
        dShort = np.zeros((4,8,2), dtype="float64")
        for quad in range(4):

            for sec in range(8):

                # corner positions (in micrometers)
                input_x = metrology[quad,sec,sec_coord_order[sec],0].reshape(2,2)
                input_y = metrology[quad,sec,sec_coord_order[sec],1].reshape(2,2)
                input_z = metrology[quad,sec,sec_coord_order[sec],2].reshape(2,2)

                # interpolate coordinates over to the pixel map
                self.x_coordinates[quad,sec] = interp.map_coordinates(input_x, sec_coords)
                self.y_coordinates[quad,sec] = interp.map_coordinates(input_y, sec_coords)
                self.z_coordinates[quad,sec] = interp.map_coordinates(input_z, sec_coords)

                # ! in micrometers! Need to convert to pixel units
                dL = np.array([ abs(input_x[0,1]-input_x[0,0])/391, 
                                abs(input_x[1,1]-input_x[1,0])/391,
                                abs(input_y[0,0]-input_y[0,1])/391,
                                abs(input_y[1,0]-input_y[1,1])/391 ])
                dLong[quad,sec] = dL[dL>100] # filter out the nonsense ones

                dS = np.array([ abs(input_y[0,0]-input_y[1,0])/185,
                                abs(input_y[0,1]-input_y[1,1])/185, 
                                abs(input_x[0,0]-input_x[1,0])/185,
                                abs(input_x[0,1]-input_x[1,1])/185 ])
                dShort[quad,sec] = dS[dS>100] # filter out the nonsense ones

        dTotal = np.concatenate( (dLong.ravel(), dShort.ravel() ))
        
        if verbose:
            print "Pixel-size:"
            print "     long side average:    %.2f +- %.2f "%( dLong.mean(), dLong.std())
            print "     short side average:   %.2f +- %.2f "%( dShort.mean(), dShort.std())
            print "     all sides average:    %.2f +- %.2f "%( dTotal.mean(), dTotal.std())

        # use the total to convert it all to pixel units
        self.x_coordinates = self.x_coordinates / dTotal.mean()
        self.y_coordinates = self.y_coordinates / dTotal.mean()
        self.z_coordinates = self.z_coordinates / dTotal.mean()

        origin = np.array([ [834.0, 834.0], [834.0, 834.0], 
                            [834.0, 834.0], [834.0, 834.0] ])
        
        for quad in range(4):
            # For each quad, rotate and shift into the image coordinate system
            if quad==0 :
                savex = np.array( self.x_coordinates[quad] )
                self.x_coordinates[quad] = origin[quad][0] - self.y_coordinates[quad]
                self.y_coordinates[quad] = origin[quad][1] - savex
            if quad==1 :
                self.x_coordinates[quad] = origin[quad][0] + self.x_coordinates[quad]
                self.y_coordinates[quad] = origin[quad][1] - self.y_coordinates[quad]
            if quad==2 :
                savex = np.array( self.x_coordinates[quad] )                
                self.x_coordinates[quad] = origin[quad][0] + self.y_coordinates[quad]
                self.y_coordinates[quad] = origin[quad][1] + savex
            if quad==3 :
                self.x_coordinates[quad] = origin[quad][0] - self.x_coordinates[quad]
                self.y_coordinates[quad] = origin[quad][1] + self.y_coordinates[quad]        
        
        return
        
        
    def _apply_geom_corrections(self):
        """
        Apply specific geometry corrections, defined by the psana parameters,
        to the metrology pixel map.
        
        Returns
        -------
        x, y, z : np.ndarrays of shape (4,8,185,388)
            The x,y,z locations of each pixel in real-space, corrected.
        """
        
        x = self.x_coordinates.copy()
        y = self.y_coordinates.copy()
        z = self.z_coordinates.copy()
        
        for i in range(4):
            x[i,:,:,:] += self.offset_coor[0,i]
            y[i,:,:,:] += self.offset_coor[1,i]
            z[i,:,:,:] += self.offset_coor[2,i]
        
        return x, y, z
    
        
class BasisGrid(object):
    """
    A class representing a set of rectangular grids in space -- specifically,
    x-ray scattering detectors.
    
    Note that the geometry below is definied in "slow" and "fast" scan 
    dimensions. These are simply the two dimensions that define the plane
    a single rectangular pixel grid lives in. They may also be called the x and
    y dimensions without any loss of generality.
    
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
        
        self.num_grids  = 0
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
        self.num_grids += 1
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
        pix_size = (np.linalg.norm(s), np.linalg.norm(f))
        center = (np.array(shape) - 1) * np.array(pix_size) / 2.
        assert center.shape == (2,)
        p = p_center.copy()
        p[:2] = p[:2] - center
        
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
    
        
    def to_explicit(self):
        """
        Return the entire grid as an n x 3 array, defining the x,y,z positions
        of each pixel.
                
        Returns
        -------
        xyz : np.ndarray, float
            An N x 3 array of the x,y,z positions of each pixel
        """
        xyz = np.concatenate([ self.grid_as_explicit(i) for i in range(self.num_grids) ])
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
            An N x 3 array of the x,y,z positions of each pixel
        """
        
        p, s, f, shape = self.get_grid(grid_number)
        
        # xyz = i * s + j * f, where i,j are ints running over range `shape`
        mg = np.mgrid[0:shape[0]-1:1j*shape[0], 0:shape[1]-1:1j*shape[1]]
        xyz = np.outer(mg[0].flatten(), s) + np.outer(mg[1].flatten(), f)
        
        # translate
        xyz += p
        
        return xyz
