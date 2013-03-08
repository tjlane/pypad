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
from autogeom import default


# the expected size of each parameter array
_array_sizes = { 'center' :         (12, 8),
                 'center_corr' :    (12, 8),
                 'common_mode' :    (3,),
                 'filter' :         (3,),
                 'marg_gap_shift' : (3, 4),
                 'offset' :         (3, 4),
                 'offset_corr' :    (3, 4),
                 'offset_corr_xy' : (2, 4),        # our addition
                 'pedestals' :      (5920, 388),
                 'pixel_status' :   (5920, 388),
                 'quad_rotation' :  (4,),
                 'quad_tilt' :      (4,),
                 'rotation' :       (4, 8),
                 'tilt' :           (4, 8),
                 'beam_location' :  (2,),  # new in psana
                 'beam_vector' :    (3,) } # new in psana
                 

# the parameters that we can expect any psana project to serve up
_psana_params = [ 'center',
                  'center_corr',
                  'common_mode',
                  'filter',
                  'marg_gap_shift',
                  'offset',
                  'offset_corr',
                  'pedestals',
                  'pixel_status',
                  'quad_rotation',
                  'quad_tilt',
                  'rotation',
                  'tilt']


# we have some additions of our own, enumerated here, that we don't expect
# a psana project to have at the beginning
_autogeom_additional_parameters = ['beam_location', 'offset_corr_xy']


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
        self.pixel_size = 0.10992
        
        self._check_params(param_dict)
        
        # inject each parameter as an attribute into the class
        for k,v in param_dict.items():
            self.__setattr__(k,v)
        
        self._process_parameters()
                
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
        if param_name == 'offset_corr_xy':
            if value.shape == _array_sizes[param_name]:
                self.offset_corr[:2,:] = value
        elif param_name in self._param_list:
            if value.shape == _array_sizes[param_name]:
                self.__dict__[param_name] = value
            else:
                raise ValueError('`value` has wrong shape for: %s' % param_name)
                
        else:
            raise ValueError('No parameter with name: %s' % param_name)
        if process:
            self._process_parameters()
    
    
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
    
        
    def _check_params(self, param_dict):
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
        Compute and return the x,y,z positions of each pixel. The format is
        shape-(3, 4, 8, 185, 388), indexed as: (x/y/z, quad, 2x1, slow, fast)
        """
        
        bg = self._generate_basis()
        pix_pos = np.zeros((3, 4, 8, 185, 388))
        for i in range(4):
            for j in range(8):
                grid = bg.grid_as_explicit(i*8 + j)
                
                if grid.shape == (388,185,3):
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
    
        
    def _process_parameters(self):
        """
        Alignment calibrations as defined for psana. Injects
        
        -- self.tilt
        -- self.sec_offset
        -- self.section_centers
        
        into the local namespace.
        
        This method adapted from XtcExplorer's cspad.py
        """
        
        # angle of each section = rotation (nx90 degrees) + tilt (small angles)
        # ... (4 rows (quads) x 8 columns (sections))
        self.section_angles = self.rotation + self.tilt
        
        # position of sections in each quadrant (quadrant coordinates)
        # ... (3*4 rows (4 quads, 3 xyz coordinates) x 8 columns (sections))
        self.section_centers = np.reshape( self.center + self.center_corr, (3,4,8) )
        
        # quadrant offset parameters (w.r.t. image 0,0 in upper left coner)
        # ... (3 rows (xyz) x 4 columns (quads) )
        quad_position = self.offset + self.offset_corr

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
        
        # inject offset_corr_xy
        self.offset_corr_xy = self.offset_corr[:2,:]
        
        return     
    
    
    def _rotate_xy(self, vector, degrees_cw):
        """
        Perform a rotation in the x-y plane of a 3-vector
        """
        ct = np.cos( np.deg2rad(degrees_cw) )
        st = np.sin( np.deg2rad(degrees_cw) )
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
        # but the geometry specification here outputs these quads rotated 90-deg
        # counter-clockwise:
        #
        #     Q1     Q2
        # 
        #         x
        # 
        #     Q0     Q3
        # 
        # if you use range(4) to loop over the quads I am not 100% sure right 
        # now why this is the case, but the modified loop appears to give the
        # correct geometry. This should be triple-checked.
        # -- TJL 28.2.13

        bg = BasisGrid()
        
        # assemble each 2x1 on the quad
        for quad_index in [1,2,3,0]:
            for i in range(8):
                
                # in the below, the following convention is correct:
                # slow == row == x
                # fast == col == y
                
                shape = (185, 388) # slow, fast

                # we have to re-orient each 2x1 -- the geometrical positions
                # must also match the way intensity data are layed out in
                # memory...
                
                # re-orient quads 0,1 & 4,5, which are rotated
                if (i==0 or i==1):
                    # reverse slow dim, switch slow/fast
                    s = np.array([  0.0, 0.10992, 0.0 ])
                    f = np.array([ -0.10992, 0.0, 0.0 ])
                elif (i==4 or i==5):
                    # reverse fast dim, switch slow/fast
                    s = np.array([ 0.0, -0.10992, 0.0 ])
                    f = np.array([ 0.10992,  0.0, 0.0 ])
                else:
                    s = np.array([ 0.10992, 0.0, 0.0 ])
                    f = np.array([ 0.0, 0.10992, 0.0 ])

                # now, apply `tilt` correction - a small rotation in x-y
                s = self._rotate_xy(s, self.tilt[quad_index][i])
                f = self._rotate_xy(f, self.tilt[quad_index][i])

                # find the center of the 2x1
                cx = self.sec_offset[0] + self.section_centers[0][quad_index][i]
                cy = self.sec_offset[1] + self.section_centers[1][quad_index][i]
                cz = self.sec_offset[2] + self.section_centers[2][quad_index][i]
                                
                center =  np.array([850-cx, 850-cy, cz])
                center *= 0.10992 # convert to mm
                
                # convert to p -- the natural BasisGrid convention. This is
                # necessary to get the rotation below correct.
                x = (np.array(shape) - 1)
                center_correction = ((x[0] * s) + (x[1] * f)) / 2.
                p = center.copy()
                p -= center_correction
                
                # rotate each quad by the appropriate amount ( n * 90-deg CW 
                # rotations for quads n = { 0, 1, 2, 3 } ), then translate
                # them so that the top-left corners of an 850 x 850 pixel box
                # overlap. This remains true to psana convention.
                                
                # perform the rotation
                s = self._rotate_xy( s, 90*(4-quad_index) ) #+ self.quad_rotation[quad_index])
                f = self._rotate_xy( f, 90*(4-quad_index) ) #+ self.quad_rotation[quad_index])
                p = self._rotate_xy( p, 90*(4-quad_index) ) #+ self.quad_rotation[quad_index])
                
                # now translate so that the top-left corners of an 850 x 850 box
                # overlap
                if quad_index == 0:
                    translation = np.array([0.0,     0.0, 0.0])
                elif quad_index == 1:
                    translation = np.array([  0.0, 850.0, 0.0])
                elif quad_index == 2:
                    translation = np.array([850.0, 850.0, 0.0])    
                elif quad_index == 3:
                    translation = np.array([850.0,   0.0, 0.0])    
                
                p += translation * 0.10992 # convert to microns
                
                # add the quad offset, which defines the relative spatial
                # orientations of each quad
                p += self.quad_offset[:,quad_index] * 0.10992
                
                # finally, add these to our basis grid
                bg.add_grid(p, s, f, shape)
                
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

               pair = interp.rotate(pair, self.tilt[quad_index][i],
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

            if quad_index>0:
                # reorient the quad_index_image as needed
                quad_index_image = np.rot90( quad_index_image, 4-quad_index )

            qoff_x = int( self.quad_offset[0,quad_index] )
            qoff_y = int( self.quad_offset[1,quad_index] )
            
            if (qoff_x < 0) or (qoff_x >= bounds):
                raise Exception('qoff_x: %d out of bounds [0,%d)' % (qoff_x, bounds))
            if (qoff_y < 0) or (qoff_y >= bounds):
                raise Exception('qoff_y: %d out of bounds [0,%d)' % (qoff_y, bounds))
            
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
            
        Returns
        -------
        self : cspad.CSPad
            A CSPad object.
        """
        
        defaults = { 'beam_location'  : np.array([900.0, 870.0]),
                     'beam_vector'    : np.array([0.0, 0.0, 100.0]),
                     'offset_corr_xy' : np.zeros((2, 4))       # never accessed
                   }
        
        param_dict = {}
    
        # look for and load the parameters that we "guarentee" are there
        for p in _array_sizes.keys():
            
            # scan through the possible files and find the right one matching
            # our run number
            files = glob( pjoin(path, p, '*') )
            
            # if theres nothing in the dir, complain, unless we're looking for
            if len(files) == 0:
                if p in defaults.keys():
                    print "Could not locate dir '%s', using default value" % p
                    param_dict[p] = defaults[p]
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
        
        
    @classmethod    
    def default(cls):
        """
        Generate a default CSPad geometry.
        
        Returns
        -------
        self : cspad.CSPad
            A CSPad object.
        """
        defaults = default.DefaultCSPadParams()
        return cls(defaults.__dict__)
    

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
        print "Generating metrology from: %s" % metrology_file
        self.load(metrology_file)
        self.quad_offset = np.zeros((3,4))
        self.pixel_size = 0.10992
        
    
    @property
    def _flat_pixel_positions(self):
        """
        Returns all the pixels as an (x,y,z) / N x 3 array.
        """
        return np.vstack(( self.x_coordinates.flatten(),
                           self.y_coordinates.flatten(),
                           self.z_coordinates.flatten() )).T
                           

    @property           
    def pixel_positions(self):
        return np.array((self.x_coordinates,
                         self.y_coordinates,
                         self.z_coordinates))
                         
    
    # the following methods are an attempt to provide an interface similar
    # to that of CSPad -- the only parameter we really have, though, is
    # "offset", which is the quad offset
    
    def get_param(self, param_name):
        """
        Parameter getter function.
        """
        if param_name == 'offset_corr':
            return self.quad_offset
        elif param_name == 'offset_corr_xy':
            xy = self.quad_offset[:2,:]
            assert xy.shape == (2,4)
            return xy
        else:
            raise ValueError('Metrology objects only has `offset_corr` and '
                             '`offset_corr_xy` parameters, got: %s' % param_name)

        
    def set_param(self, param_name, value):
        """
        Parameter setter.
        """
        if param_name == 'offset_corr':
            if value.shape == (3,4):
                self.quad_offset = value
            else:
                raise ValueError('`offset_corr` should have shape (3,4),'
                                 ' got: %s' % str(value.shape))
        elif param_name == 'offset_corr_xy':
            if value.shape == (2,4):
                self.quad_offset[:2,:] = value
            else:
                raise ValueError('`offset_corr_xy` should have shape (2,4),'
                                 ' got: %s' % str(value.shape))
        else:
            raise ValueError('Metrology objects only has `offset_corr` and '
                             '`offset_corr_xy` parameters, got: %s' % param_name)


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
        
        # ----------------------------------------------------------------------
        # this is supremely stupid, but to ensure that (at least for now) I am
        # consistent w/Mikhail, I am going to re-parse the metrology file and
        # apply exactly his scheme. Why there is no unified codebase in pyana,
        # I don't know...

        self.points_for_quadrants = [ [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26] ]

        #Base index for 2x1:
        #             0  1   2  3   4   5   6   7
        self.ibase = [5, 1, 13, 9, 17, 21, 29, 25]

        self.metrology_array = np.zeros( (4,33,4), dtype=np.int32 )

        f = open(metrology_file, 'r')
        # Print out 7th entry in each line.
        for line in f:

            # ignore empty lines
            if len(line) == 1:
                continue 

            list_of_fields = line.split()

            # ignore quad header lines
            if list_of_fields[0] == 'Quad': 
                quad = int(list_of_fields[1])
                continue

            # ignore the title lines
            if list_of_fields[0] == 'Sensor':
                continue
            
            # Ignore lines with non-expected number of fields
            if len(list_of_fields) != 4:
                if verbose:
                    print 'WARNING: len(list_of_fields) =', len(list_of_fields)
                    print 'RECORD IS IGNORED due to unexpected format of the line:',line
                continue              

            point = int(list_of_fields[0])
            X = int(list_of_fields[1])
            Y = int(list_of_fields[2])
            Z = int(list_of_fields[3])

            self.metrology_array[quad,point,0] = point
            self.metrology_array[quad,point,1] = X
            self.metrology_array[quad,point,2] = Y
            self.metrology_array[quad,point,3] = Z

        f.close()
        
        return
        
        
    def apply_quad_offset(self, quad_offset):
        """
        Apply specific geometry corrections, defined by the psana parameters,
        to the metrology pixel map.
        
        Parameters
        ----------
        quad_offset : np.ndarray
            A shape-(3,4) array describing the (x,y) offset to apply to each
            quad.
        
        Returns
        -------
        x, y, z : np.ndarrays of shape (4,8,185,388)
            The x,y,z locations of each pixel in real-space, corrected.
        """
        
        if not quad_offset.shape == (3,4):
            raise ValueError('`quad_offset` is wrong shape, must be (3,4), got'
                             ' %s' % str(quad_offset.shape))
        
        # keep track of the total offsets applied so far
        self.quad_offset += quad_offset
        
        for i in range(4):
            self.x_coordinates[i] += quad_offset[0,i]
            self.y_coordinates[i] += quad_offset[1,i]
            self.z_coordinates[i] += quad_offset[2,i]
        
        return
        
        
    def to_basisgrid(self):
        """
        Convert the optical metrology to a basis grid representation.
        
        Returns
        -------
        bg : cspad.BasisGrid
            A basis representation of the metrology.
        """
        
        bg = BasisGrid()
        
        for i in range(4):
            for j in range(8):
                p = np.array([ self.x_coordinates[i,j,0,0],
                               self.y_coordinates[i,j,0,0],
                               self.z_coordinates[i,j,0,0] ])
                s = np.array([ self.x_coordinates[i,j,1,0] - self.x_coordinates[i,j,0,0],
                               self.y_coordinates[i,j,1,0] - self.x_coordinates[i,j,0,0],
                               self.z_coordinates[i,j,1,0] - self.x_coordinates[i,j,0,0] ])
                f = np.array([ self.x_coordinates[i,j,0,1] - self.x_coordinates[i,j,0,0],
                               self.y_coordinates[i,j,0,1] - self.x_coordinates[i,j,0,0],
                               self.z_coordinates[i,j,0,1] - self.x_coordinates[i,j,0,0] ])
                shape = (185, 388)
                bg.add_grid(p, s, f, shape)
        
        return bg
    
        
    def _compute_center_coordinates(self):
        """
        Compute the center positions of each 2x1
        """

        # this is *not* pixel units (why on earth is dtype int)?!?
        arrXmu = np.zeros( (4,8), dtype=np.int32 )
        arrYmu = np.zeros( (4,8), dtype=np.int32 )
        arrZmu = np.zeros( (4,8), dtype=np.int32 )
        
        # store the pixel unit computation
        self.centers = np.zeros( (3,4,8), dtype=np.float )

        ix = 1
        iy = 2
        iz = 3

        for quad in range(4) :
            for pair in range(8) :

                icor1 = self.ibase[pair]   
                icor2 = self.ibase[pair] + 1
                icor3 = self.ibase[pair] + 2
                icor4 = self.ibase[pair] + 3

                X = 0.25 * ( self.metrology_array[quad,icor1,ix]
                           + self.metrology_array[quad,icor2,ix]
                           + self.metrology_array[quad,icor3,ix]
                           + self.metrology_array[quad,icor4,ix] )

                Y = 0.25 * ( self.metrology_array[quad,icor1,iy]
                           + self.metrology_array[quad,icor2,iy]
                           + self.metrology_array[quad,icor3,iy]
                           + self.metrology_array[quad,icor4,iy] ) 

                Z = 0.25 * ( self.metrology_array[quad,icor1,iz]
                           + self.metrology_array[quad,icor2,iz]
                           + self.metrology_array[quad,icor3,iz]
                           + self.metrology_array[quad,icor4,iz] ) 

                arrXmu[quad][pair] = X
                arrYmu[quad][pair] = Y
                arrZmu[quad][pair] = Z

                xyz = np.array((X, Y, Z)).astype(np.float)
                self.centers[:,quad,pair] = xyz / self.pixel_size
                
        self.centers = self.centers.reshape(12,8)
                
        return
    

    def _compute_length_width_angle(self):
        """
        Really all this is used for right now is to compute the `tilt` parameter
        """
        
        # MORE supreme stupidity -- this again duplicates a lot of the
        # computation already done, but I'm copying it to generate a test case

        S1  = np.zeros( (4,8), dtype=np.int32 )
        S2  = np.zeros( (4,8), dtype=np.int32 )

        dS1 = np.zeros( (4,8), dtype=np.int32 )
        dS2 = np.zeros( (4,8), dtype=np.int32 )

        L1  = np.zeros( (4,8), dtype=np.int32 )
        L2  = np.zeros( (4,8), dtype=np.int32 )

        dL1 = np.zeros( (4,8), dtype=np.int32 )
        dL2 = np.zeros( (4,8), dtype=np.int32 )

        D1  = np.zeros( (4,8), dtype=np.int32 )
        D2  = np.zeros( (4,8), dtype=np.int32 )
        dD  = np.zeros( (4,8), dtype=np.int32 )

        ddS = np.zeros( (4,8), dtype=np.int32 )
        ddL = np.zeros( (4,8), dtype=np.int32 )

        self.tilt = np.zeros( (4,8), dtype=np.float32 )

        # no idea on the below...
        ix = 1
        iy = 2

        for quad in range(4):
            for pair in range(8):

                icor1 = self.ibase[pair]   
                icor2 = self.ibase[pair] + 1
                icor3 = self.ibase[pair] + 2
                icor4 = self.ibase[pair] + 3

                if pair == 0 or pair == 1 or pair == 4 or pair == 5:

                    S1[quad][pair] = self.metrology_array[quad,icor2,iy] - self.metrology_array[quad,icor1,iy]
                    S2[quad][pair] = self.metrology_array[quad,icor3,iy] - self.metrology_array[quad,icor4,iy]

                    dS1[quad][pair] = self.metrology_array[quad,icor4,iy] - self.metrology_array[quad,icor1,iy]
                    dS2[quad][pair] = self.metrology_array[quad,icor3,iy] - self.metrology_array[quad,icor2,iy]

                    L1[quad][pair] = self.metrology_array[quad,icor4,ix] - self.metrology_array[quad,icor1,ix]
                    L2[quad][pair] = self.metrology_array[quad,icor3,ix] - self.metrology_array[quad,icor2,ix]

                    dL1[quad][pair] = self.metrology_array[quad,icor2,ix] - self.metrology_array[quad,icor1,ix]
                    dL2[quad][pair] = self.metrology_array[quad,icor3,ix] - self.metrology_array[quad,icor4,ix]


                else:

                    S1[quad][pair] =   self.metrology_array[quad,icor4,ix] - self.metrology_array[quad,icor1,ix]
                    S2[quad][pair] =   self.metrology_array[quad,icor3,ix] - self.metrology_array[quad,icor2,ix]
                                                                                       
                    dS1[quad][pair] = -(self.metrology_array[quad,icor2,ix] - self.metrology_array[quad,icor1,ix]) # sign is chosen 
                    dS2[quad][pair] = -(self.metrology_array[quad,icor3,ix] - self.metrology_array[quad,icor4,ix]) # for positive phi

                    L1[quad][pair] =   self.metrology_array[quad,icor2,iy] - self.metrology_array[quad,icor1,iy]
                    L2[quad][pair] =   self.metrology_array[quad,icor3,iy] - self.metrology_array[quad,icor4,iy]
                                                                                       
                    dL1[quad][pair] =   self.metrology_array[quad,icor4,iy] - self.metrology_array[quad,icor1,iy]
                    dL2[quad][pair] =   self.metrology_array[quad,icor3,iy] - self.metrology_array[quad,icor2,iy]


                diag1x = float(self.metrology_array[quad,icor1,ix] - self.metrology_array[quad,icor3,ix])
                diag2x = float(self.metrology_array[quad,icor2,ix] - self.metrology_array[quad,icor4,ix])
                diag1y = float(self.metrology_array[quad,icor1,iy] - self.metrology_array[quad,icor3,iy])
                diag2y = float(self.metrology_array[quad,icor2,iy] - self.metrology_array[quad,icor4,iy])

                D1[quad][pair] = int( np.sqrt(diag1x*diag1x + diag1y*diag1y) )
                D2[quad][pair] = int( np.sqrt(diag2x*diag2x + diag2y*diag2y) )
                dD[quad][pair] = D1[quad][pair] - D2[quad][pair]

                ddS[quad][pair] = dS1[quad][pair] - dS2[quad][pair]
                ddL[quad][pair] = dL1[quad][pair] - dL2[quad][pair]

                ang1 = 0
                ang2 = 0
                if dL1[quad][pair] != 0:
                    ang1 = float(dS1[quad][pair]) / L1[quad][pair]
                if dL2[quad][pair] != 0:
                    ang2 = float(dS2[quad][pair]) / L2[quad][pair]

                angav = (ang1 + ang2) * 0.5
                ang_deg = 180.0 / 3.1415927 * angav

                self.tilt[quad][pair] = ang_deg
        
        return
    
        
    def to_cspad(self):
        """
        Convert the optical metrology to a set of psana alignment parameters,
        in the form of a cspad object.
        
        Returns
        -------
        cs : cspad.CSPad
            A CSPad object with the tilt/center parameters changed to reflect
            the metrology. The rest of the parameters are set to the pyana
            default values.
        """
        
        # compute the necessary quantities for the parameters: centers, tilt
        self._compute_center_coordinates()
        self._compute_length_width_angle()
        
        # get a default CSPad object and set centers, tilt
        cs = CSPad.default()
        cs.set_param('center', self.centers)
        cs.set_param('tilt', self.tilt)
        
        return cs
    
        
    def to_dir(self, dirname):
        """
        A quick and dirty save method. Probably should write h5s.
        
        todo
        """
        
        if not os.path.isdir(dirname):
            print "Creating directory: %s" % dirname
            os.mkdir(dirname)
            
        np.savez('x_coordinates.npz', self.x_coordinates)
        np.savez('y_coordinates.npz', self.y_coordinates)
        np.savez('z_coordinates.npz', self.z_coordinates)
        
        np.savetxt('offset_corr.dat', self.offset_corr)
        
        print "Wrote metrology specification to: %s" % dirname
        
        return
        
    @classmethod
    def from_dir(self, dirname):
        raise NotImplementedError()
        

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
