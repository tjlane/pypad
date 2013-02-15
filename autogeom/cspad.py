#!/usr/bin/env python

"""
A high-level interface to pyana for assembling images from psana geometry
parameters. Main function to pay attention to is `assemble_image`.
"""

import sys
import os
from glob import glob
from os.path import join as pjoin

import numpy as np
import scipy.ndimage.interpolation as interp
import matplotlib.pyplot as plt

# PyAna -- testing only
# import PyCSPadImage.CalibParsDefault as cald
# import PyCSPadImage.CalibPars as calp
# import PyCSPadImage.CalibParsEvaluated as cpe
# import PyCSPadImage.CSPadConfigPars as ccp
# import PyCSPadImage.CSPadImageProducer as cip
# import PyCSPadImage.GlobalMethods as gm # getCSPadArrayFromFile for pedestal subtraction


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
                 'tilt' :           (4, 8) }


class CSPad(object):
    """
    This is a container class for saving, loading, and interacting with the 
    parameter set that defines the CSPad geometry.
    
    This class is largely based on XtcExplorer's cspad.py
    """
    
    def __init__(self, param_dict):
        """
        doc
        """
        
        # --- INTERNALS ------------------------------------------------------
        # There are some constants that this class needs to know about : what
        # parameters to expect, what their sizes should be, etc.

        self._param_list = _array_sizes.keys()

        # --------------------------------------------------------------------
        
        
        # we'll keep all the parameters in a dictionary called self.params -- 
        # we could expose them as attributes, but this allows us to be a
        # little more dynamic. Never know what to expect from PyAna...
        
        self._check_param_dict(param_dict)
        self.params = param_dict
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
        return
        
    
    @property
    def pixel_positions(self):
        """
        Compute and return the x,y,z positions of each pixel.
        """
        
        
        return
    
        
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
        self.rotation_array = self.params['rotation']
        self.tilt_array     = self.params['tilt']
        self.section_angles = self.rotation_array + self.tilt_array
        
        # position of sections in each quadrant (quadrant coordinates)
        # ... (3*4 rows (4 quads, 3 xyz coordinates) x 8 columns (sections))
        centers            = self.params['center']
        center_corrections = self.params['center_corr']
        self.section_centers = np.reshape( centers + center_corrections, (3,4,8) )
        
        # quadrant offset parameters (w.r.t. image 0,0 in upper left coner)
        # ... (3 rows (xyz) x 4 columns (quads) )
        quad_pos      = self.params['offset']
        quad_pos_corr = self.params['offset_corr']
        quad_position = quad_pos + quad_pos_corr

        # ... (3 rows (x,y,z) x 4 columns (section offset, quad offset, 
        # quad gap, quad shift)
        marg_gap_shift = self.params['marg_gap_shift']

        # break it down (extract each column, make full arrays to be added to 
        # the above ones)
        self.sec_offset = marg_gap_shift[:,0]

        quad_offset = marg_gap_shift[:,1]
        quad_gap = marg_gap_shift[:,2]
        quad_shift = marg_gap_shift[:,3]

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
    
    
    def _enforce_raw_img_shape(self, raw_image):
        """
        Make sure that the `raw_image` has shape: (4,8,185,388).
        
        Which is (quad, 2x1, fast, slow). This function will attempt to get
        the image into that form, and if it can't throw an error.
        """
        
        # XtcExporter format
        if raw_image.shape == (4,8,185,388):
            new_image = raw_image # we're good
            
        # PyCSPad format
        elif raw_image.shape == (32,185,388):
            new_image = np.zeros((4, 8, 185, 388))
            for i in range(8):
                for j in range(4):
                    psind = i + j * 8
                    new_image[j,i,:,:] = raw_image[psind,:,:]
            
        # Cheetah format
        elif raw_image.shape == (1480, 1552):
            new_image = np.zeros((4, 8, 185, 388))
            for i in range(8):
                for j in range(4):
                    x_start = 185 * i
                    x_stop  = 185 * (i+1)
                    y_start = 388 * j
                    y_stop  = 388 * (j+1)
                    psind = i + j * 8
                    new_image[j,i,:,:] = raw_image[x_start:x_stop,y_start:y_stop]
        
        else:
            raise ValueError("Cannot understand `raw_image`: does not have any"
                             " known dimension structure")
        
        return new_image
    
        
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
                   pair = interp.rotate(pair, self.tilt_array[quad_index][i])

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
        
        raw_image = self._enforce_raw_img_shape(raw_image)
        
        assembled_image = np.zeros((2*850+100, 2*850+100), dtype=np.float64)
        
        # iterate over quads
        for quad_index in range(4):

            quad_index_image = self._assemble_quad( raw_image[quad_index], quad_index )

            if quad_index>0:
                # reorient the quad_index_image as needed
                quad_index_image = np.rot90( quad_index_image, 4-quad_index)

            qoff_x = self.quad_offset[0,quad_index]
            qoff_y = self.quad_offset[1,quad_index]
            assembled_image[qoff_x:qoff_x+850, qoff_y:qoff_y+850] = quad_index_image

        return assembled_image
    
        
    def coordinate_map(self, metrology_file="cspad_2011-08-10-Metrology.txt",
                       verbose=False):
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

        origin = [[834,834],[834,834],[834,834],[834,834]]
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
        
        return self.x_coordinates, self.y_coordinates, self.z_coordinates
    

    def to_dir(self, dir_name):
        """
        Write the parameters to `dir_name`, in the standard psana/pyana format.
        
        This method will create a 2-level directory tree under `dir_name` to
        hold these parameters.
        
        Parameters
        ----------
        dir_name : str
            The name of the parent dir (parameter set) to place these parameters
            under.
        """
    
        check_param_dict(param_dict)
    
        if os.path.exists(dir_name):
            print "WARNING: rm'ing %s..." % dir_name
            os.system('rm -r %s' % dir_name)
    
        os.mkdir(dir_name)
    
        for key in param_dict.keys():
            os.mkdir( pjoin( dir_name, key ))
            fname = pjoin( dir_name, key, '0-end.data' )
            
            # the file-format for these two parameters is slightly different
            # from all the others... try to maintain that.
            if key == 'center' or key == 'center_corr':
            
                v = param_dict[key]
                v = v.reshape(3,4,8)
            
                f = open(fname, 'w')
                for i in range(3):
                    for j in range(4):
                        f.write( '\t'.join([ str(x) for x in v[i,j,:] ]) + '\n' )
                    f.write('\n')
            
            else:
                np.savetxt(fname, param_dict[key], fmt='%.2f')
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
        """
        
        '0-end.data'
    
        param_dict = {}
    
        to_load = ['center', 'center_corr', 'common_mode', 'filter',
                   'marg_gap_shift', 'offset', 'offset_corr', 'pedestals',
                   'pixel_status', 'quad_rotation' 'quad_tilt', 'rotation', 
                   'tilt']
    
        for p in _array_sizes.keys():
            
            # scan through the possible files and find the right one matching
            # our run number
            files = glob( pjoin(path, p, '*') )
            
            for f in files:
                
                start, end = os.path.basename(f).split('.')[0].split('-')
                start = int(start)
                if end == 'end':
                    end = 9999
                else:
                    end = int(end)
                    
                if (run_number >= start) and (run_number <= end):
                    filename = f
            
            print "Loading parameters in:", filename
            param_dict[p] = np.loadtxt(filename).reshape(_array_sizes[p])
    
        return cls(param_dict)
        
# ------------------------------------------------------------------------------
# TESTING
#
#
    

def get_avg_from_hdf5(hd5_image_file, calibration_path, dsname, start_event,
                      n_events_to_avg=5):
    """
    Extract an average image from a psana hdf5 file.    
    """

    print 'Getting raw CSPad event %d from file %s \ndataset %s' % (event, fname, dsname)
    ds1ev = hm.getAverageCSPadEvent( fname, dsname, start_event, nevents=n_events_to_avg )
    if not ds1ev.shape == (32, 185, 388):
        print 'WARNING: ds1ev.shape =', ds1ev.shape, "should be (32, 185, 388)"

    return ds1ev


def get_event_from_npz(npz_image_file):
    return np.load(npz_image_file)['arr_0']
    
    
def show_assembled_image(image):
    plt.imshow(image.T)
    plt.show()
    return

    
def pyana_assembly(raw_image, calibration_path, run_number=0):
    
    print 'Loading calibration parameters from: %s' % calibration_path
    calp.calibpars.setCalibParsForPath( run=run_number, path=calibration_path )
    calp.calibpars.printCalibPars()
    cpe.cpeval.printCalibParsEvaluatedAll()
    
    print 'Constructing the CSPad image from raw array'
    cspadimg = cip.CSPadImageProducer(rotation=0, tiltIsOn=True )
    image = cspadimg.getCSPadImage( raw_image )
    
    return image


def test_assembly_from_dir():
    
    # this one
    raw_image = get_event_from_npz('npz_examples/cxi64813_r58_evt1.npz')
    d = CSPad.from_dir('example_calibration_dir')
    show_assembled_image( d(raw_image) )
    
    # should be the same as this one
    ai = pyana_assembly(raw_image, 'example_calibration_dir')
    show_assembled_image(ai)
    
    return
    
    
def test_metrology():
    raw_image = get_event_from_npz('npz_examples/cxi64813_r58_evt1.npz')
    d = CSPad.from_dir('example_calibration_dir')
    x,y,z = d.coordinate_map(metrology_file="CSPad/cspad_2011-08-10-Metrology.txt")
    
    plt.imshow(x[0,0,:,:])
    plt.show()
    
    plt.imshow(x[0,1,:,:])
    plt.show()
    
    return
    

if __name__ == "__main__":
    #test_assembly_from_dir()
    #test_metrology()
