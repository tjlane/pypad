

"""
export.py

Library for exporting autogeom geometries to other software packages.
"""

import numpy as np
import h5py

    
def to_cheetah(geometry, filename="pixelmap-cheetah-raw.h5"):
    """
    Convert a CSPad or Metrology object to a Cheetah h5 pixel map.
    
    Writes: 'pixelmap-cheetah-raw.h5'
    
    Parameters
    ----------
    geometry : cspad.CSPad OR cspad.Metrology
        The detector geometry to write to disk
	filename : string
		The file name for the output pixel map
    """
    
    print "Exporting to Cheetah format..."
    
    coordinates = ['x', 'y', 'z']
    
    pp = geometry.pixel_positions
    assert pp.shape == (3, 4, 16, 185, 194)
    
    # write an h5
    f = h5py.File(filename, 'w')
    
    # iterate over x/y/z
    for x in range(len(coordinates)):
        
        cheetah_image = np.zeros((1480, 1552), dtype=np.float32)
        
        # iterate over each 2x1/quad (note switch)
        for i in range(8):
            for j in range(4):
                x_start = 185 * i
                x_stop  = 185 * (i+1)
                y_start = 388 * j
                y_stop  = 388 * (j+1)
                quad = (j + 2) % 4 # cheetah's quads are not quite the same...
                # jas: ensure cheetah geometry is in [m] and not [mm]
                two_by_one = np.hstack(( pp[x,quad,i,:,:], pp[x,quad,i*2+1,:,:] )) / 1000.0
                cheetah_image[x_start:x_stop,y_start:y_stop] = two_by_one
        
        f['/%s' % coordinates[x]] = cheetah_image
    
    print "Wrote: %s" % (filename)
    f.close()
    
    return

    
def to_odin(geometry, filename):
    """
    Generate an ODIN detector object and write it to disk.
    
    Parameters
    ----------
    geometry : cspad.CSPad OR cspad.Metrology
        The detector geometry to write to disk
        
    filname : str
        The name of file to write. Will end in '.dtc'
    """
    
    try:
        from odin import xray
    except ImportError as e:
        raise ImportError('Cannot find ODIN. You must have ODIN installed to '
                          'export to ODIN. Download and install ODIN from '
                          'https://github.com/tjlane/odin')
        
    raise NotImplementedError('Waiting for hooks on ODIN side.')
    
    return
