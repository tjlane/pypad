
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
export.py

Library for exporting autogeom geometries to other software packages.
"""

import numpy as np
import h5py
import math
    
def to_cheetah(geometry, filename="pixelmap-cheetah-raw.h5"):
    """
    Convert a CSPad or Metrology object to a Cheetah h5 pixel map.
    
    Writes: 'pixelmap-cheetah-raw.h5'
    
    Parameters
    ----------
    geometry : cspad.CSPad
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
    for xyz in range(len(coordinates)):
        
        cheetah_image = np.zeros((1480, 1552), dtype=np.float32)
        
        # iterate over each 2x1/quad (note switch)
        
        # new
        for q in range(4):
            for twoXone in range(8):
                
                x_start = 388 * q
                x_stop  = 388 * (q+1)
                
                y_start = 185 * twoXone
                y_stop  = 185 * (twoXone + 1)
                
                # convert mm --> m
                cheetah_image[y_start:y_stop,x_start:x_stop] = np.hstack(( pp[xyz,q,twoXone*2,:,:],
                                                                           pp[xyz,q,twoXone*2+1,:,:] )) / 1000.0
        
        
        f['/%s' % coordinates[xyz]] = cheetah_image
    
    print "Wrote: %s" % (filename)
    f.close()
    
    return

    
def to_thor(geometry, energy, distance_offset, filename=None):
    """
    Generate an ODIN detector object and write it to disk.
    
    Parameters
    ----------
    geometry : cspad.CSPad
        The detector geometry to write to disk
        
    energy : float
        The energy of the beam, in eV.
        
    distance_offset : float
        A distance offset (along the z-direction), in mm. This can be measured
        directly using a pypad.autogeom.score.PowderReference instance.
        
    filname : str
        The name of file to write. Will end in '.dtc'. If `None`, don't
        save the file.
    """
    
    try:
        from thor import xray
    except ImportError as e:
        print e
        raise ImportError('Cannot find Thor. You must have Thor installed to '
                          'export to Thor. Download and install Thor from '
                          'https://github.com/tjlane/thor')
        
    # we have to generate an thor-type basis grid -- maintains code separation
    pypad_bg = geometry.basis_repr
    thor_bg  = xray.BasisGrid()
    for i in range(pypad_bg.num_grids):
        p, s, f, shp = pypad_bg.get_grid(i)
        p[2] += distance_offset # translate to make the origin the interaction site
        thor_bg.add_grid(p, s, f, shp)
    
        
    # recall that pypad assumes the beam is along the z-direction (CXI conven.)
    energy /= 1000.0 # convert to keV : pypad is eV, Thor is keV
    b = xray.Beam(1e11, energy=energy) # 1e11 photons per shot
    d = xray.Detector(thor_bg, b)
    if filename:
        d.save(filename)
    
    return d
    
    
def to_crystfel(geometry, filename, intensity_file_type='cheetah'):
    """
    Write a CSPAD geometry to disk in CrystFEL format. Note that some fields
    will be written but left blank -- these are fields you probably should
    fill in before performing any computations in CrystFEL, but are information
    that pypad has no handle on (e.g. detector gain).
    
    Thanks to Rick Kirian & Tom White for assistance with this function.
    
    Parameters
    ----------
    geometry : cspad.CSPad
        The detector geometry to write to disk
        
    filname : str
        The name of file to write. Will end in '.dtc'
        
    Optional Parameters
    -------------------
    intensity_file_type : str, {'cheetah'}
        The kind of file this geometry file will be used with. Necessary to tell
        CrystFEL how intensity data map onto the detector
    """
    
    pixel_size = 0.10992 # in mm
    bg = geometry.basis_repr
    
    def get_sign(v):
        if v >= 0:
            s = '+'
        else:
            s = '-'
        return s
    
    
    if intensity_file_type == 'cheetah':
        
        # this is complex, so I went the lazy route and copied an
        # existing file
        intensity_map = crystfel_cheetah_intensities.split('-')
        assert len(intensity_map) == 64
        
    else:
        raise ValueError('Cannot write geometries for '
                        '`intensity_file_type`: %s, only currently '
                        'have implemented writers for '
                        '{"cheetah"}' % intensity_file_type)
    
    
    
    with open(filename, 'w') as of:
    
        print >> of, "; This file contains a CSPAD geometry generated by PyPad"
        print >> of, "; http://www.github.com/tjlane/pypad"
    
        header = """
; --- VALUES YOU MAY WANT TO FILL IN MANUALLY ---
; pypad cannot fill in these values for you :(

;clen = 
;coffset = 

;adu_per_eV = 
;max_adu = 

;mask = 
;mask_good = 
;mask_bad = 

;bad_jet/min_x = 
;bad_jet/max_x = 
;bad_jet/min_y = 
;bad_jet/max_y = 
; -----------------------------------------------


"""
        
        print >> of, header
        
        # iterate over each basis grid object
        for quad in range(4):
            for asic in range(16):
                
                grid_index = quad * 16 + asic
                p, s, f, sp = bg.get_grid(grid_index)
    
                panel_name = "q%da%d" % (quad, asic)
                
                # tell crystFEL how read intensity values in a file
                print >> of, intensity_map[grid_index].strip()
                
                # these are constant for the CSPAD
                print >> of, "%s/badrow_direction = -" % panel_name
                print >> of, "%s/res = 9097.525473" % panel_name
                
                # write the basis vectors           
                sqt = math.sqrt(f[0]**2 + f[1]**2) 
                print >> of, "%s/fs = %s%fx %s%fy" % ( panel_name,
                                                       get_sign(f[0]/sqt), abs(f[0]/sqt), 
                                                       get_sign(f[1]/sqt), abs(f[1]/sqt) )
                sqt = math.sqrt(s[0]**2 + s[1]**2)
                print >> of, "%s/ss = %s%fx %s%fy" % ( panel_name,
                                                       get_sign(s[0]/sqt), abs(s[0]/sqt), 
                                                       get_sign(s[1]/sqt), abs(s[1]/sqt) )
                
                # write the corner positions
                tagcx = "%s/corner_x" % panel_name
                tagcy = "%s/corner_y" % panel_name
            
                # CrystFEL measures the corner from the actual *corner*, and not
                # the center of the corner pixel, hence the 0.5's
                print >> of, "%s = %f" % (tagcx, float(p[0])/pixel_size - 0.5 )
                print >> of, "%s = %f" % (tagcy, float(p[1])/pixel_size - 0.5 )
                
                # this tells CrystFEL to use this panel
                print >> of, "%s/no_index = 0" % panel_name
                
                print >> of, "" # new line
    
    return
    
    
def to_text(geometry, filename):
    """
    Write a CSPAD geometry to disk in the following format:
    
        p_x p_y p_z     s_x s_y s_z     f_x f_y f_z
        ...
        
    and include some comments.    
    
    
    Parameters
    ----------
    geometry : cspad.CSPad
        The detector geometry to write to disk
        
    filname : str
        The name of file to write. Will end in '.dtc'
    """
    
    # generate a preamble
    preamble = """
# This file contains a CSPAD geometry generated by PyPad: 
# http://www.github.com/tjlane/pypad
#
# The following is a basis grid representation with the following vectors
#
#   p : position vector for an ASIC
#   s : slow-scan pixel vector
#   f : fast-scan pixel vector
#
# all units are mm. Each ASIC is 185 x 194 pixels.
# See the PyPad documentation for more information.


#            p_x        p_y        p_z           s_x        s_y        s_z           f_x        f_y        f_z
"""

    # loop over each grid element and add it to the file
    bg = geometry.basis_repr
    body = ""
    
    
    def format(s, total_len=10):
        """
        A little formatting function
        """
        sf = '%.5f' % s
        pad = total_len - len(sf)
        if pad > 0:
            sf = ' ' * pad + sf
        return sf
        
    
    for i in range(bg.num_grids):
        
        # if we're starting a new quad, note that in the file
        if i % 16 == 0:
            body += ('\n# QUAD %d\n' % (i/16))
        
        # add the basis grid
        p, s, f, shp = bg.get_grid(i)
        strp = ' '.join( [ format(x) for x in p ] )
        strs = ' '.join( [ format(x) for x in s ] )
        strf = ' '.join( [ format(x) for x in f ] )
        
        tb = ' ' * 4
        asic = str(i)
        if len(asic) == 1:
            asic = ' ' + asic
        
        body += (asic + tb + strp + tb + strs + tb + strf + '\n')
        
    f = open(filename, 'w')
    f.write(preamble + body)
    f.close()
    
    print "Wrote CSPAD to text at: %s" % filename
    
    return



#  ---------- REFERENCE DATA ---------------------------------------------------


crystfel_cheetah_intensities = """q0a0/min_fs = 0
q0a0/min_ss = 0
q0a0/max_fs = 193
q0a0/max_ss = 184
-
q0a1/min_fs = 194
q0a1/min_ss = 0
q0a1/max_fs = 387
q0a1/max_ss = 184
-
q0a2/min_fs = 0
q0a2/min_ss = 185
q0a2/max_fs = 193
q0a2/max_ss = 369
-
q0a3/min_fs = 194
q0a3/min_ss = 185
q0a3/max_fs = 387
q0a3/max_ss = 369
-
q0a4/min_fs = 0
q0a4/min_ss = 370
q0a4/max_fs = 193
q0a4/max_ss = 554
-
q0a5/min_fs = 194
q0a5/min_ss = 370
q0a5/max_fs = 387
q0a5/max_ss = 554
-
q0a6/min_fs = 0
q0a6/min_ss = 555
q0a6/max_fs = 193
q0a6/max_ss = 739
-
q0a7/min_fs = 194
q0a7/min_ss = 555
q0a7/max_fs = 387
q0a7/max_ss = 739
-
q0a8/min_fs = 0
q0a8/min_ss = 740
q0a8/max_fs = 193
q0a8/max_ss = 924
-
q0a9/min_fs = 194
q0a9/min_ss = 740
q0a9/max_fs = 387
q0a9/max_ss = 924
-
q0a10/min_fs = 0
q0a10/min_ss = 925
q0a10/max_fs = 193
q0a10/max_ss = 1109
-
q0a11/min_fs = 194
q0a11/min_ss = 925
q0a11/max_fs = 387
q0a11/max_ss = 1109
-
q0a12/min_fs = 0
q0a12/min_ss = 1110
q0a12/max_fs = 193
q0a12/max_ss = 1294
-
q0a13/min_fs = 194
q0a13/min_ss = 1110
q0a13/max_fs = 387
q0a13/max_ss = 1294
-
q0a14/min_fs = 0
q0a14/min_ss = 1295
q0a14/max_fs = 193
q0a14/max_ss = 1479
-
q0a15/min_fs = 194
q0a15/min_ss = 1295
q0a15/max_fs = 387
q0a15/max_ss = 1479
-
q1a0/min_fs = 388
q1a0/min_ss = 0
q1a0/max_fs = 581
q1a0/max_ss = 184
-
q1a1/min_fs = 582
q1a1/min_ss = 0
q1a1/max_fs = 775
q1a1/max_ss = 184
-
q1a2/min_fs = 388
q1a2/min_ss = 185
q1a2/max_fs = 581
q1a2/max_ss = 369
-
q1a3/min_fs = 582
q1a3/min_ss = 185
q1a3/max_fs = 775
q1a3/max_ss = 369
-
q1a4/min_fs = 388
q1a4/min_ss = 370
q1a4/max_fs = 581
q1a4/max_ss = 554
-
q1a5/min_fs = 582
q1a5/min_ss = 370
q1a5/max_fs = 775
q1a5/max_ss = 554
-
q1a6/min_fs = 388
q1a6/min_ss = 555
q1a6/max_fs = 581
q1a6/max_ss = 739
-
q1a7/min_fs = 582
q1a7/min_ss = 555
q1a7/max_fs = 775
q1a7/max_ss = 739
-
q1a8/min_fs = 388
q1a8/min_ss = 740
q1a8/max_fs = 581
q1a8/max_ss = 924
-
q1a9/min_fs = 582
q1a9/min_ss = 740
q1a9/max_fs = 775
q1a9/max_ss = 924
-
q1a10/min_fs = 388
q1a10/min_ss = 925
q1a10/max_fs = 581
q1a10/max_ss = 1109
-
q1a11/min_fs = 582
q1a11/min_ss = 925
q1a11/max_fs = 775
q1a11/max_ss = 1109
-
q1a12/min_fs = 388
q1a12/min_ss = 1110
q1a12/max_fs = 581
q1a12/max_ss = 1294
-
q1a13/min_fs = 582
q1a13/min_ss = 1110
q1a13/max_fs = 775
q1a13/max_ss = 1294
-
q1a14/min_fs = 388
q1a14/min_ss = 1295
q1a14/max_fs = 581
q1a14/max_ss = 1479
-
q1a15/min_fs = 582
q1a15/min_ss = 1295
q1a15/max_fs = 775
q1a15/max_ss = 1479
-
q2a0/min_fs = 776
q2a0/min_ss = 0
q2a0/max_fs = 969
q2a0/max_ss = 184
-
q2a1/min_fs = 970
q2a1/min_ss = 0
q2a1/max_fs = 1163
q2a1/max_ss = 184
-
q2a2/min_fs = 776
q2a2/min_ss = 185
q2a2/max_fs = 969
q2a2/max_ss = 369
-
q2a3/min_fs = 970
q2a3/min_ss = 185
q2a3/max_fs = 1163
q2a3/max_ss = 369
-
q2a4/min_fs = 776
q2a4/min_ss = 370
q2a4/max_fs = 969
q2a4/max_ss = 554
-
q2a5/min_fs = 970
q2a5/min_ss = 370
q2a5/max_fs = 1163
q2a5/max_ss = 554
-
q2a6/min_fs = 776
q2a6/min_ss = 555
q2a6/max_fs = 969
q2a6/max_ss = 739
-
q2a7/min_fs = 970
q2a7/min_ss = 555
q2a7/max_fs = 1163
q2a7/max_ss = 739
-
q2a8/min_fs = 776
q2a8/min_ss = 740
q2a8/max_fs = 969
q2a8/max_ss = 924
-
q2a9/min_fs = 970
q2a9/min_ss = 740
q2a9/max_fs = 1163
q2a9/max_ss = 924
-
q2a10/min_fs = 776
q2a10/min_ss = 925
q2a10/max_fs = 969
q2a10/max_ss = 1109
-
q2a11/min_fs = 970
q2a11/min_ss = 925
q2a11/max_fs = 1163
q2a11/max_ss = 1109
-
q2a12/min_fs = 776
q2a12/min_ss = 1110
q2a12/max_fs = 969
q2a12/max_ss = 1294
-
q2a13/min_fs = 970
q2a13/min_ss = 1110
q2a13/max_fs = 1163
q2a13/max_ss = 1294
-
q2a14/min_fs = 776
q2a14/min_ss = 1295
q2a14/max_fs = 969
q2a14/max_ss = 1479
-
q2a15/min_fs = 970
q2a15/min_ss = 1295
q2a15/max_fs = 1163
q2a15/max_ss = 1479
-
q3a0/min_fs = 1164
q3a0/min_ss = 0
q3a0/max_fs = 1357
q3a0/max_ss = 184
-
q3a1/min_fs = 1358
q3a1/min_ss = 0
q3a1/max_fs = 1551
q3a1/max_ss = 184
-
q3a2/min_fs = 1164
q3a2/min_ss = 185
q3a2/max_fs = 1357
q3a2/max_ss = 369
-
q3a3/min_fs = 1358
q3a3/min_ss = 185
q3a3/max_fs = 1551
q3a3/max_ss = 369
-
q3a4/min_fs = 1164
q3a4/min_ss = 370
q3a4/max_fs = 1357
q3a4/max_ss = 554
-
q3a5/min_fs = 1358
q3a5/min_ss = 370
q3a5/max_fs = 1551
q3a5/max_ss = 554
-
q3a6/min_fs = 1164
q3a6/min_ss = 555
q3a6/max_fs = 1357
q3a6/max_ss = 739
-
q3a7/min_fs = 1358
q3a7/min_ss = 555
q3a7/max_fs = 1551
q3a7/max_ss = 739
-
q3a8/min_fs = 1164
q3a8/min_ss = 740
q3a8/max_fs = 1357
q3a8/max_ss = 924
-
q3a9/min_fs = 1358
q3a9/min_ss = 740
q3a9/max_fs = 1551
q3a9/max_ss = 924
-
q3a10/min_fs = 1164
q3a10/min_ss = 925
q3a10/max_fs = 1357
q3a10/max_ss = 1109
-
q3a11/min_fs = 1358
q3a11/min_ss = 925
q3a11/max_fs = 1551
q3a11/max_ss = 1109
-
q3a12/min_fs = 1164
q3a12/min_ss = 1110
q3a12/max_fs = 1357
q3a12/max_ss = 1294
-
q3a13/min_fs = 1358
q3a13/min_ss = 1110
q3a13/max_fs = 1551
q3a13/max_ss = 1294
-
q3a14/min_fs = 1164
q3a14/min_ss = 1295
q3a14/max_fs = 1357
q3a14/max_ss = 1479
-
q3a15/min_fs = 1358
q3a15/min_ss = 1295
q3a15/max_fs = 1551
q3a15/max_ss = 1479"""

