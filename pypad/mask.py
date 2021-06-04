
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
mask.py

Provides a "mask" object for CSPads.
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.widgets import Button
from matplotlib.path import Path

from pypad import utils
from pypad import read
from pypad.plot import ToggleButton


class PadMask(object):
    """
    An mask for a CSPad object.
    """
    
    def __init__(self):
        """
        Initialize a CSPad mask object.
        """
        self._masks = {}
        self._masks['base'] = self._blank_mask()
        return
    
    
    @property
    def mask(self):
        m = np.product( np.array(list(self._masks.values())), axis=0 )
        assert m.shape == (4,16,185,194)
        return m
    
        
    @property
    def mask2d(self):
        return utils.flatten_2x1s( self.mask )
    
    
    @property
    def inverted(self):
        """
        Invert the mask. Usually "True" is a good pixel, and "False" is a bad
        one, but this flips them.
        
        Returns
        -------
        """
        return np.logical_not(self.mask)
        
    
    @property
    def num_masked(self):
        """
        Returns the number of masked pixels
        """
        inv_mask = 1.0 - self.mask
        return np.sum(inv_mask)
        
        
    @property
    def types_applied(self):
        return list(self._masks.keys())
    
    
    def remove_mask(self, mask_name):
        """
        Remove a mask that has been applied.
        """
        
        if mask_name == 'base': # this one is special
            self._masks['base'] = self._blank_mask()
        elif not mask_name in list(self._masks.keys()):
            raise KeyError('Mask: %s not applied' % mask_name)
        else:
            x = self._masks.pop(mask_name)
            
        print("Removed mask: %s" % mask_name)
            
        return
    
        
    def _inject_mask(self, mask_name, mask, override_previous=False):
        """
        Add a new kind of mask to this mask object. Provides some typechecking.
        
        All this really does is deposit `mask` into dict self._masks with key
        `mask_name`.
        """
        
        assert type(mask_name) == str
        assert mask.shape == (4, 16, 185, 194)
        
        if not mask.dtype == np.bool:
            mask = mask.astype(np.bool)
            
        if (not mask_name in list(self._masks.keys())) or override_previous:
            self._masks[mask_name] = mask
        else:
            raise KeyError('Mask object already has `%s` mask.' % mask_name)
        
        return
        
        
    def _check_image(self, image):
        """
        Sanity check on `image`.
        """
        if not image.shape == (4, 16, 185, 194):
            raise ValueError('`image` must be shape (4, 16, 185, 194), got '
                             '%s' % str(image.shape))
        return
        
        
    def _blank_mask(self):
        """
        Utility function that just returns a blank mask.
        """
        return np.ones((4, 16, 185, 194), dtype=np.int32)
    
        
    # ----------
    # below we provide many methods of the form PadMask.mask_*(), which allow
    # one to mask pixels via various criteria (indicated by the *)
    
    # to add a new kind of mask, make a new method here. Follow mask_threshold
    # as a template
    
    def mask_pixel(self, quad, asic, x, y):
        """
        Mask a single pixel, or series of pixels. To do the latter, pass arrays
        as the arguments (even though the below says int).
        
        Parameters
        ----------
        quad : int
            [0,1,2,3], denoting the quad.
            
        asic : int
            Int in [0,7], denoting asic.
            
        x : int
            Int in [0,184], denoting x position.
        
        y : int
            Int in [0,194], denoting x position.
        """
        self._masks['base'][quad, asic, x, y] = 0
        return
    
        
    def unmask_pixel(self, quad, asic, x, y):
        """
        Mask a single pixel, or series of pixels. To do the latter, pass arrays
        as the arguments (even though the below says int).
        
        Parameters
        ----------
        quad : int
            [0,1,2,3], denoting the quad.
            
        asic : int
            Int in [0,7], denoting asic.
            
        x : int
            Int in [0,184], denoting x position.
        
        y : int
            Int in [0,194], denoting x position.
        """
        self._masks['base'][quad, asic, x, y] = 1
        return
    
        
    def mask_threshold(self, image, upper=None, lower=None):
        """
        Mask pixels by threshold values, in ADU. You must supply either `upper`
        or `lower`, but one one is required.
        
        Parameters
        ----------
        image : np.ndarray
            A shape (4, 8, 185, 388) array describing a CSPad image.
        
        upper : int
            Values greater than this are masked.
            
        lower : int
            Values lower than this are masked.
        """
        
        print("Masking pixels outside of [%s,%s]" % (str(lower), str(upper)))
        
        if (upper == None) and (lower == None):
            raise ValueError('Either `upper` or `lower` (or both) must be specified')
            
        if upper and lower:
            if upper <= lower:
                raise ValueError('Must have: `upper` > `lower` to threshold')
        
        self._check_image(image)
        
        m = self._blank_mask()
        ind = (image > upper) + (image < lower)
        m[ind] = 0
        
        self._inject_mask('threshold', m, override_previous=True)
        
        return
    
        
    def mask_nonbonded(self, nearest_neighbours=True):
        """
        Mask pixels on the CSPad that were never bonded.
        
        Optional Parameters
        -------------------
        nearest_neighbours : bool
            Also mask four of their nearest neighbours, which give anomoulous 
            responses.
        """
        
        print("Masking nonbonded pixels")
        
        m = self._blank_mask()
        
        for i in range(4):
            for j in range(16):
                for p in range(0, 185, 10):
                    
                    m[i,j,p,p] = 0
                    if nearest_neighbours:
                        if p == 0:
                            m[i,j,p+1,p] = 0
                            m[i,j,p,p+1] = 0
                        else:
                            m[i,j,p-1:p+2,p] = 0
                            m[i,j,p,p-1:p+2] = 0
                    
        self._inject_mask('nonbonded', m, override_previous=True)
        
        return
    
        
    def mask_borders(self, num_pixels=1):
        """
        Mask the border of each ASIC, to a width of `num_pixels`.
        
        Parameters
        ----------
        num_pixels : int
            The size of the border region to mask.
        """
        
        print("Masking %d pixels around the border of each 2x1" % num_pixels)
        
        n = int(num_pixels)        
        m = self._blank_mask()
        
        if (num_pixels < 0) or (num_pixels > 194):
            raise ValueError('`num_pixels` must be >0, <194')
        
        for i in range(4):
            for j in range(16):
                
                # mask along the y-dim
                m[i,j,:,0:n] = np.bool(False)
                m[i,j,:,194-n:194] = np.bool(False)
                
                # mask along the x-dim
                m[i,j,0:n,:] = np.bool(False)
                m[i,j,185-n:185,:] = np.bool(False)
                
                # # mask a bar along y in the middle of the 2x1
                # m[i,j,:,194-n:194+n] = np.bool(False)
        
        self._inject_mask('border', m, override_previous=True)
        
        return
    
        
    def mask_row13(self):
        
        print("Masking row 13")
        
        #raise NotImplementedError()
        print("Warning: row 13 masking is untested, tell the dev team if you need it.")
        
        #this is for masking out row13 of the CSPAD
        col=181
        for i in range(8):
            self.automask[:,col]=1
            col+= 194
            
        self._inject_mask('row13', m, override_previous=True)
    
    # ----------
        
    def merge(self, *args):
        """
        Merge two or more masks, masking with an OR operator for masked pixels.
        """
        
        for mask in args:
            for mtype in list(mask._masks.keys()):
                if mtype in list(self._masks.keys()):
                    self._masks[mtype] = np.logical_not( np.logical_or(self._masks[mtype], 
                                                           mask._masks[mtype]) )
                else:
                    self._masks[mtype] = mask._masks[mtype]
        
        return
    
        
    def save(self, filename, fmt='pypad'):
        """
        Save the PadMask object to one of many possible formats:
        
        -- pypad : An hdf5 format that includes all metadata associated with the
                   mask. Not read by other software (suffix: .mask).
                   
        -- cheetah : Stores the mask as a two-dimensional array in an HDF5 format.
                     Easily read into Cheetah and Thor (suffix: .h5).
        
        Parameters
        ----------
        filename : str
            The name of the file to write. This function will append an 
            appropriate suffix if none is provided.
            
        fmt : str, {'pypad', 'cheetah', 'thor'}
            The format to save in. See above for documentation.
        """
        
        if fmt == 'pypad':
            if not filename.endswith('.mask'):
                filename += '.mask'
            
            f = h5py.File(filename, 'w')
            for k in list(self._masks.keys()):
                f['/' + k] = self._masks[k]
            f.close()
            
        
        elif fmt == 'thor':
            if not filename.endswith('.h5'):
                filename += '.h5'


            # super ghetto, but it was the easy way out. sorry.
            # this converts first to cheetah 2d format, then
            # to Thor 1d format. I used it b/c I could be sure it
            # worked just by c/p code...

            itx = self.mask2d

            if not itx.shape == (1480, 1552):
                 raise ValueError('`itx` argument array incorrect shape! Must be:'
                                 ' (1480, 1552), got %s.' % str(itx.shape))

            flat_itx = np.zeros(1480 * 1552, dtype=itx.dtype)

            for q in range(4):
                for twoXone in range(8):

                    # extract the cheetah itx
                    x_start = 388 * q
                    x_stop = 388 * (q+1)

                    y_start = 185 * twoXone
                    y_stop = 185 * (twoXone + 1)

                    # each sec is a ASIC, both belong to the same 2x1
                    sec1, sec2 = np.hsplit(itx[y_start:y_stop,x_start:x_stop], 2)

                    # determine the positions of the flat array to put intens data in
                    n_ASIC_pixels = 185 * 194
                    flat_start = (q * 8 + twoXone) * (n_ASIC_pixels * 2) # 2x1 index X px in 2x1

                    # inject them into the thor array
                    flat_itx[flat_start:flat_start+n_ASIC_pixels] = sec1.flatten()
                    flat_itx[flat_start+n_ASIC_pixels:flat_start+n_ASIC_pixels*2] = sec2.flatten()
            
            f = h5py.File(filename, 'w')
            f['/mask'] = flat_itx
            f.close()

            
        elif fmt in ['cheetah', 'twod']:
            if not filename.endswith('.h5'):
                filename += '.h5'
            
            f = h5py.File(filename, 'w')
            f['/data/data'] = self.mask2d
            f.close()
            
            
        else:
            raise IOError('Unrecognized format for PadMask: %s. Should be one of'
                          ' {"pypad", "thor", "cheetah", "twod"}' % fmt)
        
        print("Wrote: %s" % filename)
        return
    
    
    @classmethod    
    def load(cls, filename):
        """
        Load a saved mask. Can be one of many formats:
        
            -- pypad .mask
            -- cheetah .h5
        
        Parameters
        ----------
        filename : str
            The name of the file to read.
        """
        
        m = cls()
        
        if filename.endswith('.mask'):
            f = h5py.File(filename, 'r')
            for k in f:
                m._masks[k] = np.array(f[k])
            f.close()
            
        elif filename.endswith('.h5'):
            
            try:
                f = h5py.File(filename, 'r')
                d = np.array( f['/data/data'] )
                assert d.shape == (1480, 1552)
                f.close()
            except:
                raise IOError('Cannot read data inside: %s. Either data is '
                              'corrupt or not in cheetah format [in /data/data'
                              ' and shape (1480, 1552)]' % filename)
            
            m._masks['cheetah'] = np.array( read.enforce_raw_img_shape(d) )
        
        else:
            raise IOError('Can only read files with {.mask, .h5} format -- got: %s' % filename)
            
        return m
    
    
class MaskGUI(object):

    def __init__(self, raw_image, mask=None, filename='my_mask', fmt='pypad'):
        """
        Instantiate an interactive masking session.
        
        Parameters
        ----------
        raw_image : np.ndarray
            A shape (4, 8, 185, 388) array containing a reference image that
            the user will use to guide their masking.
            
        mask : padmask.PadMask
            A PadMask object to modify. If `None` (default), generate a new 
            mask.
            
        filename : str
            The name of the file to generate at the end of the session.
            
        fmt : str
            The file format of `filename` to write.
        """
        
        
        self.print_gui_help()
        
        self.filename = filename
        self.file_fmt = fmt
        
        if not raw_image.shape == (4, 16, 185, 194):
            raise ValueError("`raw_image` must have shape: (4, 16, 185, 194)")
            
        if mask == None:
            self.mask = PadMask()            
        elif isinstance(mask, PadMask):
            self.mask = mask
            
        else:
            raise TypeError('`mask` argument must be a pypad.mask.PadMask object')
        
        
        # inject a new mask type into our PadMask obj
        m = self.mask._blank_mask()
        
        if not 'manual' in list(self.mask._masks.keys()):
            self.mask._inject_mask('manual', m)

        # deal with negative values
        if not 'negatives' in list(self.mask._masks.keys()):
            self.mask._inject_mask('negatives', m.copy())
            self.mask._masks['negatives'][raw_image <= 0.0] = 0
            print("Masked: %d negative pixels" % np.sum(np.logical_not(self.mask._masks['negatives'])))
        
        
        # we're going to plot the log of the image, so do that up front
        self.raw_image_4d = raw_image
        self.raw_image = utils.flatten_2x1s(raw_image)
        
        self.log_image = self.raw_image.copy()
        self.log_image[self.log_image < 0.0] = 0.0
        self.log_image = np.log10(self.log_image + 1.0)
        
        
        # populate an array containing the indices of all pixels in the image
        mg = np.meshgrid( np.arange(self.raw_image.shape[0]),
                          np.arange(self.raw_image.shape[1]) )
        self.points = np.vstack((mg[0].flatten(), mg[1].flatten())).T
        
        
        # create a colormap with masked pixels clearly highlighted
        self.palette = plt.cm.PuOr_r # reversed purple-orange -- base cm
        self.palette.set_under(color='green')
        
                
        # draw the main GUI, which is an image that can be interactively masked
        plt.figure(figsize=(9,6))
        self.ax = plt.subplot(111)
        
        self.im = self.ax.imshow( (self.log_image * self.mask.mask2d) - 1e-10, cmap=self.palette,
                                  origin='lower', interpolation='nearest', vmin=1e-10, aspect=1,
                                  extent=[0, self.log_image.shape[0], 0, self.log_image.shape[1]] )
        
        self.lc, = self.ax.plot((0,0),(0,0),'-+m', linewidth=1, markersize=8, markeredgewidth=1)
        self.lm, = self.ax.plot((0,0),(0,0),'-+m', linewidth=1, markersize=8, markeredgewidth=1)
        
        self.line_corner = (0,0)
        self.xy = None
        self.lines_xy = None
        self.single_px = None # for masking single pixels
        
        self.colorbar = plt.colorbar(self.im, pad=0.01)
        self.colorbar.set_label(r'$\log_{10}$ Intensity')

        cidb = plt.connect('button_press_event',  self.on_click)
        cidk = plt.connect('key_press_event',     self.on_keypress)
        cidm = plt.connect('motion_notify_event', self.on_move)
        
        plt.xlim([0, self.log_image.shape[0]])
        plt.ylim([0, self.log_image.shape[1]])
        
        self.ax.get_xaxis().set_ticks([])
        self.ax.get_yaxis().set_ticks([])

        
        # add toggle buttons that allow the user to turn on and off std masks
        # I used to have this in its own nice function, but MPL didn't like 
        # that for some reason... there is probably a better way, I just dont
        # know the innerds of MPL enough --TJL
        
        
        axcolor = 'lightgoldenrodyellow'

        ax1 = plt.axes([0.04, 0.7, 0.12, 0.08])                       
        self.b1 = ToggleButton(ax1, 'nonbonded', color=axcolor, hovercolor='0.975')
        self.b1.on_turned_on(self.mask.mask_nonbonded)
        self.b1.on_turned_off(self.mask.remove_mask, 'nonbonded')
        self.b1.on_turned_on(self.update_image)
        self.b1.on_turned_off(self.update_image)
        
        ax2 = plt.axes([0.04, 0.6, 0.12, 0.08])                       
        self.b2 = ToggleButton(ax2, 'row 13', color=axcolor, hovercolor='0.975')
        self.b2.on_turned_on(self.mask.mask_row13)
        self.b2.on_turned_off(self.mask.remove_mask, 'row13')
        self.b2.on_turned_on(self.update_image)
        self.b2.on_turned_off(self.update_image)
        
        ax3 = plt.axes([0.04, 0.5, 0.12, 0.08])                       
        self.b3 = ToggleButton(ax3, 'borders', color=axcolor, hovercolor='0.975')
        self.b3.on_turned_on(self._set_borderwidth)
        self.mask_border_cid = self.b3.on_turned_on(self.mask.mask_borders)
        self.b3.on_turned_off(self.mask.remove_mask, 'border')
        self.b3.on_turned_on(self.update_image)
        self.b3.on_turned_off(self.update_image)
        
        ax4 = plt.axes([0.04, 0.4, 0.12, 0.08])                       
        self.b4 = ToggleButton(ax4, 'threshold', color=axcolor, hovercolor='0.975')
        self.b4.on_turned_on(self._set_threshold)
        self.mask_threshold_cid = self.b4.on_turned_on(self.mask.mask_threshold, self.raw_image_4d, None, None)
        self.b4.on_turned_off(self.mask.remove_mask, 'threshold')
        self.b4.on_turned_on(self.update_image)
        self.b4.on_turned_off(self.update_image)
                           
        plt.show()
        
        return
        
        
    def _set_threshold(self):
        print("\n --- Enter threshold values --- ")
        self.lower_thld = float( input('Enter lower threshold: ') )
        self.upper_thld = float( input('Enter upper threshold: ') )
        self.b4.onstate_exargs[ self.mask_threshold_cid ] = (self.raw_image_4d, self.upper_thld, self.lower_thld)
        return
    

    def _set_borderwidth(self):
        print("\n --- Enter the desired border width --- ")
        raw_in = input('Size of border (in pixels) [1]: ')
        if raw_in == '':
            self.borderwidth = 1
        else:
            self.borderwidth = int( raw_in )
        self.b3.onstate_exargs[ self.mask_border_cid ] = (self.borderwidth,)
        return
    
    
    def update_image(self):
        self.im.set_data( (self.log_image * self.mask.mask2d) - 1e-10 )
        return

    
    def on_click(self, event):
        
        # for WHATEVER reason, the imshow drawing is stretched incorrectly
        # such that pixel positions for x/y are off by the ratio used below
        # ... this is likely due to me not understanding MPL, hence this hack
        # -- TJL
        ratio = float(self.log_image.shape[0]) / float(self.log_image.shape[1])
        x_coord = event.xdata / ratio
        y_coord = event.ydata * ratio
         
        # if a button that is *not* the left click is pressed
        if event.inaxes and (event.button is not 1):

            # save the points for masking in pixel coordinates
            if self.xy != None:
                self.xy = np.vstack(( self.xy, np.array([int(x_coord), 
                                                         int(y_coord)]) ))
            else:
                self.xy = np.array([int(x_coord), int(y_coord)])

            # save the points for drawing the lines in MPL coordinates
            if self.lines_xy != None:
                self.lines_xy = np.vstack(( self.lines_xy, np.array([int(event.xdata), 
                                                                     int(event.ydata)]) ))
            else:
                self.lines_xy = np.array([int(event.xdata), int(event.ydata)])
            
            self.lc.set_data(self.lines_xy.T) # draws lines
            self.line_corner = (int(event.xdata), int(event.ydata))
            
        # if the left button is pressed
        elif event.inaxes and (event.button is 1):
            self.single_px = (int(x_coord), int(y_coord))
            print("Selected: (%s, %s)" % self.single_px)        
        
        return


    def on_keypress(self, event):

        # mask or unmask
        if event.key in ['m', 'u']:
            
            if self.xy == None:
                print("No area selected, mask not changed.")
            else:
            
                # print "Masking region inside:"
                # print self.xy
           
                # wrap around to close polygon
                self.xy = np.vstack(( self.xy, self.xy[0,:] ))
                path = Path(self.xy)
                in_area = path.contains_points(self.points+0.5)
                inds = self.points[in_area]
            
                #print self.xy
                #print inds

                # if we're going to mask, mask
                if event.key == 'm':
                    print('Masking convex area...')
                    x = self._conv_2dinds_to_4d(inds)
                    self.mask._masks['manual'][x[:,0],x[:,1],x[:,2],x[:,3]] = 0
                
                # if we're unmasking, unmask
                elif event.key == 'u':
                    print('Unmasking convex area...')
                    x = self._conv_2dinds_to_4d(inds)
                    self.mask._masks['manual'][x[:,0],x[:,1],x[:,2],x[:,3]] = 1
            
                # draw and reset
                self.update_image()
                self._reset()
                plt.draw()
        

        # reset all masks
        elif event.key == 'r':
            print('Unmasking all')
            self.mask._masks['manual'] = self.mask._blank_mask()
            self.update_image()
            self._reset()
            #self.im.autoscale()
            plt.draw()
        
        # toggle selection    
        elif event.key == 't':
            
            if self.single_px != None:            
                x = self._conv_2dinds_to_4d( np.array(self.single_px)[None,:] )
                x = x.flatten()
                if self.mask.mask[x[0],x[1],x[2],x[3]] == 0:
                    print("Unmasking single pixel:", self.single_px)
                    self.mask._masks['manual'][x[0],x[1],x[2],x[3]] = 1
                else:
                    print("Masking single pixel:", self.single_px)
                    self.mask._masks['manual'][x[0],x[1],x[2],x[3]] = 0
                    
            else:
                print("No single pixel selected to toggle. Click a pixel and")
                print("    press `t` to toggle the mask on that pixel.")
            
            self.update_image()
            self._reset()
            #self.im.autoscale()
            plt.draw()
            
        # clear mouse selection
        elif event.key == 'x':
            print("Reset selections")
            self._reset()
            
           
        # save and exit
        elif event.key == 'w':
            self.mask.save(self.filename, fmt=self.file_fmt)
            plt.close()
            return
          
        # exit w/o saving
        elif event.key == 'q':
            print('Exiting without saving...')
            plt.close()
            return
          
        # else:
        #     print "Could not understand key: %s" % event.key
        #     print "Valid options: {m, u, r, k, q}"
            
        return
    

    def on_move(self, event):
        if not event.inaxes: return
        xm, ym = int(event.xdata), int(event.ydata)
        
        # update the line positions
        if self.line_corner != (0,0):
            self.lm.set_data((self.line_corner[0],xm), (self.line_corner[1],ym))
            plt.draw()
            
        return
    
        
    def _reset(self):
        self.single_pixel = None
        self.xy = None
        self.lines_xy = None
        self.lc.set_data([], [])
        self.lm.set_data([], [])
        self.line_corner = (0, 0)
        return
    
        
    def _conv_2dinds_to_4d(self, inds):
        """
        Convert indices in a 2d Cheetah array to (4,16,185,194).
        
        Parameters
        ----------
        inds : np.ndarray, int
            An N x 2 array, where the first column indexes x on the 2d image,
            and the second column indexes y.
            
        Returns
        -------
        inds_4d : np.ndarray, int
            An N x 4 array, with each column indexing quads/asics/y/x,
        """
        
        inds_4d = np.zeros((inds.shape[0], 4), dtype=np.int32)
        
        # index each asic, in the correct order
        of64 = (inds[:,1] / 185) * 2 + (inds[:,0] / 388) * 16 + (inds[:,0] / 194) % 2
        assert np.all(of64 < 64)
        
        # quads / asics
        # print 'masking in ASICs:', inds, of64
        inds_4d[:,0] = of64 / 16
        inds_4d[:,1] = of64 % 16
        
        # x / y : note the image is displayed transposed
        inds_4d[:,2] = (inds[:,1] % 185)
        inds_4d[:,3] = (inds[:,0] % 194)
        
        return inds_4d
    
        
    def print_gui_help(self):
        
        print()
        print()
        print("   --- WELCOME TO PYPAD's INTERACTIVE MASKING ENVIRONMENT --- ")
        print()
        print(" Green pixels are masked.")
        print()
        print(" Keystrokes")
        print(" ----------")
        print(" m : mask               u : unmask            r : reset ")
        print(" x : clear selection    w : save & exit       t : toggle pixel")
        print(" q : exit w/o saving")
        print()
        print(" Mouse")
        print(" -----")
        print(" Right click on three or more points to draw a polygon around a")
        print(" set of pixels. Then press `m` or `u` to mask or unmask that area.")
        print()
        print(" You can also mask/unmask single pixels by clicking on them with")
        print(" the mouse and pressing `t` to toggle the mask state.")
        print()
        print(" Toggle Buttons (left)")
        print(" ---------------------")
        print(" nonbonded : Mask nonbonded pixels, and their nearest neighbours.")
        print("             These pixels aren't even connected to the detector.")
        print() 
        print(" row 13    : In some old experiments, row 13 on one ASIC was ")
        print("             busted -- mask that (will be clear if this is needed)")
        print() 
        print(" threshold : You will be prompted for an upper and lower limit --")
        print("             pixels outside that range are masked. Units are ADUs")
        print("             and you can set only a lower/upper limit by passing")
        print("             'None' for one option.")
        print() 
        print(" borders   : Mask the borders of each ASIC. These often give")
        print("             anomoulous responses. Recommended to mask one pixel")
        print("             borders at least.")
        print("")
        print("                          ----- // -----")
        print()
        
        
        
        
        

