
"""
Provides a "mask" object for CSPads
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.nxutils import points_inside_poly
from matplotlib.widgets import Button

from pypad import utils


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
        m = np.product( np.array(self._masks.values()), axis=0 )
        assert m.shape == (4,8,185,388)
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
        return self._masks.keys()
    
    
    def remove_mask(self, mask_name):
        """
        Remove a mask that has been applied.
        """
        
        if mask_name == 'base': # this one is special
            self._masks['base'] = self._blank_mask()
        elif not mask_name in self._masks.keys():
            raise KeyError('Mask: %s not applied' % mask_name)
        else:
            x = self._masks.pop(mask_name)
            
        print "Removed mask: %s" % mask_name
            
        return
    
        
    def _inject_mask(self, mask_name, mask):
        """
        Add a new kind of mask to this mask object. Provides some typechecking.
        
        All this really does is deposit `mask` into dict self._masks with key
        `mask_name`.
        """
        
        assert type(mask_name) == str
        assert mask.shape == (4, 8, 185, 388)
        
        if not mask.dtype == np.bool:
            mask = mask.astype(np.bool)
            
        self._masks[mask_name] = mask
        
        return
        
        
    def _check_image(self, image):
        """
        Sanity check on `image`.
        """
        if not image.shape == (4, 8, 185, 388):
            raise ValueError('`image` must be shape (4, 8, 185, 388), got '
                             '%s' % str(image.shape))
        return
        
        
    def _blank_mask(self):
        """
        Utility function that just returns a blank mask.
        """
        return np.ones((4, 8, 185, 388), dtype=np.int32)
    
        
    # ----------
    # below we provide many methods of the form PadMask.mask_*(), which allow
    # one to mask pixels via various criteria (indicated by the *)
    
    # to add a new kind of mask, make a new method here. Follow mask_threshold
    # as a template
    
    def mask_pixel(self, quad, two_by_one, x, y):
        """
        Mask a single pixel, or series of pixels. To do the latter, pass arrays
        as the arguments (even though the below says int).
        
        Parameters
        ----------
        quad : int
            [0,1,2,3], denoting the quad.
            
        two_by_one : int
            Int in [0,7], denoting 2x1.
            
        x : int
            Int in [0,184], denoting x position.
        
        y : int
            Int in [0,388], denoting x position.
        """
        self._masks['base'][quad, two_by_one, x, y] = 0
        return
    
        
    def unmask_pixel(self, quad, two_by_one, x, y):
        """
        Mask a single pixel, or series of pixels. To do the latter, pass arrays
        as the arguments (even though the below says int).
        
        Parameters
        ----------
        quad : int
            [0,1,2,3], denoting the quad.
            
        two_by_one : int
            Int in [0,7], denoting 2x1.
            
        x : int
            Int in [0,184], denoting x position.
        
        y : int
            Int in [0,388], denoting x position.
        """
        self._masks['base'][quad, two_by_one, x, y] = 1
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
        
        print "Masking pixels outside of [%s,%s]" % (str(lower), str(upper))
        
        if (upper == None) and (lower == None):
            raise ValueError('Either `upper` or `lower` (or both) must be specified')
        
        self._check_image(image)
        
        m = self._blank_mask()
        ind = (image > upper) + (image < lower)
        m[ind] = np.bool(False)
        
        self._inject_mask('threshold', m)
        
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
        
        print "Masking nonbonded pixels"
        
        m = self._blank_mask()
        
        for i in range(4):
            for j in range(8):
                for p in range(0, 185, 10):
                    
                    if nearest_neighbours:                      
                        xl = max(0, p-1)
                        xh = min(184, p+1)
                        yl = max(0, p-1)
                        yh = min(387, p+1)
                        m[xl:xh,yl:yh] = 0
                        
                    else:
                        m[p,p] = 0
                    
        self._inject_mask('nonbonded', m)
        
        return
    
        
    def mask_borders(self, num_pixels=1):
        """
        Mask the border of each ASIC, to a width of `num_pixels`.
        
        Parameters
        ----------
        num_pixels : int
            The size of the border region to mask.
        """
        
        print "Masking %d pixels around the border of each 2x1" % num_pixels
        
        n = int(num_pixels)        
        m = self._blank_mask()
        
        if (num_pixels < 0) or (num_pixels > 194):
            raise ValueError('`num_pixels` must be >0, <194')
        
        for i in range(4):
            for j in range(8):
                
                # mask along the y-dim
                m[i,j,:,0:n] = np.bool(False)
                m[i,j,:,388-n:388] = np.bool(False)
                
                # mask along the x-dim
                m[i,j,0:n,:] = np.bool(False)
                m[i,j,185-n:185,:] = np.bool(False)
                
                # mask a bar along y in the middle of the 2x1
                m[i,j,:,194-n:194+n] = np.bool(False)
        
        self._inject_mask('border', m)
        
        return
    
        
    def mask_row13(self):
        
        print "Masking row 13"
        
        raise NotImplementedError()
        
        #this is for masking out row13 of the CSPAD
        col=181
        for i in range(8):
            self.automask[:,col]=1
            col+= 194
            
        self._inject_mask('row13', m)
    
    # ----------
        
    def merge(self, *args):
        """
        Merge two or more masks, masking with an OR operator for masked pixels.
        """
        
        for mask in args:
            for mtype in mask._masks.keys():
                if mtype in self._masks.keys():
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
                   
        -- twod  : Stores the mask as a two-dimensional array in an HDF5 format.
                   Easily read into Cheetah and Odin (suffix: .h5).
        
        Parameters
        ----------
        filename : str
            The name of the file to write. This function will append an 
            appropriate suffix if none is provided.
            
        fmt : str, {'pypad', 'twod'}
            The format to save in. See above for documentation.
        """
        
        if fmt == 'pypad':
            if not filename.endswith('.mask'):
                filename += '.mask'
            
            f = h5py.File(filename, 'w')
            for k in self._masks.keys():
                f['/' + k] = self._masks[k]
            f.close()
            
            
        elif fmt == 'twod':
            if not filename.endswith('.h5'):
                filename += '.h5'
                
                # need jonas to dbl check this is right for Cheetah
                f = h5py.File(filename, 'w')
                f['/data'] = self.mask
                f.close()
            
            
        else:
            raise IOError('Unrecognized format for PadMask: %s. Should be one of'
                          ' {"pypad", "twod"}' % fmt)
        
        print "Wrote: %s" % filename
        return
    
    
    @classmethod    
    def load(cls, filename):
        """
        Load a saved mask
        
        Parameters
        ----------
        filename : str
            The name of the file to read.
        """
        
        if not filename.endswith('.mask'):
            raise IOError('Can only read files with .mask format -- got: %s' % filename)
            
        m = cls()
        
        f = h5py.File(filename, 'r')
        for k in f:
            m._masks[k] = np.array(f[k])
        f.close()
        
        return
        
        
class ToggleButton(Button):
    
    def __init__(self, ax, label, default_state=False, **kwargs):
        """
        Create a matplotlib Button that can "toggle".
        
        When pressed, it flips state from either "on" to "off" of vice versa, 
        calling `on_func` or `off_func` as appropriate.
        
        Parameters
        ----------
        ax : matplotlib.axes
            The axes to draw the button on.
            
        label : str
            The name that appears in the middle of the button.
        
        default_state : bool
            Start off (False, default) or on (True).
        """
        
        # call Button(ax, label)
        super( ToggleButton, self ).__init__(ax, label, **kwargs)
        
        if not type(default_state) == bool:
            raise TypeError('`default_state` must be type: bool')
        self.on = default_state # keep track of state
        
        # we split observers into onstate/offstate_observers
        del self.cnt
        del self.observers
        
        self.cnt_on = 0
        self.onstate_observers  = {}
        
        self.cnt_off = 0
        self.offstate_observers = {}
        
        self.onstate_exargs  = {}
        self.offstate_exargs = {}
        
        return
    
        
    # override Button's release
    def _release(self, event):
        if self.ignore(event):
            return
        if event.canvas.mouse_grabber != self.ax:
            return
        event.canvas.release_mouse(self.ax)
        if not self.eventson:
            return
        if event.inaxes != self.ax:
            return
            
        # call either the on or off function
        if not self.on:
            for cid, func in self.onstate_observers.iteritems():
                if len(self.onstate_exargs[cid]) > 0:
                    func(*self.onstate_exargs[cid])
                else:
                    func()
        else:
            for cid, func in self.offstate_observers.iteritems():
                if len(self.offstate_exargs[cid]) > 0:
                    func(*self.offstate_exargs[cid])
                else:
                    func()
                
        # and toggle!
        self.on = not(self.on)

        
    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            if cid in self.onstate_observers:
                del self.onstate_observers[cid]
            else:
                del self.offstate_observers[cid]
        except KeyError:
            pass
    
        
    def on_turned_on(self, func, *args):
        """
        Call function `func` in the on state.
        """
        cid = self.cnt_on
        self.onstate_observers[cid] = func
        self.onstate_exargs[cid] = args
        self.cnt_on += 1
        return cid
    
        
    def on_turned_off(self, func, *args):
        """
        Call function `func` in the off state.
        """
        cid = self.cnt_off
        self.offstate_observers[cid] = func
        self.offstate_exargs[cid] = args
        self.cnt_off += 1
        return cid
    
    
class MaskGUI(object):

    def __init__(self, raw_image, mask=None, filename='my_mask', fmt='pypad'):
        
        if not raw_image.shape == (4, 8, 185, 388):
            raise ValueError("`raw_image` must have shape: (4, 8, 185, 388)")
            
        if mask == None:
            self.mask = PadMask()
        elif isinstance(mask, PadMask):
            self.mask = mask
        else:
            raise TypeError('`mask` argument must be a pypad.padmask.PadMask object')
        
        
        # inject a new mask type into our PadMask obj
        m = self.mask._blank_mask()
        self.mask._inject_mask('manual', m)
        
        
        # deal with negative values
        self.mask._inject_mask('negatives', m.copy())
        self.mask._masks['negatives'][raw_image < 0.0] = 0
        print "Masked: %d negative pixels" % np.sum(np.logical_not(self.mask._masks['negatives']))
        
        
        # we're going to plot the log of the image, so do that up front
        self.raw_image = utils.flatten_2x1s(raw_image)
        
        self.log_image = self.raw_image.copy()
        self.log_image[self.log_image < 0.0] = 0.0
        self.log_image = np.log10(self.log_image + 1.0)
        
        
        # populate an array containing the indices of all pixels in the image
        mg = np.meshgrid( np.arange(self.raw_image.shape[0]),
                          np.arange(self.raw_image.shape[1]) )
        self.points = np.vstack((mg[0].flatten(), mg[1].flatten())).T
        
        
        palette = plt.cm.jet
        palette.set_bad('w',1.0)
                
                
        # draw the main GUI, which is an image that can be interactively masked
        plt.figure()
        
        self.ax = plt.subplot(111)
        self.im = plt.imshow( (self.log_image * self.mask.mask2d).T, cmap=palette,
                              origin='lower', interpolation='nearest', vmin=0, 
                              extent=[0, self.log_image.shape[0], 0, self.log_image.shape[1]] )

        plt.title('Press m : mask || u : unmask || r : reset || k : save & exit || q : exit w/o saving')

        self.lc, = self.ax.plot((0,0),(0,0),'-+m', linewidth=1, markersize=8, markeredgewidth=1)
        self.lm, = self.ax.plot((0,0),(0,0),'-+m', linewidth=1, markersize=8, markeredgewidth=1)
        
        self.line_corner = (0,0)
        self.xy = None

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
        b1 = ToggleButton(ax1, 'nonbonded', color=axcolor, hovercolor='0.975')
        b1.on_turned_on(self.mask.mask_nonbonded)
        b1.on_turned_off(self.mask.remove_mask, 'nonbonded')
        b1.on_turned_on(self.update_image)
        b1.on_turned_off(self.update_image)
        
        ax2 = plt.axes([0.04, 0.6, 0.12, 0.08])                       
        b2 = ToggleButton(ax2, 'row 13', color=axcolor, hovercolor='0.975')
        b2.on_turned_on(self.mask.mask_row13)
        b2.on_turned_off(self.mask.remove_mask, 'row13')
        b2.on_turned_on(self.update_image)
        b2.on_turned_off(self.update_image)
        
        ax3 = plt.axes([0.04, 0.5, 0.12, 0.08])                       
        b3 = ToggleButton(ax3, 'borders', color=axcolor, hovercolor='0.975')
        b3.on_turned_on(self.mask.mask_borders)
        b3.on_turned_off(self.mask.remove_mask, 'border')
        b3.on_turned_on(self.update_image)
        b3.on_turned_off(self.update_image)
        
        ax4 = plt.axes([0.04, 0.4, 0.12, 0.08])                       
        b4 = ToggleButton(ax4, 'threshold', color=axcolor, hovercolor='0.975')
        
        # interactive
        def mt():
            print " --- Enter threshold values --- "
            upper = raw_input('')
            return self.mask.mask_threshold, (upper, lower)
        
        b4.on_turned_on(mt)
        b4.on_turned_off(self.mask.remove_mask, 'threshold')
        b4.on_turned_on(self.update_image)
        b4.on_turned_off(self.update_image)
                           
        plt.show()
        
        return
    
    
    def update_image(self):
        self.im.set_data( (self.log_image * self.mask.mask2d).T )
        return

    
    def on_click(self, event):
         
        if not event.inaxes: return
        
        if self.xy != None:
            self.xy = np.vstack(( self.xy, np.array([int(event.xdata), int(event.ydata)]) ))
        else:
            self.xy = np.array([int(event.xdata), int(event.ydata)])
            
        self.lc.set_data(self.xy.T) # draws lines
        self.line_corner = (int(event.xdata), int(event.ydata))
        
        return


    def on_keypress(self, event):

        if event.key in ['m', 'u']:
           
            # wrap around to close polygon
            self.xy = np.vstack(( self.xy, self.xy[0,:] ))
            inds = self.points[points_inside_poly(self.points, self.xy)]

            # if we're going to mask, mask
            if event.key == 'm':
                print 'Masking convex area...'
                x = self._conv_2dinds_to_4d(inds)
                self.mask._masks['manual'][x[:,0],x[:,1],x[:,2],x[:,3]] = 0
                
            # if we're unmasking, unmask
            elif event.key == 'u':
                print 'Unmasking convex area...'
                x = self._conv_2dinds_to_4d(inds)
                self.mask._masks['manual'][x[:,0],x[:,1],x[:,2],x[:,3]] = 1
            
            # draw and reset
            self.update_image()

            self._reset()
            self.im.autoscale()
            plt.draw()

        elif event.key == 'r':
            print 'Unmasking all'
            
            self.mask._masks['manual'] = self.mask._blank_mask()
            
            self.update_image()
            
            self._reset()
            self.im.autoscale()
            plt.draw()
           
        elif event.key == 'k':
            self.mask.save(self.filename, fmt=self.file_fmt)
            plt.close()
            return
          
        elif event.key == 'q':
            print 'Exiting without saving...'
            plt.close()
            return
          
        else:
            print "Could not understand key: %d" % event.key
            print "Valid options: {m, u, r, k, q}"
            
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
        self.xy = None
        self.lc.set_data([], [])
        self.lm.set_data([], [])
        self.line_corner = (0, 0)
        return
    
        
    def _conv_2dinds_to_4d(self, inds):
        """
        Convert indices in a Cheetah array to (4,8,185,388).
        
        Parameters
        ----------
        inds : np.ndarray, int
            An N x 2 array, where the first column indexes x on the 2d image,
            and the second column indexes y.
            
        Returns
        -------
        inds_4d : np.ndarray, int
            An N x 4 array, with each column indexing quads/2x1/x/y,
        """
        
        inds_4d = np.zeros((inds.shape[0], 4), dtype=np.int32)
        
        # abs 2x1 index = x / num_x + y / num_y * 2x1s-in-x
        of32 = (inds[:,0] / 185) + (inds[:,1] / 388) * 8
        assert np.all(of32 < 32)
        
        # quads / 2x1s
        inds_4d[:,0] = of32 % 4
        inds_4d[:,1] = of32 / 4
        
        # x / y
        inds_4d[:,2] = inds[:,0] % 185
        inds_4d[:,3] = inds[:,1] % 388
        
        return inds_4d
        
        
