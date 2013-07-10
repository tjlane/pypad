

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
A library for plotting pretty images of all kinds.
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.patches as plt_patches

quad_colors = ['k', 'g', 'purple', 'b']

def sketch_2x1s(pixel_positions, mpl_axes=None):
    """
    Draw a rough sketch of the layout of the CSPAD

    Parameters
    ----------
    pixel_positions : np.ndarray
        The x,y,z coordinates of the pixels on the CSPAD
    """
    
    if pixel_positions.shape not in [(3,4,16,185,194), (2,4,16,185,194)]:
        raise ValueError('`pixel_positions` has incorrect shape: '
                         '%s' % str(pixel_positions.shape))

    if not mpl_axes:
        plt.figure()
        ax = plt.subplot(111)
    else:
        ax = mpl_axes

    for i in range(4):
        for j in range(pixel_positions.shape[2]):
            x = pixel_positions[0,i,j,:,:]
            y = pixel_positions[1,i,j,:,:]
            corners = np.zeros((5,2))

            corners[0,:] = np.array([ x[0,0],   y[0,0] ])     # bottom left
            corners[1,:] = np.array([ x[0,-1],  y[0,-1] ])    # bottom right
            corners[3,:] = np.array([ x[-1,0],  y[-1,0] ])    # top left
            corners[2,:] = np.array([ x[-1,-1], y[-1,-1] ])   # top right
            corners[4,:] = np.array([ x[0,0],   y[0,0] ])     # make rectangle

            ax.plot(corners[:,0], corners[:,1], lw=2, color=quad_colors[i])
            ax.scatter(x[0,0], y[0,0])
            
    beam_center = plt_patches.Circle((0, 0), 2, fill=True, lw=1, color='orange')
    ax.add_patch(beam_center)
            
    # mirror x axis for CXI convention
    if not ax.xaxis_inverted():
        ax.invert_xaxis()

    if mpl_axes:
        return ax
    else:
        plt.show()


def imshow_cspad(image, vmin=0, vmax=None, ax=None, show=False, scrollable=False):
    """
    Show an assembled image (e.g. from CSPad(raw_image) ) as it would be seen
    when viewed from upstream at CXI. CXI convention is that the plus-x direction
    is towards the hutch door, plus-y is upwards, and plus-z is the direction
    of the beam.
    
    Parameters
    ----------
    image : np.ndarray
        A two-dimensional assembled image
    
    Returns
    -------
    im : axes.imshow
        The imshow instance.
    """
    
    # to be consistent with CXI convention, we want +x going left, and +y up
    
    if ax == None:
        ax = plt.subplot(111)

    im = ax.imshow( image, origin='lower', vmin=vmin, vmax=vmax,
                    interpolation='nearest' )
                    
    if not ax.xaxis_inverted():
        ax.invert_xaxis()
        
    def on_scroll(event):
        lims = im.get_clim()
        speed = 1.2

        if event.button == 'up':
            colmax = lims[1] / speed
        elif event.button == 'down':
            colmax = lims[1] * speed

        im.set_clim(lims[0], colmax)
        plt.draw()
        
    if scrollable:
        ax.figure.canvas.mpl_connect('scroll_event', on_scroll)
        
    if show:
        plt.show()
    
    return im
    
    
class ToggleButton(Button, object):
    """
    A quick subclass of MPL's button that has a state, such that it can toggle.
    """
    
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


    def ignore(self, *args):
        """
        This function is necessary to ensure matplotlib compatability
        for versions less than 1.2 (tested for 1.0.0 and 1.1.1).
        """
        if hasattr( super( ToggleButton, self), 'ignore' ):
            return super( ToggleButton, self ).ignore(*args)
        else:
            return False
        
