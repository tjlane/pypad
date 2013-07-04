

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
            
            
    # mirror x axis for CXI convention
    if not ax.xaxis_inverted():
        ax.invert_xaxis()

    if mpl_axes:
        return ax
    else:
        plt.show()


def imshow_cspad(image, vmin=0, vmax=None, ax=None, show=False):
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


class InteractiveImshow(object):
    """
    A brief extension to matplotlib's imshow that puts a colorbar next to 
    the image that you can click on to scale the maximum numeric value
    displayed.
    
    Based on code from pyana_misc by Ingrid Ofte.
    
    Parameters
    ----------
    inarr : np.ndarray
        The array to imshow()
        
    filename : {str, None}
        The filename to call the file if it is saved. If `None`, disable saving
        ability.
    """
    
    def __init__(self, inarr, filename=None, fig=None, ax=None):
        """
        Parameters
        ----------
        inarr : np.ndarray
            The array to imshow()

        filename : {str, None}
            The filename to call the file if it is saved. If `None`, disable saving
            ability.
            
        fig : pyplot.figure
            A figure object to draw on.
            
        ax : pyplot.axes
            An axes canvas to draw on.
        """
        self.inarr = inarr
        self.filename = filename
        self.fig = fig
        self.ax = ax
        self.cmax = self.inarr.max()
        self.cmin = self.inarr.min()
        self._draw_img()
        

    def _on_keypress(self, event):
        
        if event.key == 's':
            if not self.filename:
                self.filename = raw_input('Saving. Enter filename: ')
            plt.savefig(self.filename)
            print "Saved image: %s" % self.filename
            
        elif event.key == 'r':
            print "Reset plot"
            colmin, colmax = self.orglims
            plt.clim(colmin, colmax)
            plt.draw()
            

    def _on_click(self, event):
        if event.inaxes:
            lims = self.im.get_clim()
            colmin = lims[0]
            colmax = lims[1]
            rng = colmax - colmin
            value = colmin + event.ydata * rng
            if event.button is 1:
                if value > colmin and value < colmax :
                    colmax = value
            elif event.button is 2:
                colmin, colmax = self.orglims
            elif event.button is 3:
                if value > colmin and value < colmax:
                    colmix = value
            self.im.set_clim(colmin, colmax)
            plt.draw()
            
            
    def _on_scroll(self, event):
        lims = self.im.get_clim()
        speed = 1.1
        
        if event.button == 'up':
            colmax = lims[1] / speed
        elif event.button == 'down':
            colmax = lims[1] * speed
            
        self.im.set_clim(lims[0], colmax)
        plt.draw()
            

    def _draw_img(self):
        
        if not self.fig:
            self.fig = plt.figure()
            
        cid1 = self.fig.canvas.mpl_connect('key_press_event', self._on_keypress)
        cid2 = self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        cid3 = self.fig.canvas.mpl_connect('scroll_event', self._on_scroll)
        
        if not self.ax:
            self.ax = self.fig.add_subplot(111)
        
        self.im = self.ax.imshow(self.inarr, vmin=0.0, vmax=self.cmax, origin='lower')
        self.colbar = plt.colorbar(self.im, pad=0.01)
        self.orglims = self.im.get_clim()
        
