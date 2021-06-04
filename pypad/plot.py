

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
from matplotlib.widgets import Button, Slider
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
            for cid, func in self.onstate_observers.items():
                if len(self.onstate_exargs[cid]) > 0:
                    func(*self.onstate_exargs[cid])
                else:
                    func()
        else:
            for cid, func in self.offstate_observers.items():
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
        

class TwoPanelCSPAD(object):
    """
    This is a base for a class that displays the CSPAD in the left hand panel,
    and the radial projections on the right hand side.
    """
    
    def __init__(self, image, cspad):
        
        self.image = image
        self.cspad = cspad
        
        # initialize the plot
        self.fig = plt.figure(figsize=(12,7))
        self.axL = plt.subplot(121)
        self.axR = plt.subplot(122)
        
        self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        
        # actually show the image!
        self.im = imshow_cspad( self.cspad(self.image), ax=self.axL, 
                                scrollable=True )
        self.update_image()
        
        return
    

    def show(self):
        """ display the image and wait for the user to quit out of it """
        plt.show()
        return
        
        
    def draw_cspad(self, image):
        """
        draw the lefthand panel
        """
        
        self.im.set_data(self.cspad(self.image))
        
        # plot beam center
        beam_center = plt_patches.Circle((1000,1000), 2, fill=True, lw=1, color='w')
        self.axL.add_patch(beam_center)
    
        return
    

    def draw_projection(self, image):
        """
        Draw a radial projection of the intensities, in the rightmost panel.
        """
        
        self.axR.cla()
        
        n_bins = 800
        
        # plot all quads
        bin_centers, a = self.cspad.intensity_profile(image, n_bins=n_bins)
        self.axR.plot(bin_centers, a / a.max(), color='orange', lw=4)
        
        # plot each quad individually
        for i in range(4):
            bin_centers, a = self.cspad.intensity_profile(image, n_bins=n_bins, quad=i)
            a /= a.max()
            a += 1.0 * i + 1.0
            self.axR.plot(bin_centers, a, color=quad_colors[i], lw=2)
        
        self.axR.set_xlabel('Radius [mm]')
        self.axR.set_ylabel('Intensity')
        self.axR.get_yaxis().set_ticks([])
        
        self.axR.set_ylim([-0.3, 5.3])
        
        self.axR.text(100, -0.2, 'All Quads')
        self.axR.text(100,  0.8, 'Quad 0')
        self.axR.text(100,  1.8, 'Quad 1')
        self.axR.text(100,  2.8, 'Quad 2')
        self.axR.text(100,  3.8, 'Quad 3')
        
        return
    

    def update_image(self, val=None):
        """
        If a slider gets moved, update the panels
        """
                
        self.draw_cspad(self.image)
        self.draw_projection(self.image)
        plt.draw()
        
        return
                
        
    def on_scroll(self, event):
        """
        when we scroll, zoom in on the CSPAD intensities
        """
        
        lims = self.im.get_clim()
        speed = 1.2

        if event.button == 'up':
            colmax = lims[1] / speed
        elif event.button == 'down':
            colmax = lims[1] * speed

        self.im.set_clim(lims[0], colmax)
        plt.draw()
        
        
class ManipTwoPanelCSPAD(TwoPanelCSPAD):
    """
    A TwoPanelCSPAD object with some additional interactive features,
    specifically:
    
    -- a dilation slider
    -- the ability to set the beam center by right-clicking
    """
    
    def __init__(self, image, cspad):
        
        super(ManipTwoPanelCSPAD, self).__init__(image, cspad)
        
        plt.subplots_adjust(bottom=0.25)
        
        # initial values
        self.dilation     = 0.0
        self.center       = [0.0, 0.0]
        
        self.axcolor = 'lightgoldenrodyellow'
                
        # build the dilation slider
        self.ax_dilt = plt.axes([0.20, 0.10, 0.60, 0.03], axisbg=self.axcolor)
        self.dilation_slider = Slider(self.ax_dilt, 'Dilate [mm]', 0.0, 10.0, valinit=5.0)
        self.dilation_slider.on_changed(self.update_image)
        
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        
        self.update_image()
        
        return
    
        
    def update_image(self, val=None):
        
        # apply any dilation
        if hasattr(self, 'dilation_slider'):
            delta_dilation = self.dilation_slider.val - self.dilation
            self.cspad.dilate(delta_dilation)
            self.dilation = self.dilation_slider.val
            if np.abs(delta_dilation) > 1e-8:
                print("Dilation set to: %.2f" % self.dilation)
        
        super(ManipTwoPanelCSPAD, self).update_image(val)
        
        return
    
        
    def on_click(self, event):
        """
        If the user clicks on the image
        """
        if event.inaxes and (event.button is not 1):

            # clicks on left panel -- set beam center
            # the center is always at (1000, 1000): so translate the entire 
            # cspad to move the center
            if event.inaxes == self.axL:
                self.center[0] += 1000.0 - event.xdata
                self.center[1] += 1000.0 - event.ydata
                delta_center = (1000.0 - event.xdata, 1000.0 - event.ydata)
                print("Shifting center in x/y by: (%.2f, %.2f)" % delta_center)
                
                offset = np.array(delta_center)[None,:] * 0.10992 # pxl --> mm
                
                if np.any( np.abs(self.cspad.quad_offset + offset) > 15.0 ):
                    print("Warning: center move placed image outside reasonable")
                    print("bounds: try again with.")
                else:
                    self.cspad.quad_offset += offset                    
                    self.update_image(event)
                
            # clicks on right panel -- set opt limits
            elif event.inaxes == self.axR:
                v_line = event.xdata
                print("Set peak limit at: %f" % v_line)
                self.regions.append(v_line)
                self.update_image(event)
                    
        return

  
