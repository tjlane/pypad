
"""
A library for plotting pretty images of all kinds.
"""

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)

import numpy as np
import matplotlib.pyplot as plt


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
    
    quad_color = ['k', 'g', 'purple', 'b']

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

            ax.plot(corners[:,0], corners[:,1], lw=2, color=quad_color[i])
            ax.scatter(x[0,0], y[0,0])
            
            
    # mirror x axis for CXI convention
    ax.invert_xaxis()

    if mpl_axes:
        return ax
    else:
        plt.show()


def imshow_cspad(image, vmin=0, vmax=None, ax=None):
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
    ax : pyplot.axes
    im : axes.imshow
    """
    
    if ax == None:
        ax = plt.subplot(111)

    im = ax.imshow( image, origin='lower', vmin=vmin, vmax=vmax,
                    interpolation='nearest' )
    ax.invert_xaxis()
    return im

                                                    
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
            logger.info("Saved image: %s" % self.filename)
            
        elif event.key == 'r':
            logger.info("Reset plot")
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
        
