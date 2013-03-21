
## my favorate functions just to simplify life

from numpy import *

# a 2D array where each pixel is set to its x coordinate
def make_xarr(nx,ny): 
	xarr = outer(arange(0-nx/2,nx-nx/2,1),ones(ny))
	return xarr

# a 2D array where each pixel is set to its y coordinate 
def make_yarr(nx,ny):    	
	yarr = outer(ones(nx),arange(0-ny/2,ny-ny/2,1))
	return yarr

# shift - a 2D version of numpy's roll
def array_shift(array,xshift=0,yshift=0):
	array = roll(array,xshift,0)
	array = roll(array,yshift,1)
	return array

## make an array with a circle set to one
def circle(nx, ny, rad=None, cenx=None, ceny=None, invert=0 ): 
    #set defaults
    if rad is None: rad = min(nx,ny)/2
    if cenx is None: cenx = nx/2
    if ceny is None: ceny = ny/2

    # define the circle
    x = outer(arange(0-cenx,nx-cenx,1),ones(ny))
    print x.size, x.shape
    y = outer(ones(nx),arange(0-ceny,ny-ceny,1))
    print y.size, y.shape
    dist = sqrt(x**2 + y**2)
    a = zeros([nx,ny])
    icirc = where(dist <= rad)
    a[icirc] = 1.
    return a
#end circle

