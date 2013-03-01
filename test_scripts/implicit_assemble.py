
import numpy as np
import matplotlib.pyplot as plt

from autogeom import cspad
from autogeom.utils import _assemble_implicit, load_raw_image, sketch_2x1s

# ----------------

print "loading raw image"
raw_img = load_raw_image('data/test_images/gold/gold_avg.npz')

# ----------------

print "assembling on alignment params"

cs1 = cspad.CSPad.from_dir('data/ex_params/')

bg = cs1.basis_repr

plt.figure()

quad_color = ['k', 'g', 'purple', 'b']

for i in range(32):
    twoXone = bg.grid_as_explicit(i)
    print twoXone.shape
    x = twoXone[:,:,0]
    y = twoXone[:,:,1]
    corners = np.zeros((5,2))
    
    corners[0,:] = np.array([ x[0,0],   y[0,0] ])     # bottom left
    corners[1,:] = np.array([ x[0,-1],  y[0,-1] ])    # bottom right
    corners[3,:] = np.array([ x[-1,0],  y[-1,0] ])    # top left
    corners[2,:] = np.array([ x[-1,-1], y[-1,-1] ])   # top right
    corners[4,:] = np.array([ x[0,0],   y[0,0] ])     # make rectangle
        
    plt.plot(corners[:,0], corners[:,1], lw=2, color=quad_color[i/8])
    plt.scatter(x[0,0], y[0,0])

plt.show()

# ------------------

print "pos", cs1.pixel_positions

print "assembling on CSPad basis grid"
z = _assemble_implicit(cs1.pixel_positions, raw_img)
print z

plt.imshow(z.T, interpolation='nearest', vmin=0, vmax=100)
plt.show()

# ------------------

# print "testing quad rendering is correct"
# 
# test_img = np.zeros((4, 8, 185, 388))
# test_img[2,0,:,:] = 1.0
# 
# plt.figure()
# plt.imshow( cs1(test_img) )
# plt.show()

# ------------------

# print "assembling on metrology..."
# 
# met = cspad.Metrology('data/metrology/cspad_2011-08-10-Metrology.txt')
# print met.pixel_positions.shape
#
# 
# z = _assemble_implicit(pixels, raw_img)
# print z
# 
# plt.imshow(z.T, interpolation='nearest')
# plt.show()
# 
# xy = np.array((met.x_coordinates, met.y_coordinates))
# sketch_2x1s(xy)

# ---------------------