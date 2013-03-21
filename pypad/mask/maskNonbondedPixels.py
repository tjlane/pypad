#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

import os
import sys
import string
import numpy as N
import h5py as H
import matplotlib
import matplotlib.pyplot as P
from optparse import OptionParser

# Usage:
# First argument: input HDF5 file containing pixel mask in /data/data
# Second argument: output HDF5 file containing pixel mask with masked borders in /data/data

#inputfile = sys.argv[1]
#outputfile = sys.argv[2]

parser = OptionParser()
parser.add_option("-f", "--file", action="store", type="string", dest="fileName",
                                        help="input mask you wish to mask non-bonded pixels", metavar="FILENAME", default="")
parser.add_option("-o", "--output", action="store", type="string", dest="maskName",
                                        help="output mask you wish to save", metavar="MASKNAME", default="mask+borders.h5")
parser.add_option("-n", "--neighbors", action="store", type="int", dest="nearestNeighbors",
                                        help="number of nearest neighbors to mask for each non-bonded pixel (default: 0)", metavar="NPIXELS", default="0")

(options, args) = parser.parse_args()

inputfile = options.fileName
outputfile = options.maskName
write_dir = os.getcwd()

ROWS = 194
COLS = 185
RAW_DATA_LENGTH = [8*COLS,8*ROWS]

# Open input mask
f = H.File(inputfile,"r")
mask = N.array(f['/data/data'])
f.close()

masked = N.where(mask < 1)
print "%s contains %s masked pixels." % (inputfile, len(mask[masked]))

########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, filename):
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		self.cmax = self.inarr.max()
		self.cmin = self.inarr.min()
	
	def on_keypress(self,event):
		if event.key == 'p':
			if not os.path.exists(write_dir):
				os.mkdir(write_dir)
			pngtag = write_dir + "/%s.png" % (self.filename)	
			print "saving image as " + pngtag 
			P.savefig(pngtag)
		if event.key == 'r':
			colmin, colmax = self.orglims
			P.clim(colmin, colmax)
			P.draw()


	def on_click(self, event):
		if event.inaxes:
			lims = self.axes.get_clim()
			colmin = lims[0]
			colmax = lims[1]
			range = colmax - colmin
			value = colmin + event.ydata * range
			if event.button is 1 :
				if value > colmin and value < colmax :
					colmin = value
			elif event.button is 2 :
				colmin, colmax = self.orglims
			elif event.button is 3 :
				if value > colmin and value < colmax:
					colmax = value
			P.clim(colmin, colmax)
			P.draw()
				

	def draw_img(self):
		fig = P.figure()
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin = 'lower', vmax = self.cmax)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		P.show()


print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

maskImg = img_class(mask,inputfile)
maskImg.draw_img()

########################################################
# Function to recursively mask out nearest neighbors
########################################################
def maskNeighbors(matrix, rowpos, colpos, nn):
	matrix[rowpos,colpos] = 0
	if (nn > 0):
		# CHECK LOWER EDGE OF ASIC
		if (rowpos % COLS != 0):
			maskNeighbors(matrix, rowpos-1, colpos, nn-1)
		# CHECK UPPER EDGE OF ASIC
		if (rowpos % COLS != COLS-1):
			maskNeighbors(matrix, rowpos+1, colpos, nn-1)
		# CHECK LEFT EDGE OF ASIC
		if (colpos % ROWS != 0):
			maskNeighbors(matrix, rowpos, colpos-1, nn-1)
		# CHECK RIGHT EDGE OF ASIC
		if (colpos % ROWS != ROWS-1):
			maskNeighbors(matrix, rowpos, colpos+1, nn-1)
	return


# Create output array of type 16-bit integer with non-bonded pixels (located on the diagonal every 10n pixels)
outputmask = mask
for arow in N.arange(8):
	for acol in N.arange(8):
		for n in N.arange(0,COLS,10):
			outputmask[arow*COLS+n,acol*ROWS+n] = 0
			maskNeighbors(outputmask, arow*COLS+n, acol*ROWS+n, options.nearestNeighbors)

# Save output array to HDF5 file
h5output = H.File(outputfile,'w')
datagroup = h5output.create_group("data")
dataset = datagroup.create_dataset("data",RAW_DATA_LENGTH,dtype="int16")
dataset[...] = outputmask[:,:]
h5output.close()

maskedo = N.where(outputmask < 1)
print "Saved new mask as %s containing %s masked pixels." % (outputfile, len(outputmask[maskedo]))

outputmaskImg = img_class(outputmask, outputfile)
outputmaskImg.draw_img()
