#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./viewAssembledSum.py -rxxxx -a (optional) -m 10000
# For details, type 
#	 python viewAssembledSum --help
# where rxxxx is the run number of hits and nonhits found using the hitfinder executable. 
# By default, this script looks into the h5 files that are in the appropriate rxxxx directory
#

import os
import sys
import string
import numpy as N
import h5py as H
import matplotlib
import matplotlib.pyplot as P
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", action="store", type="string", dest="fileName",
		                        help="input gain file you wish to mask", metavar="FILENAME", default="")
parser.add_option("-o", "--output", action="store", type="string", dest="maskName",
		                        help="output mask you wish to save (default = 'pixel_mask_FILENAME.h5')", metavar="MASKNAME", default="")
parser.add_option("-l", "--lowthreshold", action="store", type="float", dest="lowThreshold",
		                        help="threshold (ADU) below which masking is performed (default = 0.75)", metavar="THRESHOLD", default="0.75")
parser.add_option("-u", "--highthreshold", action="store", type="float", dest="highThreshold",
		                        help="threshold (ADU) above which masking is performed (default = 1.25)", metavar="THRESHOLD", default="1.25")

(options, args) = parser.parse_args()

########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cheetah_input_files/"
#write_dir = "/reg/neh/home/sellberg/NML-2013/analysis/cheetah_scripts/figures/"
write_dir = "/reg/neh/home/sellberg/NML-2013/analysis/masking_script/"

ROWS = 194
COLS = 185
RAW_DATA_LENGTH = [8*COLS,8*ROWS]

if (options.maskName == ""):
	outputfilename = "pixel_mask_%s-%sADUs_%s" % (options.lowThreshold, options.highThreshold, options.fileName)
else:
	outputfilename = options.maskName

if (os.path.exists(options.fileName)):
	print options.fileName
	f = H.File(options.fileName,"r")
	d = N.array(f['/data/data'])
	f.close()
elif (os.path.exists(source_dir+options.fileName)):
	print source_dir+options.fileName
	f = H.File(source_dir+options.fileName,"r")
	d = N.array(f['/data/data'])
	f.close()
else:
	print "%s not found, aborting."
	sys.exit(1)

# Find pixels to mask
negindex = N.where(d < options.lowThreshold)
posindex = N.where(d > options.highThreshold)

# Create output array of type 16-bit integer with ones, fill in zeros
outputmask = N.ones(RAW_DATA_LENGTH, dtype="int16")
outputmask[negindex] = 0
outputmask[posindex] = 0
print "Masked out " + str(len(d[negindex])+len(d[posindex])) + " pixels."

# Save output array to HDF5 file
h5output = H.File(outputfilename, 'w')
datagroup = h5output.create_group("data")
dataset = datagroup.create_dataset("data",RAW_DATA_LENGTH,dtype="int16")
dataset[...] = outputmask[:,:]
h5output.close()
print "Saved mask as %s" % (outputfilename)

########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, filename):
		self.inarr = inarr
		self.filename = filename
		self.cmax = self.inarr.max()
		self.cmin = self.inarr.min()
	
	def on_keypress(self,event):
		if event.key == 'p':
			pngtag = write_dir + "%s.png" % (self.filename)
			print "saving image as " + pngtag 
			P.savefig(pngtag)
		if event.key == 'r':
			colmin = self.cmin
			colmax = self.cmax
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
		P.rc('image',origin='lower')
		self.axes = P.imshow(self.inarr, vmin = options.lowThreshold, vmax = options.highThreshold)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		#P.show()
		pngtag = write_dir + "%s.png" % (self.filename)
		print "saving image as " + pngtag 
		P.savefig(pngtag)


print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

currImg = img_class(d*outputmask, outputfilename)
currImg.draw_img()
