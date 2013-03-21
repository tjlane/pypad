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
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", 
					help="run number you wish to view", metavar="rXXXX", default="")
parser.add_option("-o", "--output", action="store", type="string", dest="maskName",
		                        help="output mask you wish to save (default = 'pixel_mask_rXXXX_variance.h5')", metavar="MASKNAME", default="")
parser.add_option("-t", "--threshold", action="store", type="float", dest="threshold",
		                        help="threshold (ADU^2) above which masking is performed (default = 16)", metavar="THRESHOLD", default="16")
parser.add_option("-d", "--dead_thresh", action="store", type="float", dest="dead",
		                        help="threshold (ADU^2) below which pixel is considered dead (default = 0)", metavar="THRESHOLD", default="0")

(options, args) = parser.parse_args()

########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
#source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/darkcal/"
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/darkcal/"
#write_dir = "/reg/neh/home/sellberg/NML-2013/analysis/cheetah_scripts/figures/"
write_dir = "/reg/neh/home/sellberg/NML-2013/analysis/masking_script/"

ROWS = 194
COLS = 185
RAW_DATA_LENGTH = [8*COLS,8*ROWS]

runtag = "r%s"%(options.runNumber)
if (options.maskName == ""):
	outputfilename = "pixel_mask_cxi74613_%s_variance_%s-%sADUs.h5" % (runtag, str(options.dead), str(options.threshold))
else:
	outputfilename = options.maskName

print source_dir+runtag+"-darkcal_variance.h5"
f = H.File(source_dir+runtag+"-darkcal_variance.h5","r")
d = N.array(f['/data/data'])
f.close()

# Find pixels to mask
negindex = N.where(d < 0)
if (len(d[negindex]) > 0):
	print "WARNING: the indices below have negative variance, double check the algorithm!"
	print negindex
deadindex = N.where(d <= options.dead)
if (len(d[deadindex]) > 0):
	print len(d[deadindex])-len(d[negindex]), "pixels are dead (variance <= %s ADU)" % (str(options.dead))
posindex = N.where(d > options.threshold)

# Create output array of type 16-bit integer with ones, fill in zeros
outputmask = N.ones(RAW_DATA_LENGTH, dtype="int16")
outputmask[deadindex] = 0
outputmask[posindex] = 0
print "Masked out " + str(len(d[deadindex]) + len(d[posindex])) + " pixels."

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
			if not os.path.exists(write_dir + runtag):
				os.mkdir(write_dir + runtag)
			pngtag = write_dir + "pixel_mask_cxi74613_%s_variance_%s-%sADUs.png" % (self.filename, str(options.dead), str(options.threshold))
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
		canvas.set_title('pixel_mask_'+self.filename+'_variance')
		P.rc('image',origin='lower')
		self.axes = P.imshow(self.inarr, vmin = 0, vmax = options.threshold)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		#P.show()
		pngtag = write_dir + "pixel_mask_cxi74613_%s_variance_%s-%sADUs.png" % (self.filename, str(options.dead), str(options.threshold))
		print "saving image as " + pngtag 
		P.savefig(pngtag)
		P.close()
			

print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

currImg = img_class(d*outputmask, runtag)
currImg.draw_img()
