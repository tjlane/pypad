#!/usr/bin/python

# Usage:
# First argument: input HDF5 file containing pixel mask in /data/data
# Second argument: output HDF5 file containing pixel mask with masked borders in /data/data

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

(options, args) = parser.parse_args()

write_dir = os.getcwd()

ROWS = 194
COLS = 185
RAW_DATA_LENGTH = [8*COLS,8*ROWS]

if (options.maskName == ""):
        outputfilename = "pixel_mask_%s" % (options.fileName)
else:
	outputfilename = options.maskName

# Open input auto mask
f = H.File(options.fileName,"r")
mask = N.array(f['/masks/auto'])
f.close()

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

maskImg = img_class(mask,options.fileName)
maskImg.draw_img()

# Create output array of type 16-bit integer with increased automask
outputmask = N.ones(RAW_DATA_LENGTH, dtype="int16")
maskindices = N.where(mask > 0)
outputmask[maskindices] = 0
counter = 0
for i in range(len(mask[maskindices])):
	row = maskindices[0][i]
	col = maskindices[1][i]
	# make an exception for row13, it does not have to be increased in size
	skip = 0
	for n in range(8):
		if (col == n*ROWS-13):
			skip = 1
	if (skip == 1):
		continue
	if ((row+1) % COLS != 0):
		if (outputmask[row+1,col] == 1):
			counter += 1
		outputmask[row+1,col] = 0
		if ((col+1) % ROWS != 0):
			if (outputmask[row+1,col+1] == 1):
				counter += 1
			outputmask[row+1,col+1] = 0
		if (((col-1)-(ROWS-1)) % ROWS != 0):
			if (outputmask[row+1,col-1] == 1):
				counter += 1
			outputmask[row+1,col-1] = 0
	if (((row-1)-(COLS-1)) % COLS != 0):
		if (outputmask[row-1,col] == 1):
			counter += 1
		outputmask[row-1,col] = 0
		if ((col+1) % ROWS != 0):
			if (outputmask[row-1,col+1] == 1):
				counter += 1
			outputmask[row-1,col+1] = 0
		if (((col-1)-(ROWS-1)) % ROWS != 0):
			if (outputmask[row-1,col-1] == 1):
				counter += 1
			outputmask[row-1,col-1] = 0
	if ((col+1) % ROWS != 0):
		if (outputmask[row,col+1] == 1):
			counter += 1
		outputmask[row,col+1] = 0
	if (((col-1)-(ROWS-1)) % ROWS != 0):
		if (outputmask[row,col-1] == 1):
			counter += 1
		outputmask[row,col-1] = 0

print "Input mask contained %s masked pixels, masked out an additional %s pixels." % (len(mask[maskindices]), counter)

# Save output array to HDF5 file
h5output = H.File(outputfilename,'w')
datagroup = h5output.create_group("data")
dataset = datagroup.create_dataset("data",RAW_DATA_LENGTH,dtype="int16")
dataset[...] = outputmask[:,:]
h5output.close()

outputmaskImg = img_class(outputmask, outputfilename)
outputmaskImg.draw_img()
