#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# Two first arguments: input HDF5 files containing pixels masks in /data/data
# Third argument: output HDF5 file containing merged pixel mask in /data/data

import os
import sys
import string
import numpy as N
import h5py as H
import matplotlib
import matplotlib.pyplot as P

inputfile1 = sys.argv[1]
inputfile2 = sys.argv[2]
outputfile = sys.argv[3]

ROWS = 194
COLS = 185
RAW_DATA_LENGTH = [8*COLS,8*ROWS]

# Open masks to merge
f1 = H.File(inputfile1,"r")
mask1 = N.array(f1['/data/data'])
f1.close()
masked1 = N.where(mask1 < 1)
print "%s contains %s masked pixels." % (inputfile1, len(mask1[masked1]))

f2 = H.File(inputfile2,"r")
mask2 = N.array(f2['/data/data'])
f2.close()
masked2 = N.where(mask2 < 1)
print "%s contains %s masked pixels." % (inputfile2, len(mask2[masked2]))

# Create output array of type 16-bit integer
outputmask = mask1*mask2
maskedo = N.where(outputmask < 1)

# Save output array to HDF5 file
h5output = H.File(outputfile,'w')
datagroup = h5output.create_group("data")
dataset = datagroup.create_dataset("data",RAW_DATA_LENGTH,dtype="int16")
dataset[...] = outputmask[:,:]
h5output.close()
print "Saved new mask as %s containing %s masked pixels." % (outputfile, len(outputmask[maskedo]))

#########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
#########################################################
class img_class (object):
	def __init__(self, inarr, filename):
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		self.cmax = self.inarr.max()
		self.cmin = self.inarr.min()
	
	def on_keypress(self,event):
		if event.key == 'p':
			pngtag = "%s.png" % (self.filename)	
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
		#P.show() 
		pngtag = "%s.png" % (self.filename)	
		print "saving image as " + pngtag 
		P.savefig(pngtag)
		P.close()

print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

currImg = img_class(outputmask, outputfile)
currImg.draw_img()
