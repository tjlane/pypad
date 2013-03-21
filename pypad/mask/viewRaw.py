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
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", action="store", type="string", dest="fileName", 
					help="input file you wish to view", metavar="FILENAME", default="")
parser.add_option("-m", "--max", action="store", type="float", dest="max",
		                        help="upper viewing limit", metavar="MAX_VALUE", default="0")
parser.add_option("-n", "--min", action="store", type="float", dest="min",
		                        help="lower viewing limit", metavar="MIN_VALUE", default="0")

(options, args) = parser.parse_args()

import numpy as N
import h5py as H

import matplotlib
import matplotlib.pyplot as P

########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
source_dir = os.getcwd()+'/'
write_dir = source_dir

filename = options.fileName
print source_dir+filename
f = H.File(source_dir+filename,"r")
d = N.array(f['/data/data'])
f.close()

########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, filename):
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		if (options.max == options.min):
			self.cmax = 0.1*self.inarr.max()
			self.cmin = self.inarr.min()
		else:
			self.cmax = options.max
			self.cmin = options.min
	
	def on_keypress(self,event):
		if event.key == 'p':
			if not os.path.exists(write_dir):
				os.mkdir(write_dir)
			pngtag = write_dir + "/%s.png" % (self.filename)	
			print "saving image as " + pngtag 
			P.savefig(pngtag)
		if event.key == 'r':
			self.cmax = self.inarr.max()
			self.cmin = self.inarr.min()
			P.clim(self.cmin, self.cmax)
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
		#self.axes = P.imshow(self.inarr, vmin = self.cmin, vmax = self.cmax)
		self.axes = P.imshow(self.inarr, vmin = 0, vmax = 2)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		P.show()

print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

currImg = img_class(d, filename)
currImg.draw_img()
