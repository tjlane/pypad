#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./viewDarkcalVariance.py -rxxxx
# For details, type 
#	 python viewDarkcalVariance.py --help
# where rxxxx is the run number of hits and nonhits found using the hitfinder executable. 
# By default, this script looks into the h5 files that are in the appropriate rxxxx directory
#

import os
import sys
import string
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", 
					help="run number you wish to view", metavar="rxxxx", default="")

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
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
#source_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/test_runs/"
#write_dir = "/reg/d/psdm/cxi/cxi74613/ftc/sellberg/figures/"
write_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/figures/"

runtag = "r%s"%(options.runNumber)
print source_dir+runtag+"/"+runtag+"-darkcal_variance.h5"
f = H.File(source_dir+runtag+"/"+runtag+"-darkcal_variance.h5","r")
d = N.array(f['/data/data'])
f.close()

print "DARKCAL STATISTICS"
print "min = %f, max = %f, mean = %f, median = %f, std = %f" % (N.min(d), N.max(d), N.mean(d), N.median(d), N.std(d))
#negindex = N.where(d < 0)
#print len(negindex), len(d[negindex])
posindex = N.where(d > 64)
print len(d[posindex]), "pixels with variance >64"
posindex = N.where(d > 49)
print len(d[posindex]), "pixels with variance >49"
#posindex = N.where(d > 40)
#print len(d[posindex]), "pixels with variance >40"
posindex = N.where(d > 36)
print len(d[posindex]), "pixels with variance >36"
posindex = N.where(d > 25)
print len(d[posindex]), "pixels with variance >25"
posindex = N.where(d > 16)
print len(d[posindex]), "pixels with variance >16"
posindex = N.where(d > 9)
print len(d[posindex]), "pixels with variance >9"
#print len(d), len(d[0])

#hist_bins = N.arange(N.floor(d.min()), N.ceil(d.max()) + 2, 0.1) - 0.05
#hist_bins = N.arange(N.floor(d.min()), N.ceil(d.max()) + 2) - 0.5
hist_bins = N.arange(0, 50 + 2) - 0.5
hist, hist_bins = N.histogram(d, bins=hist_bins)
hist_bins = [(hist_bins[j] + hist_bins[j+1])/2 for j in range(len(hist))]

fig = P.figure()
canvas = fig.add_subplot(111)
#canvas.set_title("Histogram %s, min=%d, max=%d, mean=%.1f, std=%.1f" % (runtag, d.min(), d.max(), d.mean(), d.std()))
canvas.set_title("Histogram %s" % (runtag))
P.bar(hist_bins, hist, align='center')
P.xlabel("Intensity (ADUs)")
P.ylabel("Number of pixels")
P.show()


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
			pngtag = write_dir + runtag + "/%s-darkcal_variance.png" % (self.filename)	
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
		self.axes = P.imshow(self.inarr, vmin = 0, vmax = 50)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		P.show()

print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

currImg = img_class(d, runtag)
currImg.draw_img()

P.show()
