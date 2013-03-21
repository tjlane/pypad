#!/usr/bin/python

# Usage:
# In this directory, type:
#    ./viewCleanedSum.py -rxxxx -a (optional) -m 10000
# For details, type 
#	 python viewCleanedSum.py --help
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
parser.add_option("-a", "--avg", action="store_true", dest="ang_average", 
					help="flag to compute angular average", default=False)
parser.add_option("-m", "--max", action="store", type="float", dest="max_value", 
					help="Mask out pixels above this value", metavar="max_value", default="10000000")
					

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
source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
write_dir = ""

runtag = "r%s"%(options.runNumber)
print source_dir+runtag+"/"+runtag+"_200imtest_cleaned_average.h5"
f = H.File(source_dir+runtag+"/"+runtag+"_200imtest_cleaned_average.h5","r")
d = N.array(f['/data/data'])
d *= (d < options.max_value)

f.close()

vx = N.arange(-512,512)
vy = N.arange(-512,512)
X,Y = N.meshgrid(vx,vx)
arr = (N.sqrt(X*X + Y*Y)).astype(int)

lenOfAvg = len(set(arr.flatten()))+5
avg_count = N.zeros(lenOfAvg)
avg_vals = N.zeros(lenOfAvg)
avg_var = N.zeros(lenOfAvg)

intensMask = (d.flatten()>0.)

#blockMask=N.zeros((1758,1758))
#blockMask[:950,950:1300]=1.
#blockMask[:750,700:950]=1.
#blockMask[:375,375:700]=1.
##blockMask[1180:,600:970]=1.
##blockMask[1036:,855:1252]=1.
#intensMask *= blockMask.flatten()


########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, filename):
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		self.cmax = 0.1*self.inarr.max()
		self.cmin = self.inarr.min()
	
	def on_keypress(self,event):
		if event.key == 'p':
			if not os.path.exists(write_dir + runtag):
				os.mkdir(write_dir + runtag)
			pngtag = write_dir + runtag + "/%s.png" % (self.filename)	
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
		self.axes = P.imshow(self.inarr, vmax = self.cmax)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		P.show() 

print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

if (options.ang_average is True):
	compressedIntens = N.compress(intensMask, d.flatten())
	compressedPos = N.compress(intensMask, arr.flatten())

	for i,j in zip(compressedPos,compressedIntens):
		avg_count[i] += 1.
		avg_vals[i] += j
		avg_var[i] += j*j

	for i in range(len(avg_count)):
		if(avg_count[i] > 0.):
			avg_vals[i] /= avg_count[i]
			avg_var[i] /= avg_count[i]

	P.figure()
	P.plot(avg_vals)
	P.draw()
	#avg_vals.tofile(write_dir + runtag + "_ang_avg.dat", sep="\n",format="%e")

currImg = img_class(d, runtag + "_cleaned200")
currImg.draw_img()

P.show()
