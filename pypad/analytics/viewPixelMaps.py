#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./viewPixelMaps.py -rxxxx
# For details, type
#        python viewPixelMaps.py --help
# where rxxxx is the run number of hits the produced pixel maps using the hitfinder executable. 
# By default, this script looks into the h5 files that are in the appropriate rxxxx directory
#

import os
import sys
import string
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", 
					help="run number you wish to view", metavar="XXXX", default="")
parser.add_option("-t", "--tag", action="store", type="string", dest="fileTag",
		  			help="file tag for run (default: img)", metavar="FILETAG", default="img")

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

########################################################
# Search specified directory for *.h5 files
########################################################
searchstring = "pix"+"[A-Z]+_+[a-z]+.h5"
h5pattern = re.compile(searchstring)
h5files = [h5pattern.findall(x) for x in os.listdir(source_dir+runtag)]
h5files = [items for sublists in h5files for items in sublists]

searchstring_asm = "pix"+"[A-Z]+_asm.h5"
h5pattern_asm = re.compile(searchstring_asm)
h5files_asm = [h5pattern_asm.findall(x) for x in os.listdir(source_dir+runtag)]
h5files_asm = [items for sublists in h5files_asm for items in sublists]
sh5files_asm = set(h5files_asm)

colmax = 100000
colmin = -100000

########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, filename, asm):
		self.inarr = inarr
		if (asm):
			for i in range(len(inarr)):
				self.inarr[i] = self.inarr[i][::-1]
		self.filename = filename
		self.cmax = self.inarr.max()
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
				

	def save_img(self):
		fig = P.figure()
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin = 'lower', vmin = self.cmin, vmax = self.cmax)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		if not os.path.exists(write_dir + runtag):
			os.mkdir(write_dir + runtag)
		pngtag = write_dir + runtag + "/%s.png" %(self.filename)
		print "saving image as " + pngtag
		P.savefig(pngtag)
		P.close()

print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

########################################################
# Loop to display all H5 files found.
########################################################
for fname in h5files:
	print source_dir+runtag+"/"+fname
	f = H.File(source_dir+runtag+"/"+fname, 'r')
	d = N.array(f['/data/data'])
	f.close()
	currImg = img_class(d, fname, fname in sh5files_asm)
	currImg.save_img()

#P.show()
