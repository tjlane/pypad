#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./analyzeDarkcals.py -f FILENAME (optional)
# For details, type 
#	 python analyzeDarkcals.py --help
# where FILENAME is the name of the filelist where all the darkcals are listed
# By default, this script looks into the hdf5 files that are in the appropriate darkcal directory
#

import numpy as N
import h5py as H
import glob as G
import matplotlib
import matplotlib.pyplot as P
from pylab import *
import scipy
import scipy.interpolate as I
from scipy import *
from scipy import optimize
import sys, os, re, shutil, subprocess, time
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", action="store", type="string", dest="filelist",
		  help="input filelist you wish to analyze (default: filelist.dat)", metavar="FILENAME", default="filelist.dat")
parser.add_option("-t", "--tag", action="store", type="string", dest="filetag",
		                        help="optional file tag for saved files", metavar="FILETAG", default="")
parser.add_option("-n", "--noise_thresh", action="store", type="float", dest="threshold",
		  help="variance threshold (ADU^2) above which pixel is considered noisy (default = 16)", metavar="THRESHOLD", default="16")
parser.add_option("-d", "--dead_thresh", action="store", type="float", dest="dead",
		  help="variance threshold (ADU^2) below which pixel is considered dead (default = 0)", metavar="THRESHOLD", default="0")

(options, args) = parser.parse_args()

########################################################
# Edit this variable accordingly
# Files are read for source_dir and
# written to write_dir.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/darkcal/"
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/darkcal/"
#write_dir = "/reg/d/psdm/cxi/cxi74613/ftc/sellberg/figures/"
write_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/figures/"

runs = []
runnumber = []
pedmean = []
pedmedian = []
pedstd = []
varmean = []
varmedian = []
varstd = []
varmasked = []

print "Reading darkcals from %s ..." % (source_dir+options.filelist)
f = open(source_dir+options.filelist, 'r')
for line in f:
	if line[0] != '#':
		line = re.sub("\n", '', line)
		run = int(re.sub("r", '', line))
		runs.append(line)
		runnumber.append(run)

f.close()

for runtag in runs:
	
	print source_dir+runtag+"-darkcal.h5"
	f = H.File(source_dir+runtag+"-darkcal.h5","r")
	d = N.array(f['/data/data'])
	pedmean.append(d.mean())
	pedmedian.append(N.median(d))
	pedstd.append(d.std())
	f.close()
	
	print source_dir+runtag+"-darkcal_variance.h5"	
	f = H.File(source_dir+runtag+"-darkcal_variance.h5","r")
	d = N.array(f['/data/data'])
	varmean.append(d.mean())
	varmedian.append(N.median(d))
	varstd.append(d.std())
	deadindex = N.where(d <= options.dead)
	posindex = N.where(d > options.threshold)
	varmasked.append(len(d[deadindex]) + len(d[posindex]))
	f.close()

print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
#print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

#def on_keypress(event):
#    if event.key == 'p':

#cid1 = fig.canvas.mpl_connect('key_press_event', on_keypress)

fig = P.figure(num=None, figsize=(13.5, 11), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(221)
canvas.set_title("Pedestal mean", fontsize='medium')
P.plot(runnumber, pedmean)
P.xlabel("run")
P.ylabel("mean intensity (ADU/pixel)")

canvas = fig.add_subplot(222)
canvas.set_title("Pedestal median", fontsize='medium')
P.plot(runnumber, pedmedian)
P.xlabel("run")
P.ylabel("median intensity (ADU/pixel)")

canvas = fig.add_subplot(223)
canvas.set_title("Pedestal standard deviation", fontsize='medium')
P.plot(runnumber, pedstd)
P.xlabel("run")
P.ylabel("pixel std (ADU)")

if (options.filetag == ""):
	pngtag = write_dir + "darkcals_pedestal.png"
else:
	pngtag = write_dir + "darkcals-%s_pedestal.png" % (options.filetag)
print "saving image as " + pngtag
P.savefig(pngtag)

fig = P.figure(num=None, figsize=(13.5, 11), dpi=100, facecolor='w', edgecolor='k')
canvas = fig.add_subplot(221)
canvas.set_title("Variance mean", fontsize='medium')
P.plot(runnumber, varmean)
P.xlabel("run")
P.ylabel("mean variance (ADU^2/pixel)")

canvas = fig.add_subplot(222)
canvas.set_title("Variance Median", fontsize='medium')
P.plot(runnumber, varmedian)
P.xlabel("run")
P.ylabel("median variance (ADU^2/pixel)")

canvas = fig.add_subplot(223)
canvas.set_title("Variance standard deviation", fontsize='medium')
P.plot(runnumber, varstd)
P.xlabel("run")
P.ylabel("pixel std (ADU^2)")

canvas = fig.add_subplot(224)
canvas.set_title("Pixels with variance outside %s-%s ADU^2" % (options.dead, options.threshold), fontsize='medium')
P.plot(runnumber, varmasked)
P.xlabel("run")
P.ylabel("masked pixels")

if (options.filetag == ""):
	pngtag = write_dir + "darkcals_variance.png"
else:
	pngtag = write_dir + "darkcals-%s_variance.png" % (options.filetag)
print "saving image as " + pngtag
P.savefig(pngtag)

P.show()

#print "DARKCAL STATISTICS"
#print "min = %f, max = %f, mean = %f, median = %f, std = %f" % (N.min(d), N.max(d), N.mean(d), N.median(d), N.std(d))
#
##hist_bins = N.arange(N.floor(d.min()), N.ceil(d.max()) + 2, 0.1) - 0.05
##hist_bins = N.arange(N.floor(d.min()), N.ceil(d.max()) + 2) - 0.5
#hist_bins = N.arange(0, 50 + 2) - 0.5
#hist, hist_bins = N.histogram(d, bins=hist_bins)
#hist_bins = [(hist_bins[j] + hist_bins[j+1])/2 for j in range(len(hist))]

#fig = P.figure()
#canvas = fig.add_subplot(111)
##canvas.set_title("Histogram %s, min=%d, max=%d, mean=%.1f, std=%.1f" % (runtag, d.min(), d.max(), d.mean(), d.std()))
#canvas.set_title("Histogram %s" % (runtag))
#P.bar(hist_bins, hist, align='center')
#P.xlabel("Intensity (ADUs)")
#P.ylabel("Number of pixels")
#P.show()


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


#currImg = img_class(d, runtag)
#currImg.draw_img()

#P.show()

#ROWS = 194
#COLS = 185
#RAW_DATA_LENGTH = [8*COLS,8*ROWS]
#
#if (options.maskName == ""):
#	outputfilename = "pixel_mask_cxi74613_%s_variance_%s-%sADUs.h5" % (runtag, str(options.dead), str(options.threshold))
#else:
#	outputfilename = options.maskName
#
#print source_dir+runtag+"-darkcal_variance.h5"
#f = H.File(source_dir+runtag+"-darkcal_variance.h5","r")
#d = N.array(f['/data/data'])
#f.close()
#
## Create output array of type 16-bit integer with ones, fill in zeros
#outputmask = N.ones(RAW_DATA_LENGTH, dtype="int16")
#outputmask[deadindex] = 0
#outputmask[posindex] = 0
#print "Masked out " + str(len(d[deadindex]) + len(d[posindex])) + " pixels."
#
## Save output array to HDF5 file
#h5output = H.File(outputfilename, 'w')
#datagroup = h5output.create_group("data")
#dataset = datagroup.create_dataset("data",RAW_DATA_LENGTH,dtype="int16")
#dataset[...] = outputmask[:,:]
#h5output.close()
#print "Saved mask as %s" % (outputfilename)

#print "saving data as "+source_dir+runtag+"/"+runtag+"-angavg.dat"
#d.tofile(source_dir+runtag+"/"+runtag+"-angavg.dat", sep="\n")

