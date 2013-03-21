#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./viewPixelStatistics.py -rxxxx
# For details, type
#        python viewPixelStatistics.py --help
# where rxxxx is the run number of hits and nonhits found using the hitfinder executable. 
# By default, this script looks into the h5 files that are in the appropriate rxxxx directory
#

import numpy as N
import h5py as H
import matplotlib
import matplotlib.pyplot as P
import sys
import os
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="run number you wish to view", metavar="rxxxx", default="")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="prints out the frame number as it is processed", default=False)

(options, args) = parser.parse_args()


########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
#source_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/test_runs/"
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
#write_dir = "/reg/d/psdm/cxi/cxi74613/ftc/sellberg/figures/"
write_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/figures/"

if(options.runNumber is not ""):
	print "Now examining H5 files in %sr%s/ ..."%(source_dir,options.runNumber)

runtag = "r%s"%(options.runNumber)

########################################################
# Search specified directory for *.h5 files
########################################################
searchstring="[a-zA-Z0-9\_]+"+runtag+"[a-z0-9\_]+-pixels.h5"
h5pattern = re.compile(searchstring)
h5files = [h5pattern.findall(x) for x in os.listdir(source_dir+runtag)]
h5files = [items for sublists in h5files for items in sublists]
numFiles = len(h5files)

print "found %d files." % (numFiles)

print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

########################################################
# Loop to display all H5 files found. 
########################################################

pixelIntens = []
fcounter = 0
for fname in h5files:
	f = H.File(source_dir+runtag+"/"+fname, 'r')
	d = N.array(f['/data/data'])
	f.close()
	fcounter += 1
	if (options.verbose and (round(((fcounter*10)%numFiles)/10)==0)):
		print str(fcounter) + " of " + str(numFiles) + " files read (" + str(fcounter*100/numFiles) + "%)"
	pixelIntens.append(d[0])
	
# PLOT DATA
pixelIntens = N.array(pixelIntens).transpose()
pixelHist = []
pixelHistBins = []
for i in range(len(pixelIntens)):
	hist_bins = N.arange(N.floor(pixelIntens[i].min()), N.ceil(pixelIntens[i].max()) + 2) - 0.5
	
	#nbins = int(N.ceil(pixelIntens[i].max() - pixelIntens[i].min() + 1))
	hist, hist_bins = N.histogram(pixelIntens[i], bins=hist_bins)
	hist_bins = [(hist_bins[j] + hist_bins[j+1])/2 for j in range(len(hist))]
	pixelHist.append(hist)
	pixelHistBins.append(hist_bins)
	
	fig = P.figure()
	canvas = fig.add_subplot(111)
	canvas.set_title("Pixel #%d, min=%d, max=%d, mean=%.1f, std=%.1f" % ((i+1), pixelIntens[i].min(), pixelIntens[i].max(), pixelIntens[i].mean(), pixelIntens[i].std()))
	P.bar(hist_bins, hist, align='center')
	P.xlabel("Intensity (ADUs)")
	P.ylabel("Histogram")
	
	pngtag = write_dir + runtag + "/pixel%d-%.2fADUs.png" % ((i+1), pixelIntens[i].std()*pixelIntens[i].std())
	P.savefig(pngtag)
	print "%s saved."%(pngtag)
	P.show()

