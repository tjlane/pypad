#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./viewIntensityHistograms.py -rxxxx -m 400 (optional) -x 800 (optional)
# For details, type 
#	 python viewIntensityHistograms.py --help
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
parser.add_option("-m", "--min", action="store", type="float", dest="min_value",
                                        help="ignore intensities below this value", metavar="min_value", default="-1")
parser.add_option("-x", "--max", action="store", type="float", dest="max_value",
                                        help="ignore intensities above this value", metavar="max_value", default="100")

(options, args) = parser.parse_args()

import numpy as N
import h5py as H

import matplotlib
import matplotlib.pyplot as P

import peakdetect as PD

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
searchstring="[a-zA-Z0-9\_]+"+runtag+"[a-z0-9\_]+-hist.h5"
h5pattern = re.compile(searchstring)
h5files = [h5pattern.findall(x) for x in os.listdir(source_dir+runtag)]
h5files = [items for sublists in h5files for items in sublists]

print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

########################################################
# Loop to display all H5 files found.
########################################################

#In [21]: len(d)
#Out[21]: 32
#
#In [22]: len(d[0])
#Out[22]: 131071

for fname in h5files:
	f = H.File(source_dir+runtag+"/"+fname, 'r')
	d = N.array(f['/data/data'])
	nasic = len(d)
	nhist = len(d[0])
	adu = N.arange(nhist)-65535
	f.close()
	
	for asics in N.arange(nasic):
		[max, min] = PD.peakdet(d[asics][65535-100:65535+2000], 20)
		xmax = [max[i][0]-100 for i in N.arange(len(max))]
		ymax = [max[i][1] for i in N.arange(len(max))]
		xmin = [min[i][0]-100 for i in N.arange(len(min))]
		ymin = [min[i][1] for i in N.arange(len(min))]
		#print xmax, ymax, xmin, ymin
		fig = P.figure()
		P.plot(adu, d[asics])
		P.scatter(xmax, ymax, c='g', marker='o')
		P.scatter(xmin, ymin, c='r', marker='o')
		canvas = fig.add_subplot(111)
		canvas.set_title(fname + "_asic%s"%(asics+1))
		P.xlim(options.min_value, options.max_value)
		P.xlabel("Intensity (ADUs)")
		P.ylabel("Number of Pixels")
		#P.draw()
		P.show()
	
	#P.show()


