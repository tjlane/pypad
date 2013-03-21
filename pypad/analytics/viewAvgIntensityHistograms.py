#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./viewAvgIntensityHitsograms.py -rxxxx
# For details, type 
#	 python viewAvgIntensityHitsograms.py --help
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
#source_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/test_runs/"
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
#write_dir = "/reg/d/psdm/cxi/cxi74613/ftc/sellberg/figures/"
write_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/figures/"

runtag = "r%s"%(options.runNumber)
print source_dir+runtag+"/"+runtag+"-intensity_histograms.h5"
f = H.File(source_dir+runtag+"/"+runtag+"-intensity_histograms.h5","r")
d = N.array(f['/data/data'])
intensities = d[0]
hitintensities = d[1]
xi = d[2]
f.close()

#intensities[N.isnan(intensities)] = 0
#hitintensities[N.isnan(hitintensities)] = 0

print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

fig = P.figure()
canvas = fig.add_subplot(111)
P.plot(xi, intensities, color='b', label='all events')
P.plot(xi, hitintensities, color='r', label='only hits')
canvas.set_title(runtag + "_intensities")
handles, labels = canvas.get_legend_handles_labels()
canvas.legend(handles, labels, loc='upper right')
P.xlabel("Iavg (ADU)")
P.ylabel("events")
P.show()
