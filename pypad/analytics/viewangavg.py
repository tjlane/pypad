#!/reg/g/psdm/sw/releases/ana-current/arch/x86_64-rhel5-gcc41-opt/bin/python

# Usage:
# In this directory, type:
#    ./viewangavg.py -rxxxx -m 400 (optional) -x 800 (optional)
# For details, type 
#	 python viewangavg.py --help
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
					help="ignore intensities below this q-value", metavar="min_value", default="0")
parser.add_option("-x", "--max", action="store", type="float", dest="max_value", 
					help="ignore intensities above this q-value", default="100000")
parser.add_option("-t", "--tag", action="store", type="string", dest="fileTag",
		                        help="optional file tag for run", metavar="FILETAG", default="")

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
#source_dir = "/reg/d/psdm/cxi/cxi74613/scratch/cleaned_hdf5/"
source_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/test_runs/"
#source_dir = "/reg/d/psdm/cxi/cxi74613/ftc/cleaned_hdf5/"
#write_dir = "/reg/d/psdm/cxi/cxi74613/ftc/sellberg/figures/"
write_dir = "/reg/neh/home3/sellberg/NML-2013/analysis/cheetah_scripts/figures/"


runtag = "r%s"%(options.runNumber)
print source_dir+runtag+"/"+runtag+"-angavg.h5"
f = H.File(source_dir+runtag+"/"+runtag+"-angavg.h5","r")
d = N.array(f['/data/data'])
q = d[0]
saxs = d[1]
f.close()

print "saving data as "+source_dir+runtag+"/"+runtag+"-angavg.dat"
d.tofile(source_dir+runtag+"/"+runtag+"-angavg.dat", sep="\n")

#for i in N.arange(len(saxs)):
#    print "Q: %f, I: %f"%(q[i], saxs[i])

print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

def on_keypress(event):
    if event.key == 'p':
        if not os.path.exists(write_dir + runtag):
            os.mkdir(write_dir + runtag)
        if (options.fileTag == ""):
            pngtag = write_dir + runtag + "/%s_angavg.png" % (runtag)
        else:
            pngtag = write_dir + runtag + "/%s-%s_angavg.png" % (runtag, options.fileTag)
        print "saving image as " + pngtag
        P.savefig(pngtag)

fig = P.figure()
cid1 = fig.canvas.mpl_connect('key_press_event', on_keypress)
P.plot(q, saxs)
canvas = fig.add_subplot(111)
canvas.set_title(runtag + "_angavg")
if N.amax(q) < 10:
    P.xlabel("Q (A-1)")
else:
    P.xlabel("Q (pixels)")
P.ylabel("I(Q)")
P.show()
#if (options.fileTag == ""):
#    pngtag = write_dir + runtag + "/%s_angavg.png" % (runtag)
#else:
#    pngtag = write_dir + runtag + "/%s-%s_angavg.png" % (runtag, options.fileTag)
#print "saving image as " + pngtag
#P.savefig(pngtag)
#P.close()
