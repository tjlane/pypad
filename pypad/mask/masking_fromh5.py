#!/usr/bin/env python
import Make_mask as m
import h5py as h
import pylab as p
import numpy as n

def mask_fromh5(filename,mask_file,threshold=0):
   file=h.File(filename,'r')
   data=n.array(file.get('/data/data'))
#   p.imshow(data)
#   p.show()
#   file.close()
   m.Make_mask(data,mask_file,threshold) 



filename='/reg/d/psdm/cxi/cxi35711/scratch/bad_pixels_masking/r0007-RawSum.h5'
mask_file='/reg/d/psdm/cxi/cxi35711/scratch/bad_pixels_masking/mask_r0007.h5'
threshold=20
mask_fromh5(filename,mask_file,threshold)



