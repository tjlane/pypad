
# Code to find the water streak in cspad images. Create a mask

import os
import numpy as N
import scipy as S
import scipy.ndimage as Sndim
import h5py
from matplotlib.pyplot import imsave as pltimsave
import matplotlib.image as mpimg
import my_favorites as mf
import time

# This is the main program.
def find_streak():
    
    # Define paths and names. Make the list of files to process
    #paths for SLAC
    datapath = "/reg/neh/home/sellberg/NML-2011/analysis/cheetah_scripts/test_runs/r0125/"
    outpath = "/reg/neh/home/sellberg/NML-2011/analysis/cheetah_scripts/test_runs/streak_cleaned/r0125/"    

    #paths for local windows
    #datapath = "F:\\Work back up\\LCLSWater2011\\LCLS_2011_feb_water\\r0101\\"
    #outpath = "F:\\Work back up\\LCLSWater2011\\LCLS_2011_feb_water\\"
    description = "r0116_100imtest"
    #outname = "r0107_test_streaks_removed.h5"
    maskname = description+"_mask"
    h5sumname = description+"_h5sum"
    imagesumname = description+"_imagesum"
    cleanedaveragename = description+"_cleaned_average"

    flist = MakeFileList(datapath,".h5")
    flen = len(flist)
    print "flen =", flen

    if flist[0] == -1 : 
        print "No files found"
        return
    
    # read in the first file
    # this reads the size of the file
    h5file = h5py.File(flist[0],"r")
    h5data = h5file["data/data"]    
    image = N.empty(shape=h5data.shape,dtype=h5data.dtype)
    image = h5data[...]
    test = size(image)
    #implot = plt.imshow(image)
    #plt.draw()
    s_init = h5data.shape
    xshift_init = 0
    yshift_init = 0
    cenx_init = s_init[0]/2 + xshift_init
    ceny_init = s_init[1]/2 + yshift_init
    h5file.close()
       
    # crop the image to 1024x1024 to make it faster.
    newsize = 1024
    s = [newsize,newsize]
    cenx = 512
    ceny = 512
    xshift = 0
    yshift = 0

    #define the centre region
    radius = 60
    circ_centre = mf.circle(s[0],s[1],rad=radius,cenx=cenx,ceny=ceny) 
    circ_lower = N.zeros([s[0],s[1]])
    circ_lower[:,0:s[1]/2] = circ_centre[:,0:s[1]/2]
    circ_upper = N.zeros([s[0],s[1]])
    circ_upper[:,s[1]/2:s[1]] = circ_centre[:,s[1]/2:s[1]]


    # set the threshold and the morphology masks
    threshold = 0.01
    radius = 150
    hsize = 5
    #circ = mf.circle(s[0],s[1],rad=radius,cenx=cenx,ceny=ceny) 
    mostruct = N.ones([3,3],dtype="bool")
    mcstruct = mostruct
    mdstruct = N.ones([15,15],dtype="bool")

    h5sum = N.zeros(s)
    masksum = N.zeros(s,dtype="float")
    imagesum = N.zeros(s)

    # initialise display
    #init_display()
    
    # Loop over the file list
    count = 0
    #flist2 = [flist[2]] 
    max = 100
    for ii in N.arange(max):
        i = flist[ii]
	count += 1
	print count, "/", max
        print i

        #read the next image
        h5file = h5py.File(i,"r")
        h5data = h5file["data/data"]    
        image = h5data[...]
        image = crop_array(image,newsize,newsize,cenx_init,ceny_init)
        h5file.close()

	print image
        #define the mask using the threshold
        print N.amax(image) 
	istreak = N.where(image > N.amax(image)*threshold)   
        mask = N.zeros(image.shape,dtype="bool")
	mask[istreak] = 1
        
        
        # apply morphological operations
        mask = Sndim.morphology.binary_closing(mask,structure=mcstruct)
        mask = Sndim.morphology.binary_opening(mask,structure=mostruct)

        # compensate for missing hole
        mask[cenx-hsize:cenx+hsize,ceny-hsize:ceny+hsize] = 1.

        # dilate the mask by itself  	
        mask = dilate_fft(mask)
        mask = Sndim.morphology.binary_dilation(mask,structure=mdstruct)
        print "Done with dilation"

        # find the region at the centre
        masklabelled, num_features = Sndim.label(mask)
        print "Number of regions: ", num_features
        iregion = 0
	iregion2 = 0
        for i in N.arange(num_features):
            if i == 0: continue
            ilabel = N.where(masklabelled == i)
            print i, len(ilabel[0]) 
            mask2 = N.zeros(mask.shape)
            mask2[ilabel] = 1.
            if N.sum(mask2*circ_lower) != 0: iregion = i
            if N.sum(mask2*circ_upper) != 0: iregion2 = i

        print "iregions ", iregion, iregion2
        mask2 = N.zeros(mask.shape)

	if (iregion != 0):
            ilabel = N.where(masklabelled == iregion)
            mask2[ilabel] = 1.
            print "length 1", len(ilabel[0])

        if (iregion2 !=0):  
            ilabel2 = N.where(masklabelled == iregion2)
            mask2[ilabel2] = 1.
            print "length 2", len(ilabel2[0])
	

        if (N.sum(mask2) != 0): mask = mask2           

        mask = 1 - mask
        mask = mf.array_shift(mask,xshift,yshift)
        

        #plt.close()
	#implot = plt.imshow(mask*N.abs(image)**0.3)
	#implot = plt.imshow(masklabelled + circ_upper + circ_lower)
        #plt.draw()

        # calculate image and mask sums
        imagesum += image
	h5sum += image*mask
        masksum += mask

        print N.max(image*mask), N.min(image*mask)

        ## Outputing the array as a png.
        outpng = outpath+description+str(count)+".png"
        pltimsave(outpng,N.abs(N.abs(image*mask)**0.3))

        # I should also display the mask in the loop
    # end of the loop over files

    #plt.close()
    #implot = plt.imshow(N.abs(h5sum)**0.3)
    #plt.draw()

    imask = N.where(masksum[:,:] > 0)
    av = N.zeros(s)
    av[imask] = h5sum[imask]/masksum[imask]

    time.sleep(1)
    #plt.close()
    #implot = plt.imshow(N.abs(av)**0.3)
    #plt.draw()

    # test writing a h5
    h5file2 = h5py.File(outpath+cleanedaveragename+".h5",'w')
    datagroup = h5file2.create_group("data")
    testdataset = datagroup.create_dataset("data",s,'f')
    testdataset[...] = av[:,:] # N.abs(av[:,:])**0.3
    h5file2.close()
	
    outpng = outpath+cleanedaveragename+".png"
    pltimsave(outpng,N.abs(N.abs(av)**0.3))

    h5file2 = h5py.File(outpath+maskname+".h5",'w')
    datagroup = h5file2.create_group("data")
    testdataset = datagroup.create_dataset("data",s,'f')
    testdataset[...] = masksum[:,:]
    h5file2.close()

    outpng = outpath+maskname+".png"
    pltimsave(outpng,masksum)

    h5file2 = h5py.File(outpath+h5sumname+".h5",'w')
    datagroup = h5file2.create_group("data")
    testdataset = datagroup.create_dataset("data",s,'f')
    testdataset[...] = h5sum[:,:]
    h5file2.close()

    outpng = outpath+h5sumname+".png"
    pltimsave(outpng,N.abs(N.abs(h5sum)**0.3))

    h5file2 = h5py.File(outpath+imagesumname+".h5",'w')
    datagroup = h5file2.create_group("data")
    testdataset = datagroup.create_dataset("data",s,'f')
    testdataset[...] = imagesum[:,:]
    h5file2.close()

    outpng = outpath+imagesumname+".png"
    pltimsave(outpng,N.abs(N.abs(imagesum)**0.3))


    #Clean up the display
    time.sleep(5)
    #cleanup_display()

    print "Finished!"
   
## END of find_streak


# a simple code to make a file list from a directory
def MakeFileList(path,ext):
    dirlist = os.listdir(path)
    filelist = []
    for file in dirlist:
         fext = os.path.splitext(file)
         #print fext
         if fext[1] == ext: 
             filelist.append(os.path.join(path,file))
    
    if len(filelist) == 0: filelist = -1
    return filelist


# initialise dynamic display
#def init_display():
#    plt.ion()
#    plt.hold(False)

#def cleanup_display():
#    plt.close()
#    plt.ioff()

# dilate by performing the autocorrelation
def dilate_fft(im):
    
    fim = N.fft.fft2(im)
    im2 = N.fft.ifft2(N.abs(fim)**2)

    tol = 0.1
    imask = N.where(N.abs(im2) > tol)

    mask = N.zeros(im2.shape,dtype="bool")
    mask[imask] = 1
    mask = mf.array_shift(mask,mask.shape[0]/2,mask.shape[1]/2)
    return mask

# crop a function
def crop_array(array,newx,newy,cenx,ceny):
    new_array = N.zeros([newx,newy])
    s = size(new_array)
    print s
    s[:]= s[:]/2
    print s
    print "lengths", len(new_array), len(array[cenx-s[0]+1:cenx+s[0]+1,ceny-s[1]+1:ceny+s[1]+1])
    new_array[:,:] = array[cenx-s[0]+1:cenx+s[0]+1,ceny-s[1]+1:ceny+s[1]+1]
    return new_array

# pad an array
def pad_array(array,newx,newy,cenx,ceny,value=0):
    new_array = N.zeros([newx,newy]) + value
    s = size(array)
    new_array[cenx-s[0]+1:cenx+s[0]+1,ceny-s[1]+1:ceny+s[1]+1] = array[:,:]
    return new_array

def size(image):
    s  = N.zeros(2)
    s[:] = image.shape[:]
    return s
    
find_streak()
