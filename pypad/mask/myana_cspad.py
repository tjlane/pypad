#--------------------------------------------------------------------------
# File and Version Information:
#  $Id: myana.py 1014 2010-05-13 18:29:22Z salnikov $
#
# Description:
#  Module myana...
#
#------------------------------------------------------------------------

"""User analysis job for XTC data.

This software was developed for the LUSI project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: myana.py 1014 2010-05-13 18:29:22Z salnikov $

@author Andrei Salnikov
    
//  Each "Element" represents one quadrant of a complete detector
//  and they are arranged as follows (viewed from upstream):
//  +---+---+
//  | 0 | 1 |
//  +---+---+
//  | 3 | 2 |
//  +---+---+
//
//  Each "Element" is composed of 8 "Section"s arranged as follows:
//  +---+---+-------+
//  |   |   |   6   |
//  + 5 | 4 +-------+
//  |   |   |   7   |
//  +---+---+---+---+   (for quadrant 0)
//  |   2   |   |   |
//  +-------+ 0 | 1 |
//  |   3   |   |   |
//  +-------+---+---+
//  The layout of each successive quadrant is rotated 90 degrees clockwise
//  with respect to the previous quadrant.
//
//  Each "Section" is composed of 2*194 rows by 185 columns with the following
//  orientations (for quadrant 0):
//    Sections 0,1: row index increases from bottom to top, column index increases from left to right
//    Sections 2,3: row index increases from left to right, column index increases from top to bottom
//    Sections 4,5: row index increases from top to bottom, column index increases from right to left
//    Sections 6,7: row index increases from left to right, column index increases from top to bottom
//  Again, the orientations of the Sections for quadrant 1 are rotated 90 degrees clockwise
//  and so on for each successive quadrant.

//----------------------------------------------------------------added by JF:
//  +---+---+-------+
//  |   |   |   4   |
//  + 3 | 2 +-------+
//  |   |   |   5   |
//  +---+---+---+---+   (for quadrant 1)
//  |   0   |   |   |
//  +-------+ 6 | 7 |
//  |   1   |   |   |
//  +-------+---+---+

//  +---+---+-------+
//  |   |   |   2   |
//  + 1 | 0 +-------+
//  |   |   |   3   |
//  +---+---+---+---+   (for quadrant 2)
//  |   6   |   |   |
//  +-------+ 4 | 5 |
//  |   7   |   |   |
//  +-------+---+---+

//  +---+---+-------+
//  |   |   |   0   |
//  + 7 | 6 +-------+
//  |   |   |   1   |
//  +---+---+---+---+   (for quadrant 3)
//  |   4   |   |   |
//  +-------+ 2 | 3 |
//  |   5   |   |   |
//  +-------+---+---+

for the Nilsson beamtime, the following modules were missing
quad 1 - section 5
quad 3 - section 4
quad 3 - section 7"""





#------------------------------
#  Module's version from CVS --
#------------------------------
__version__ = "$Revision: 1014 $"
# $Source$

#--------------------------------
#  Imports of standard modules --
#--------------------------------
import sys
import logging
import numpy as n
import pylab as p
import time
import matplotlib 
import matplotlib.pyplot as plt
import pickle
from Make_mask import Make_mask
from numpy.ma import masked_array
#---------------------------------
#  Imports of base class module --
#---------------------------------

#-----------------------------
# Imports for other modules --
#-----------------------------
from pypdsdata import xtc


#----------------------------------
# Local non-exported definitions --
#----------------------------------

#------------------------
# Exported definitions --
#------------------------

#---------------------
#  Class definition --
#---------------------
class myana_cspad( object ) :
    """Example analysis class which reads waveform data and fill a 
    profile histogram with the waveforms """

    #--------------------
    #  Class variables --
    #--------------------
    
    
    #----------------
    #  Constructor --
    #----------------
    def __init__ ( self,mask_file,dark_file,threshold=0):
        """Constructor. The parameters to the constructor are passed
        from pyana.cfg file. If parameters do not have default values
        here then the must be defined in pyana.cfg. All parameters are
        passed as strings, convert to correct type before use."""
        self.n = 0
        self.mask_file=mask_file
        self.dark_file=dark_file
        self.threshold=float(threshold)
         
        #sections are the basic unit that are read out (1 section == 2 ASICs)
        self.sect_rows = 185
        self.sect_cols = 2*194
        
        #number of sections per element and number of elements
        self.num_elems = 4
        self.num_sects_per_elem = 8
        self.total_num_sects = self.num_elems * self.num_sects_per_elem
        
        #total number of rows and columns to be read out for the whole detector
        #there are 8 sections in each element, 4 quadrants in total --> 32 sections
        self.det_raw_rows = self.num_sects_per_elem * self.sect_rows    # 8 * 185
        self.det_raw_cols = self.num_elems * self.sect_cols             # 4 * 388
        
        #set up data structures for running sums over many events
        #images 
        #array dimensions (4, 8, section rows, section colums)
        self.running_sect_sum = n.zeros( (self.num_elems, self.num_sects_per_elem, self.sect_rows, self.sect_cols) )
        
        ####set up plotting 
        #turn 'interactive' plot mode on
    
    ####################################################################
    def beginjob( self, evt, env ) :
        """This method is called once at the beginning of the job"""

        logging.info( "myana.beginjob() called" )

        #####BEGIN BLOCK FROM ORIGINAL myana_cspad.py
        #self.config = env.getConfig(TypeId.Type.Id_CspadConfig, self.address)
        self.config = env.getConfig(xtc.TypeId.Type.Id_CspadConfig)
        if not self.config:
            logging.info( '*** cspad config object is missing ***')
            return
            
        if not (self.num_elems == self.config.numQuads()):
            logging.info( "*** number of elements does not match ***")
            logging.info( "***    num_elems = %d, numQuads =  %d *** ", (self.num_elems, self.config.numQuads()))
            return

        logging.info( "Cspad configuration")
        logging.info( "  N quadrants   : %d" % self.config.numQuads())
        logging.info( "  Quad mask     : %#x" % self.config.quadMask())
        logging.info( "  payloadSize   : %d" % self.config.payloadSize())
        logging.info( "  badAsicMask0  : %#x" % self.config.badAsicMask0())
        logging.info( "  badAsicMask1  : %#x" % self.config.badAsicMask1())
        logging.info( "  asicMask      : %#x" % self.config.asicMask())
        logging.info( "  numAsicsRead  : %d" % self.config.numAsicsRead())
        #####END BLOCK FROM ORIGINAL myana_cspad.py
    

    
    ####################################################################
    def beginrun( self, evt, env ) :
        """This method is called at the beginning of the run"""
        
        
        
    ####################################################################        
    def event( self, evt, env ) :
        """This method is called for every event"""
        
        #count the events
        self.normalize=1.
        self.n+=1
        if (self.n%100 == 0):
           logging.info( "\revent number: %d ", self.n)
                
        #getCsPadQuads returns list of objects of type pypdsdata.cspad.ElementV1 
        #CSPAD full name at cxi: "CxiDs1-0|Cspad-0"
        elemarray = evt.getCsPadQuads("CxiDs1-0|Cspad-0",env)
# loop on the quadrants
        for j_elem, elem in enumerate(elemarray):
            #image data as 3-dimentional array
            #looks something like this: (8, 185, 388)
            nquadr=elem.quad()
            quadr = elem.data()
            i_sec = 0
            for sect in quadr:
                while i_sec not in self.config.sections(nquadr):
                    i_sec = i_sec+1
                #add section to running section sum (for this section number)
                self.running_sect_sum[nquadr, i_sec] += sect/self.normalize
                i_sec = i_sec+1
                
    ####################################################################            
    def endrun( self, env ) :
        """This method is called at the end of the run"""
        
        #calculate average scattering pattern
        sect_avg = self.running_sect_sum / self.n*self.normalize       # 3D: (num_sect, rows, cols)
        nquad,nsect,row,col,= sect_avg.shape
        sect_mask=self.mask4D_raw_image(sect_avg)
            
    ####################################################################           
    def endjob( self, env ) :
        """This method is called at the end of the job, close your files"""
    ####################################################################               
    def mask4D_raw_image(self,image4d):
        img_raw = n.zeros((self.det_raw_rows, self.det_raw_cols))
        #go through the 32 sections, 1st level: 4 elements, 2nd level: 8 sections
        for j_quadr, quadr in enumerate(image4d):
            for i_sec, sect in enumerate(image4d[j_quadr]):

                #compose section data to a bigger matrix of 8 by 4 elements
                xlow = self.sect_rows * i_sec
                xhigh = self.sect_rows * (i_sec+1)
                ylow = self.sect_cols * j_quadr
                yhigh = self.sect_cols * (j_quadr+1)
                #write section to current raw image
                img_raw[xlow:xhigh, ylow:yhigh] = image4d[j_quadr][i_sec]
        dark=pickle.load(open(self.dark_file))
        img_raw=img_raw-dark
        mask_raw=Make_mask(img_raw,self.mask_file,self.threshold)
