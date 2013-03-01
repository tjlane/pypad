#!/usr/bin/env python

"""
This file obstensibly contains all the necessary code to convert an optical
metrology into the psana alignment parameters.

I have not merged it into the code base yet, but that is on the to do list.

-- TJL 2.27.13
"""

#----------------------------------
import sys
import numpy
import numpy as np
import math

import matplotlib.pyplot as plt
import matplotlib.lines  as lines

#----------------------------------

class Main :
    """Main"""

    def __init__(self, metrology_file):
        """Constructor."""
        print 'Start constructor'

        self.pixelSize = 109.92
        self._convert_metrology_to_cspad(metrology_file)


    def _convert_metrology_to_cspad(self, metrology_file):
        self._preprocess_file_for_conversion_to_cspad(metrology_file)
        self._compute_center_coordinates()
        self._compute_length_width_angle()


    def _preprocess_file_for_conversion_to_cspad(self, filename, verbose=False): 
        """
        Load up the metrology file and read out key quantities. This has
        """
        
        if verbose:
            print "Reading: %s" % filename

        self.points_for_quadrants = [ [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26] ]

        #Base index for 2x1:
        #             0  1   2  3   4   5   6   7
        self.ibase = [5, 1, 13, 9, 17, 21, 29, 25]

        self.arr = numpy.zeros( (4,33,4), dtype=numpy.int32 )

        f = open(filename, 'r')
        # Print out 7th entry in each line.
        for line in f:

            if len(line) == 1:
                continue # ignore empty lines

            list_of_fields = line.split()

            if list_of_fields[0] == 'Quad': # Treat quad header lines
                self.quad = int(list_of_fields[1])
                if verbose: print 'Stuff for quad', self.quad  
                continue

            if list_of_fields[0] == 'Sensor' : # Treat the title lines
                if verbose: print 'Comment line:', line  
                continue
            
            if len(list_of_fields) != 4 : # Ignore lines with non-expected number of fields
                print 'WARNING: len(list_of_fields) =', len(list_of_fields),
                print 'RECORD IS IGNORED due to unexpected format of the line:',line
                continue              

            point = int(list_of_fields[0])
            X = int(list_of_fields[1])
            Y = int(list_of_fields[2])
            Z = int(list_of_fields[3])

            self.arr[self.quad,point,0] = point
            self.arr[self.quad,point,1] = X
            self.arr[self.quad,point,2] = Y
            self.arr[self.quad,point,3] = Z

        file.close()

        print 'Array of alignment info:\n', self.arr


    def get_center_coordinates(self) :

        print 'Coordinates of the 2x1 sensor centers:'

        self.arrXmu = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.arrYmu = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.arrZmu = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.arrX   = numpy.zeros( (4,8), dtype=numpy.float )
        self.arrY   = numpy.zeros( (4,8), dtype=numpy.float )
        self.arrZ   = numpy.zeros( (4,8), dtype=numpy.float )

        ix = 1
        iy = 2
        iz = 3

        for quad in range(4) :
            for pair in range(8) :

                icor1 = self.ibase[pair]   
                icor2 = self.ibase[pair] + 1
                icor3 = self.ibase[pair] + 2
                icor4 = self.ibase[pair] + 3

                X = 0.25 * ( self.arr[quad,icor1,ix]
                           + self.arr[quad,icor2,ix]
                           + self.arr[quad,icor3,ix]
                           + self.arr[quad,icor4,ix] )

                Y = 0.25 * ( self.arr[quad,icor1,iy]
                           + self.arr[quad,icor2,iy]
                           + self.arr[quad,icor3,iy]
                           + self.arr[quad,icor4,iy] ) 

                Z = 0.25 * ( self.arr[quad,icor1,iz]
                           + self.arr[quad,icor2,iz]
                           + self.arr[quad,icor3,iz]
                           + self.arr[quad,icor4,iz] ) 

                #Xmy, Ymy, Zmy = self.convert_optic_to_my_coordinates(quad,X,Y,Z)
                Xmy, Ymy, Zmy = X, Y, Z

                print 'quad:%1d, pair:%2d,  X:%7d  Y:%7d, Z:%3d' % (quad, pair, Xmy, Ymy, Zmy)

                self.arrXmu[quad][pair] = Xmy
                self.arrYmu[quad][pair] = Ymy
                self.arrZmu[quad][pair] = Zmy

                #self.arrX[quad][pair] = int( float(Xmy) / self.pixelSize )
                #self.arrY[quad][pair] = int( float(Ymy) / self.pixelSize )
                #self.arrZ[quad][pair] = int( float(Zmy) / self.pixelSize )

                self.arrX[quad][pair] = float(Xmy) / self.pixelSize
                self.arrY[quad][pair] = float(Ymy) / self.pixelSize
                self.arrZ[quad][pair] = float(Zmy) / self.pixelSize

        print 'My X, Y, and Z coordinates of the 2x1 CENTER:'
        print 'in um (micrometer):'
        print 'X(um):\n', self.arrXmu
        print 'Y(um):\n', self.arrYmu
        print 'Z(um):\n', self.arrZmu

        print '\nIn pixels:'
        self.print_formatted_array(self.arrX,'X(pixel):')
        self.print_formatted_array(self.arrY,'Y(pixel):')
        self.print_formatted_array(self.arrZ,'Z(pixel):')

        print '\nIn pixels (w/o comma and axis titles):'
        self.print_formatted_array(self.arrX,' ', format='%7.2f ')
        self.print_formatted_array(self.arrY,' ', format='%7.2f ')
        self.print_formatted_array(self.arrZ,' ', format='%7.2f ')


    def get_length_width_angle(self) :

        self.S1  = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.S2  = numpy.zeros( (4,8), dtype=numpy.int32 )

        self.dS1 = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.dS2 = numpy.zeros( (4,8), dtype=numpy.int32 )

        self.L1  = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.L2  = numpy.zeros( (4,8), dtype=numpy.int32 )

        self.dL1 = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.dL2 = numpy.zeros( (4,8), dtype=numpy.int32 )

        self.angDegree = numpy.zeros( (4,8), dtype=numpy.float32 )

        self.D1  = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.D2  = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.dD  = numpy.zeros( (4,8), dtype=numpy.int32 )

        self.ddS = numpy.zeros( (4,8), dtype=numpy.int32 )
        self.ddL = numpy.zeros( (4,8), dtype=numpy.int32 )

        ix = 1
        iy = 2

        print 'pair:        S1      S2     dS1     dS2        L1      L2     dL1     dL2    <dS/L>  angle(deg)      D1      D2      dD   d(dS)   d(dL)'

        for quad in range(4) :

            print '\nQuad ', quad

            for pair in range(8) :

                icor1 = self.ibase[pair]   
                icor2 = self.ibase[pair] + 1
                icor3 = self.ibase[pair] + 2
                icor4 = self.ibase[pair] + 3

                if pair == 0 or  pair == 1 or  pair == 4 or  pair == 5 :

                    self. S1[quad][pair] = self.arr[quad,icor2,iy] - self.arr[quad,icor1,iy]
                    self. S2[quad][pair] = self.arr[quad,icor3,iy] - self.arr[quad,icor4,iy]

                    self.dS1[quad][pair] = self.arr[quad,icor4,iy] - self.arr[quad,icor1,iy]
                    self.dS2[quad][pair] = self.arr[quad,icor3,iy] - self.arr[quad,icor2,iy]

                    self. L1[quad][pair] = self.arr[quad,icor4,ix] - self.arr[quad,icor1,ix]
                    self. L2[quad][pair] = self.arr[quad,icor3,ix] - self.arr[quad,icor2,ix]

                    self.dL1[quad][pair] = self.arr[quad,icor2,ix] - self.arr[quad,icor1,ix]
                    self.dL2[quad][pair] = self.arr[quad,icor3,ix] - self.arr[quad,icor4,ix]


                else:

                    self. S1[quad][pair] =   self.arr[quad,icor4,ix] - self.arr[quad,icor1,ix]
                    self. S2[quad][pair] =   self.arr[quad,icor3,ix] - self.arr[quad,icor2,ix]
                                                                                           
                    self.dS1[quad][pair] = -(self.arr[quad,icor2,ix] - self.arr[quad,icor1,ix]) # sign is chosen 
                    self.dS2[quad][pair] = -(self.arr[quad,icor3,ix] - self.arr[quad,icor4,ix]) # for positive phi

                    self. L1[quad][pair] =   self.arr[quad,icor2,iy] - self.arr[quad,icor1,iy]
                    self. L2[quad][pair] =   self.arr[quad,icor3,iy] - self.arr[quad,icor4,iy]
                                                                                           
                    self.dL1[quad][pair] =   self.arr[quad,icor4,iy] - self.arr[quad,icor1,iy]
                    self.dL2[quad][pair] =   self.arr[quad,icor3,iy] - self.arr[quad,icor2,iy]


                diag1x = float(self.arr[quad,icor1,ix] - self.arr[quad,icor3,ix])
                diag2x = float(self.arr[quad,icor2,ix] - self.arr[quad,icor4,ix])
                diag1y = float(self.arr[quad,icor1,iy] - self.arr[quad,icor3,iy])
                diag2y = float(self.arr[quad,icor2,iy] - self.arr[quad,icor4,iy])

                self.D1[quad][pair] = int( math.sqrt(diag1x*diag1x + diag1y*diag1y) )
                self.D2[quad][pair] = int( math.sqrt(diag2x*diag2x + diag2y*diag2y) )
                self.dD[quad][pair] = self.D1[quad][pair] - self.D2[quad][pair]

                self.ddS[quad][pair] = self.dS1[quad][pair] - self.dS2[quad][pair]
                self.ddL[quad][pair] = self.dL1[quad][pair] - self.dL2[quad][pair]

                ang1 = ang2 = 0
                if self.dL1[quad][pair] != 0 : ang1 = float(self.dS1[quad][pair]) / self.L1[quad][pair]
                if self.dL2[quad][pair] != 0 : ang2 = float(self.dS2[quad][pair]) / self.L2[quad][pair]

                angav = (ang1 + ang2) * 0.5
                ang_deg = 180/3.1415927 * angav

                self.angDegree[quad][pair] = ang_deg

                print 'pair: %1d  %6d  %6d  %6d  %6d    %6d  %6d  %6d  %6d   %8.5f   %8.5f  %6d  %6d  %6d  %6d  %6d' % \
                    (pair, self.S1[quad][pair], self.S2[quad][pair], self.dS1[quad][pair], self.dS2[quad][pair], \
                           self.L1[quad][pair], self.L2[quad][pair], self.dL1[quad][pair], self.dL2[quad][pair], \
                           angav, ang_deg, \
                           self.D1[quad][pair], self.D2[quad][pair], self.dD[quad][pair], self.ddS[quad][pair], self.ddL[quad][pair])

        self.print_formatted_array(self.angDegree,'\ndPhi:','%8.5f,')
        self.print_formatted_array(self.angDegree,'\ndPhi (w/o comma):','%8.5f ')


