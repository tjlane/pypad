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

    def __init__(self):
        """Constructor."""
        print 'Start constructor'

        self.pixelSize = 109.92

        self.points_for_quadrants = [ [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26],
                                      [ 6, 2,14,10,18,22,30,26] ]

        #Base index for 2x1:
        #             0  1   2  3   4   5   6   7
        self.ibase = [5, 1, 13, 9, 17, 21, 29, 25] 

        self.readOpticalAlignmentFile()

        self.get_center_coordinates()

        self.get_length_width_angle()

        self.drawOpticalAlignmentFile()
        self.drawQuadsSeparately()

    #-------------------
    # methods --
    #-------------------

    def readOpticalAlignmentFile(self): 
        print 'readOpticalAlignmentFile()'


                                 # quad 0:3
                                   # point 1:32
                                      # record: point, X, Y, Z 0:3
        self.arr = numpy.zeros( (4,33,4), dtype=numpy.int32 )

        #fname = '2011-03-29-CSPAD2-Alignment-PostRun3.txt'
        #fname = '2011-06-20-CSPAD2-Alignment-Before-Run4.txt'
        #fname = '2011-08-10-Metrology.txt'
        #fname = '2011-08-DD-Run4-DSD-Metrology.txt'
        #fname = '2012-01-10-Run5-DSD-Metrology.txt'
        #fname = '2012-01-12-Run5-DSD-Metrology-corrected.txt'
        #fname = '2012-02-26-CSPAD-XPP-Metrology.txt'
        #fname = '2012-11-08-Run6-DSD-Metrology-standard.txt'
        #fname = '2013-01-24-CSPAD-XPP-Metrology-standard.txt'
        #fname = 'metrology_renumerated.txt'
        fname = 'metrology_standard.txt'

        self.fname_plot_quads = 'metrology_standard_quads.png'
        self.fname_plot_det   = 'metrology_standard_det.png'

        file = open(fname, 'r')
        # Print out 7th entry in each line.
        for line in file:

            if len(line) == 1 : continue # ignore empty lines
            #print len(line),  ' Line: ', line

            list_of_fields = line.split()

            if list_of_fields[0] == 'Quad' : # Treat quad header lines
                self.quad = int(list_of_fields[1])
                print 'Stuff for quad', self.quad  
                continue

            if list_of_fields[0] == 'Sensor' : # Treat the title lines
                print 'Comment line:', line  
                continue
            
            if len(list_of_fields) != 4 : # Ignore lines with non-expected number of fields
                print 'len(list_of_fields) =', len(list_of_fields),
                print 'RECORD IS IGNORED due to unexpected format of the line:',line
                continue              

            point = int(list_of_fields[0])
            X = int(list_of_fields[1])
            Y = int(list_of_fields[2])
            Z = int(list_of_fields[3])
            #Title = list_of_fields[4]
            
            #record = [point, X, Y, Z, Title]
            print 'ACCEPT RECORD:', point, X, Y, Z #, Title

            self.arr[self.quad,point,0] = point
            self.arr[self.quad,point,1] = X
            self.arr[self.quad,point,2] = Y
            self.arr[self.quad,point,3] = Z

        file.close()

        print 'Array of alignment info:\n', self.arr


    def convert_optic_to_my_coordinates(self,quad,Xopt,Yopt,Zopt) :
    # DEPRICATED: New version has the same orientation as in optical measurement

        self.Zmy = Zopt

        if quad == 0 :
            self.Xmy = Xopt - 10796
            self.Ymy = Yopt - 10469 

        if quad == 1 :
            self.Xmy = Xopt - 10737   
            self.Ymy = Yopt - 10452    

        if quad == 2 :
            self.Xmy = Xopt - 11369
            self.Ymy = Yopt - 10688     

        if quad == 3 :
            self.Xmy = Xopt - 10786
            self.Ymy = Yopt - 10451    

        return ( self.Xmy, self.Ymy, self.Zmy )




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


    def print_formatted_array(self, arr, title='Array', format='%7.2f,') :
        print title

        for row in range(4) :
            for col in range(8) :
                print format % (arr[row][col]),
            print ' '


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


#----------------------------------

    def drawOpticalAlignmentFile(self): 
        print 'drawOpticalAlignmentFile()'

        sizex, sizey = shape = (100,100)
        #arr   = np.arange(sizex*sizey)
        #arr.shape = shape
        #arr   = np.zeros(shape)
        fig   = plt.figure(figsize=(10,10), dpi=100, facecolor='w',edgecolor='w',frameon=True)
        axes  = fig.add_subplot(111)        
        axes.set_xlim((-50,1750))
        axes.set_ylim((-50,1750))
        #axes1 = plt.imshow(arr, origin='lower', interpolation='nearest',extent=ax_range) 

        for quad in range(4) :
            #print '\nQuad:', quad
            self.drawOneQuad(quad,axes)

        plt.show()
        fig.savefig(self.fname_plot_det)
        print 'Image saved in file:', self.fname_plot_det


    def drawOneQuad(self,quad,axes):
        print 'drawOneQuad(' + str(quad) + ')'

        line_point = 0
        self.xlp = [0,0,0,0,0]
        self.ylp = [0,0,0,0,0]
        for point in range(1,33) :
            N = self.arr[quad,point,0]
            X = self.arr[quad,point,1]
            Y = self.arr[quad,point,2]
            Z = self.arr[quad,point,3]                
            #print 'N,X,Y =', N,X,Y

            x = self.xlp[line_point] = X / self.pixelSize
            y = self.ylp[line_point] = Y / self.pixelSize
            plt.text(x, y, str(N), fontsize=7, color='k', ha='left', rotation=45)

            if N==1 :
                x, y = self.xlp[0] + 100, self.ylp[0] + 100
                plt.text(x, y, 'Quad:'+str(quad), fontsize=12, color='k', ha='left', rotation=0)

            if line_point == 3 :
                #print 'Add new line:'
                #print 'x=',self.xlp                   
                #print 'y=',self.ylp
                self.xlp[4] = self.xlp[0]
                self.ylp[4] = self.ylp[0]
                line = lines.Line2D(self.xlp, self.ylp, linewidth=1, color='r')        
                axes.add_artist(line)
                line_point = -1
                self.xlp = [0,0,0,0,0]
                self.ylp = [0,0,0,0,0]
            line_point += 1

#----------------------------------

    def drawQuadsSeparately(self): 
        print 'drawQuadsSeparately()'

        sizex, sizey = shape = (100,100)
        fig   = plt.figure(figsize=(10,10), dpi=100, facecolor='w',edgecolor='w',frameon=True)

        quadlims = (-50,870)
        
        for quad in range(4) :
            axes = fig.add_subplot(221+quad)
            axes.set_xlim(quadlims)
            axes.set_ylim(quadlims)
            self.drawOneQuad(quad,axes)

        plt.show()
        fig.savefig(self.fname_plot_quads)
        print 'Image saved in file:', self.fname_plot_quads


#----------------------------------

def main():
    run = Main()
    sys.exit()

if __name__ == '__main__':
    main()

#----------------------------------
