#!/usr/bin/env python
import pylab as p
import numpy as n
class My_imshow:
   def __init__(self,img):
      self.img=img
      self.fig=p.figure()
      cid2 = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
      cid1=self.fig.canvas.mpl_connect('key_press_event',self.keypress) 
      p.rc('image',origin='lower')
      self.axes = p.imshow(img)
      self.colbar = p.colorbar(self.axes, pad=0.01)
      self.orglims = self.axes.get_clim()
      p.show()
   def keypress(self,event):
      if event.inaxes:
         if event.key=='o':
            print 'change to logscale'
            self.axes.set_data(p.log(self.img))
            self.axes.autoscale()
            lims = self.axes.get_clim()
            p.clim(lims)
            p.draw()
         if event.key=='n':
            print 'change to linearscale'
            self.axes.set_data(self.img)
            self.axes.autoscale()
            lims = self.axes.get_clim()
            p.clim(lims)
            p.draw()
   def on_click(self,event):
      if event.inaxes:
           lims = self.axes.get_clim()
           colmin = lims[0]
           colmax = lims[1]
           range = colmax-colmin
           value = colmin+event.ydata*range
           if event.button is 1 :
                   if value > colmin and value < colmax :
                           colmin = value
           elif event.button is 2 :
                   colmin, colmax = self.orglims
           elif event.button is 3 :
                   if value > colmin and value < colmax:
                           colmax = value
           p.clim(colmin, colmax)
           p.draw()
       

   


      
