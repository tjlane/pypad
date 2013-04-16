
"""
Provides a "mask" object for CSPads
"""

import numpy as np
import matplotlib.pyplot as plt


class PadMask(object):
    """
    An mask for a CSPad object.
    """
    
    
    def __init__(self):
        """
        Initialize a CSPad mask object.
        """
        self._masks = {}
        self._masks['base'] = self._blank_mask()
        return
    
    
    @property
    def mask(self):
        return np.product( np.array(self._masks.values()) )
    
    
    @property
    def inverted(self):
        """
        Invert the mask. Usually "True" is a good pixel, and "False" is a bad
        one, but this flips them.
        
        Returns
        -------
        """
        return np.logical_not(self.mask)
        
    
    @property
    def num_masked(self):
        """
        Returns the number of masked pixels
        """
        inv_mask = 1.0 - self.mask
        return np.sum(inv_mask)
        
        
    @property
    def types_applied(self):
        return self._masks.keys()
    
    
    def remove_mask(self, mask_name):
        """
        Remove a mask that has been applied.
        """
        
        if mask_name == 'base': # this one is special
            self._masks['base'] = self._blank_mask()
        elif not mask_name in self._masks.keys():
            raise KeyError('Mask: %s not applied' % mask_name)
        else:
            x = self._masks.pop(mask_name)
            
        print "Removed mask: %s" % mask_name
            
        return
    
        
    def _inject_mask(self, mask_name, mask):
        """
        Add a new kind of mask to this mask object. Provides some typechecking.
        
        All this really does is deposit `mask` into dict self._masks with key
        `mask_name`.
        """
        
        assert type(mask_name) == str
        assert mask.shape == (4, 8, 185, 388)
        
        if not mask.dtype == np.bool:
            mask = mask.astype(np.bool)
            
        self._masks[mask_name] = mask
        
        return
        
        
    def _check_image(self, image):
        """
        Sanity check on `image`.
        """
        if not image.shape == (4, 8, 185, 388):
            raise ValueError('`image` must be shape (4, 8, 185, 388), got '
                             '%s' % str(image.shape))
        return
        
        
    def _blank_mask(self):
        """
        Utility function that just returns a blank mask.
        """
        return np.ones((4, 8, 185, 388), dtype=np.bool)
    
        
    # ----------
    # below we provide many methods of the form PadMask.mask_*(), which allow
    # one to mask pixels via various criteria (indicated by the *)
    
    # to add a new kind of mask, make a new method here. Follow mask_threshold
    # as a template
    
    def mask_pixel(self, quad, two_by_one, x, y):
        """
        Mask a single pixel, or series of pixels. To do the latter, pass arrays
        as the arguments (even though the below says int).
        
        Parameters
        ----------
        quad : int
            [0,1,2,3], denoting the quad.
            
        two_by_one : int
            Int in [0,7], denoting 2x1.
            
        x : int
            Int in [0,184], denoting x position.
        
        y : int
            Int in [0,388], denoting x position.
        """
        self._masks['base'][quad, two_by_one, x, y] = 0
        return
    
        
    def unmask_pixel(self, quad, two_by_one, x, y):
        """
        Mask a single pixel, or series of pixels. To do the latter, pass arrays
        as the arguments (even though the below says int).
        
        Parameters
        ----------
        quad : int
            [0,1,2,3], denoting the quad.
            
        two_by_one : int
            Int in [0,7], denoting 2x1.
            
        x : int
            Int in [0,184], denoting x position.
        
        y : int
            Int in [0,388], denoting x position.
        """
        self._masks['base'][quad, two_by_one, x, y] = 1
        return
    
        
    def mask_threshold(self, image, upper=None, lower=None):
        """
        Mask pixels by threshold values, in ADU. You must supply either `upper`
        or `lower`, but one one is required.
        
        Parameters
        ----------
        image : np.ndarray
            A shape (4, 8, 185, 388) array describing a CSPad image.
        
        upper : int
            Values greater than this are masked.
            
        lower : int
            Values lower than this are masked.
        """
        
        upper = int(upper)
        lower = int(lower)
        
        if (upper == None) and (lower == None):
            raise ValueError('Either `upper` or `lower` (or both) must be specified')
        
        self._check_image(image)
        
        m = self._blank_mask()
        ind = (image > upper) + (image < lower)
        m[ind] = np.bool(False)
        
        self._inject_mask('threshold', m)
        
        return
    
        
    def mask_nonbonded(self, nearest_neighbours=True):
        """
        Mask pixels on the CSPad that were never bonded.
        
        Optional Parameters
        -------------------
        nearest_neighbours : bool
            Also mask four of their nearest neighbours, which give anomoulous 
            responses.
        """
        
        m = self._blank_mask()
        
        for i in range(4):
            for j in range(8):
        		for p in range(0, 185, 10):
        		    
        		    if nearest_neighbours:
            		    xl = max(0, p-1)
            		    xh = min(184, p+1)
            		    yl = max(0, p-1)
            		    yh = min(387, p+1)
                        m[xl:xh,yl:yh] = 0
                    else:
                        m[p,p] = 0
        			
        self._inject_mask('nonbonded', m)
        
        return
    
        
    def mask_borders(self, num_pixels=1):
        """
        Mask the border of each ASIC, to a width of `num_pixels`.
        
        Parameters
        ----------
        num_pixels : int
            The size of the border region to mask.
        """
        
        n = int(num_pixels)        
        m = self._blank_mask()
        
        if (num_pixels < 0) or (num_pixels > 194):
            raise ValueError('`num_pixels` must be >0, <194')
        
        for i in range(4):
            for j in range(8):
                
                # mask along the y-dim
                m[i,j,:,0:n] = np.bool(False)
                m[i,j,:,388-n:388] = np.bool(False)
                
                # mask along the x-dim
                m[i,j,0:n,:] = np.bool(False)
                m[i,j,185-n:185,:] = np.bool(False)
                
                # mask a bar along y in the middle of the 2x1
                m[i,j,:,194-n:194+n] = np.bool(False)
        
        self._inject_mask('border', m)        
        
        return
    
    # ----------
        
    def merge(self, *args):
        """
        Merge two or more masks, masking with an OR operator for masked pixels.
        """
        
        for mask in args:
            for mtype in mask._masks.keys():
                if mtype in self._masks.keys():
                    self._masks[mtype] = np.logical_not( np.logical_or(self._masks[mtype], 
                                                           mask._masks[mtype]) )
                else:
                    self._masks[mtype] = mask._masks[mtype]
        
        return
    
        
class InteractiveMask(object):

   def __init__(self, data, mask_file, threshold=0, row13=False):
      self.key=[]  
      self.x=0
      self.y=0
      self.xy=[]
      self.xx=[]
      self.yy=[]
      self.data=data
      self.lx,self.ly=p.shape(self.data)
      self.points=[]
      for i in range(self.lx):
          for j in range(self.ly):
           self.points.append([i,j]) 
      self.mask_file=mask_file
      if os.path.exists(self.mask_file) is True:
#         if 'edf' in self.mask_file:
#            mask_f=EdfFile.EdfFile(self.mask_file)
#            self.mymask=mask_f.GetData(0)
#            num_masks=mask_f.GetNumImages()
#            if num_masks==2:
#               self.automask=mask_f.GetData(1)
#               self.anisotropic_mask=0*self.mymask
#            if num_masks==3:
#               self.automask=mask_f.GetData(1)
#               self.anisotropic_mask=mask_f.GetData(2)
#            else:
#               self.automask=0*self.mymask
#               self.anisotropic_mask=0*self.mymask
#            if p.shape(self.mymask)!=p.shape(self.data):
#               self.mymask=n.zeros((self.lx,self.ly))
#         elif 'h5' in self.mask_file:
         if 'h5' in self.mask_file:
            newfile=self.make_name
            self.open_mask()
            self.anisotropic_mask=0*self.mymask
      else:
         self.mymask=n.zeros((self.lx,self.ly))
         self.automask=n.zeros((self.lx,self.ly))
         self.anisotropic_mask=n.zeros((self.lx,self.ly))
      self.old_mymask=self.mymask
      self.old_automask=self.automask
      self.automask[n.where(self.data<=threshold)]=1
      print "automatically masking out " + str(int(self.automask.sum())) + " pixels below or equal to threshold=%s" % (threshold)
      #this is for masking out row13 of the CSPAD
      if (row13):
         print "automatically masking out row13 of CSPAD"
         col=181
         for i in range(8):
            self.automask[:,col]=1
            col+= 194
      #end of CSPAD part
      palette=p.cm.jet
      palette.set_bad('w',1.0)
      p.rc('image',origin = 'lower')
      p.rc('image',interpolation = 'nearest')
      p.figure(2)
      self.px=p.subplot(111)
      self.data=p.log(self.data+1)
      self.im=p.imshow(masked_array(self.data,self.mymask+self.automask+self.anisotropic_mask), cmap=palette)
      p.title('Select a ROI. Press m to mask it or u to unmask it. k to save/exit, q to exit without saving')
      self.lc,=self.px.plot((0,0),(0,0),'-+m',linewidth=1,markersize=8,markeredgewidth=1)
      self.lm,=self.px.plot((0,0),(0,0),'-+m',linewidth=1,markersize=8,markeredgewidth=1)
      self.px.set_xlim(0,self.ly)
      self.px.set_ylim(0,self.lx)
      self.colorbar=p.colorbar(self.im,pad=0.01)
      cidb=p.connect('button_press_event',self.on_click)
      cidk=p.connect('key_press_event',self.on_click)
      cidm=p.connect('motion_notify_event',self.on_move)
      p.show()
      
   def on_click(self,event):
       if not event.inaxes: 
           self.xy=[]
           return
       self.x,self.y=int(event.xdata), int(event.ydata)
       self.key=event.key
       self.xx.append([self.x])
       self.yy.append([self.y])
       self.xy.append([self.y,self.x])
       self.lc.set_data(self.xx,self.yy)
       if self.key=='m': 
           print 'masking'
           self.xx[-1]=self.xx[0]
           self.yy[-1]=self.yy[0]
           self.xy[-1]=self.xy[0]
           ind=p.nonzero(points_inside_poly(self.points,self.xy))
           self.mymask=self.mymask.reshape(self.lx*self.ly,1)
           self.mymask[ind]=1
           self.mymask=self.mymask.reshape(self.lx,self.ly)
           datamasked=masked_array(self.data,self.mymask+self.automask+self.anisotropic_mask)
           self.im.set_data(datamasked)
           self.xx=[]
           self.yy=[]
           self.xy=[] 
           self.lc.set_data(self.xx,self.yy)
           self.lm.set_data(self.xx,self.yy)
#           self.im.set_clim(vmax=(2*self.data.mean()))
           self.im.autoscale()
           p.draw()
           self.x=0
           self.y=0 
       if self.key=='u':
           print 'unmasking'
           self.xx[-1]=self.xx[0]
           self.yy[-1]=self.yy[0]
           self.xy[-1]=self.xy[0]
           ind=p.nonzero(points_inside_poly(self.points,self.xy))
           self.mymask=self.mymask.reshape(self.lx*self.ly,1)
           self.mymask[ind]=0
           self.mymask=self.mymask.reshape(self.lx,self.ly)
           datanew=masked_array(self.data,self.mymask+self.automask+self.anisotropic_mask)

           self.im.set_data(datanew)
           self.xx=[]
           self.yy=[]
           self.xy=[]
           self.lc.set_data(self.xx,self.yy)
           self.lm.set_data(self.xx,self.yy)
#           self.im.set_clim(vmax=(2*self.data.mean()))
           self.im.autoscale()
           p.draw()
           self.x=0
           self.y=0

       if self.key=='r':
           print 'unmasking all'
           self.mymask=0*self.mymask
           datanew=masked_array(self.data,self.mymask+self.automask+self.anisotropic_mask)
           self.im.set_data(datanew)
           self.xx=[]
           self.yy=[]
           self.xy=[] 
           self.lc.set_data(self.xx,self.yy)
           self.lm.set_data(self.xx,self.yy)

#           self.im.set_clim(vmax=(2*self.data.mean()))
           self.im.autoscale()
           p.draw()
           self.x=0
           self.y=0 
       if self.key=='k':
          print 'save and exit'
          self.save_mask()
          print 'Mask saved in file:', self.mask_file
#          mask_f=EdfFile.EdfFile(self.mask_file)
#          mask_f.WriteImage({},self.mymask,0)
#          mask_f.WriteImage({},self.automask,1)
#          mask_f.WriteImage({},self.anisotropic_mask,2)
#          del(mask_f)
          p.close()
          return self.mymask+self.automask
       if self.key=='q':
          print 'exit without saving'
          p.close()
          return self.old_mymask+self.old_automask

   def on_move(self,event):
       if not event.inaxes: return
       xm,ym=int(event.xdata), int(event.ydata)
       # update the line positions
       if self.x!=0: 
           self.lm.set_data((self.x,xm),(self.y,ym))
           p.draw()

   def save_mask(self):
       newfile=self.make_name()
       h5file = h.File(newfile,'w')
       datagroup = h5file.create_group("masks")
       dataset = datagroup.create_dataset("user",self.mymask.shape,dtype="float")
       dataset[...] = self.mymask[:,:]
       dataset2 = datagroup.create_dataset("auto",self.mymask.shape,dtype="float")
       dataset2[...] = self.automask[:,:]
       h5file.close()
       tosave=abs(n.ceil((self.mymask+self.automask)/2)-1)
       h5file = h.File(self.mask_file,'w')
       datagroup = h5file.create_group("data")
       dataset = datagroup.create_dataset("data",tosave.shape,dtype="int16")
       dataset[...] = tosave[:,:]
       h5file.close()
       
   def open_mask(self):
       newfile=self.make_name()
       h5file = h.File(newfile,'r')
       self.mymask=n.array(h5file.get("masks/user"))
       self.automask=n.array(h5file.get("masks/auto"))
       
   def make_name(self):
       dirname,filename=os.path.split(self.mask_file)
       base,ext=filename.split('.')
       newname=base+'_2masks.'+ext
       newfile=os.path.join(dirname,newname)
       return newfile