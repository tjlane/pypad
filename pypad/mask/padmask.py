
"""
Provides a "mask" object for CSPads
"""

import numpy as np


class PadMask(object):
    
    def __init__(self):
        """
        Initialize a CSPad mask object.
        """
        
        self._masks = {}
        
        return
    
    
    @property
    def mask(self):
        return np.product( np.array(self._masks.values() ) )
    
        
    @property
    def num_masked(self):
        """
        Returns the number of masked pixels
        """
        inv_mask = 1.0 - self.mask
        return np.sum(inv_mask)
    
        
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
    
        
    def mask_nonbonded(self):
        """
        Mask pixels on the CSPad that were never bonded.
        """
        
        # TJL to Jonas : I think I need your help on this one...
        #                not sure about the nearest neighbour search
        
        raise NotImplementedError()
        
        m = self._blank_mask()
        
        for arow in range(8):
        	for acol in range(8):
        		for n in range(0,COLS,10):
        			m[arow*COLS+n,acol*ROWS+n] = 0
        
        return
    
        
    def mask_borders(self, num_pixels):
        """
        Mask the border of each ASIC, to a width of `num_pixels`.
        
        Parameters
        ----------
        num_pixels : int
            The size of the border region to mask.
        """
        
        n = int(num_pixels)        
        m = self._blank_mask()
        
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
        Merge two masks, masking with an OR operator for masked pixels.
        """
        pass
        
        
    def invert(self):
        """
        Invert the mask. Usually "True" is a good pixel, and "False" is a bad
        one, but this flips them.
        
        Returns
        -------
        """
        return
        
        
        
class InteractiveMask(object):
    
    # TJL the below is just c/p'd, needs work.
    
	def __init__(self, inarr, filename):
		self.inarr = inarr
		self.filename = filename
		self.cmax = self.inarr.max()
		self.cmin = self.inarr.min()

	def on_keypress(self,event):
		if event.key == 'p':
			if not os.path.exists(write_dir + runtag):
				os.mkdir(write_dir + runtag)
			pngtag = write_dir + "pixel_mask_cxi74613_%s_variance_%s-%sADUs.png" % (self.filename, str(options.dead), str(options.threshold))
			print "saving image as " + pngtag 
			P.savefig(pngtag)
		if event.key == 'r':
			colmin = self.cmin
			colmax = self.cmax
			P.clim(colmin, colmax)
			P.draw()


	def on_click(self, event):
		if event.inaxes:
			lims = self.axes.get_clim()
			colmin = lims[0]
			colmax = lims[1]
			range = colmax - colmin
			value = colmin + event.ydata * range
			if event.button is 1 :
				if value > colmin and value < colmax :
					colmin = value
			elif event.button is 2 :
				colmin, colmax = self.orglims
			elif event.button is 3 :
				if value > colmin and value < colmax:
					colmax = value
			P.clim(colmin, colmax)
			P.draw()


	def draw_img(self):
		fig = P.figure()
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title('pixel_mask_'+self.filename+'_variance')
		P.rc('image',origin='lower')
		self.axes = P.imshow(self.inarr, vmin = 0, vmax = options.threshold)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		#P.show()
		pngtag = write_dir + "pixel_mask_cxi74613_%s_variance_%s-%sADUs.png" % (self.filename, str(options.dead), str(options.threshold))
		print "saving image as " + pngtag 
		P.savefig(pngtag)
		P.close()