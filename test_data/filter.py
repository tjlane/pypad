
import tables

import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import filters




def find_edges(image, threshold=0.0025, minf_size=1, medf_size=8):
    """
    Applies an edge filter followed by a noise reduction filter. Very good
    at locating powder rings and filtering everything else out.
    
    Parameters
    ----------
    image : ndarray
        An image to find the edges of
        
    Returns
    -------
    binary_image : ndarray, np.bool
        A binary image, with "1" where there are powder rings/strong edges
    """
    
    image = np.abs(filters.sobel(image, 0)) + np.abs(filters.sobel(image, 1))
    image -= image.min()

    assert image.min() == 0
    assert image.max() > 0

    image = (image > (image.max() * threshold)).astype(np.bool)

    image = filters.minimum_filter(image, size=minf_size)
    image = filters.median_filter(image, size=medf_size)
        
    return image.astype(np.bool)
    
    
f = tables.File('AgBe/r0003-RawSum.h5')
i = find_edges(f.root.data.data.read())

plt.imshow(i.T)
plt.show()
