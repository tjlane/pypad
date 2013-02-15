
import tables

import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import filters

import assemble, tjlane


def find_edges(image, threshold=0.01, minf_size=3, medf_size=10):
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
    
    
f = tables.File('hdf5_examples/AgBe/r0003-RawSum.h5')
i = find_edges(f.root.data.data.read())

plt.imshow(i.T)
plt.show()

i2 = tjlane.load_AgBe()
ai = assemble.assemble_image_from_dir(i2, 'example_calibration_dir')
i2f = find_edges(ai)
plt.imshow(i2f.T, interpolation='nearest')
plt.show()
