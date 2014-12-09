#!/usr/bin/local/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as plt_patches

from pypad import cspad
from pypad import utils
from pypad import read
from pypad import plot

cs1 = cspad.CSPad.load("examples/my_cspad.cspad")
cs2 = cspad.CSPad.default()
cs3 = cspad.CSPad.default()
cs3.quad_rotation[0] = 10.

raw_image = read.load_raw_image("examples/gold-minus490mm.h5")

fig = plt.figure()
ax = plt.subplot(121)
plot.imshow_cspad(cs3(raw_image), ax=ax, scrollable=True)
ax = plt.subplot(122)
plot.sketch_2x1s(cs3.pixel_positions, ax)
plt.show()
