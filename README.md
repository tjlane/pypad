autogeom
========

Algorithms for Automatically Reconstructing CSPAD Geometries
TJ Lane <tjlane@stanford.edu> | Doniach Scattering Group

--------------------------------------------------------------------------------

Autogeom is a python package that aims to:

* provide a pythonic interface to the CSPAD geometry, including assembling raw images for visualization, access to the pixel positions in a variety of representations (including basis vector representations).
* optimize CSPAD geometries based on calibration standards, correcting e.g. for quad movements at the CXI hutch.
* assess the quality of a CSPAD geometry via known calibration standards
* interact with CSPAD optical metrologies and convert those metrologies into parameters specifying a CSPAD geometry
* provide access to all the above via a clean and intuitive API and CLUI

Documentation
-------------

There are three sources of documentation for autogeom:

* This README
* The autogeom API is documented in-place in the code.
* The GOSPEL of CSPAD, provided in autogeom/doc. This is an aggregation of useful information that TJ put together, included detailed descriptions of the rather confusing CSPAD geometry, the way CSPAD gets mapped from memory into real space, and how to the pyana/psana parameter sets work (amongst other things).

You can also pass -h to the scripts to get some help.


Install
-------

Installation should be as easy as downloading this package and running:

 python setup.py install

This will automatically install the API in your python distribution and put the scripts in your path. Autogeom depends on:

* numpy
* scipy
* matplotlib
* pyyaml
* tables



Tutorial
--------

Below, I'll go through a quick example of how to use the autogeom scripts, and demonstrate their function via a relevant example: shifting the CSPAD quads to gain an optimal geometry. For this example we'll use 


(1) Generate a filter that will remove experimental noise:

  calibrate-filter tutorial/gold_avg.npz --param-dir data/ex_params/


(2) Optimize the geometry:

  optgeom tutorial/gold_avg.npz filter_params.yaml --param-dir tutorial/cspad_params/


(3) Score your optimized geometry:

  score --energy 9394.363725 --path-len 129.0 --cal-type Au --param-dir \
  my_cspad_params tutorial/gold_avg.npz

this should generate a directory `my_cspad_params`, containing your new parameters.

(4) Take a look at the result!

  assemble -param-dir my_cspad_params


