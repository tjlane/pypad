autogeom
========

Algorithms for Automatically Reconstructing CSPAD Geometries
<br>
TJ Lane  ||  <tjlane@stanford.edu>  ||  Doniach Scattering Group
<br>
Autogeom is currently v0.0.1 alpha.

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
* The GOSPEL of CSPAD, provided in `autogeom/doc`. This is an aggregation of useful information that TJ put together, included detailed descriptions of the rather confusing CSPAD geometry, the way CSPAD gets mapped from memory into real space, and how to the pyana/psana parameter sets work (amongst other things).

You can also pass -h to the scripts to get some help.


Install
-------

Installation should be as easy as downloading this package and running:

`python setup.py install`

This will automatically install the API in your python distribution and put the scripts in your path. Autogeom depends on the follow python packages:

* numpy
* scipy
* matplotlib
* pyyaml
* tables


How It Works
------------


Tutorial
--------
Below, I'll go through a quick example of how to use the autogeom scripts, and demonstrate their function via a relevant example: shifting the CSPAD quads to gain an optimal geometry. For this example we'll use an averaged run from an experiment done on gold nanoparticles at CXI. This file is provided in `tutorial/gold_avg.npz` in numpy zipped format.

Note: before you begin with any geometry optimization, you need a decently good geometry to begin with. Our example of this is contained in the set of parameters in `tutorial/cspad_params`. Providing better tools to get to this initial geometry are in the works, but for many purposes the default values should work well (and these are made available in the scripts).

(1) Generate a filter that will remove experimental noise:

`calibrate-filter tutorial/gold_avg.npz --param-dir tutorial/cspad_params/`

You'll get an interactive window with something that looks like this:

![calibrate-filter screen](https://raw.github.com/tjlane/autogeom/master/tutorial/images/calibrate-filter.png)

Follow the instructions that get printed to screen. You'll want to select regions for optimization that contain bright/sharp powder rings. This is demonstrated in the example, where I've selected the brightest ring for these gold nanoparticles.


(2) Optimize the geometry:

`optgeom tutorial/gold_avg.npz filter_params.yaml --param-dir tutorial/cspad_params/`

![assemble output](https://raw.github.com/tjlane/autogeom/master/tutorial/images/optgeom.png)



(3) Score your optimized geometry:

`score --energy 9394.363725 --path-len 129.0 --cal-type Au --param-dir  my_cspad_params tutorial/gold_avg.npz`

this should generate a directory `my_cspad_params`, containing your new parameters.

(4) Take a look at the result!

`assemble -param-dir my_cspad_params`

You should get something like this:

![assemble output](https://raw.github.com/tjlane/autogeom/master/tutorial/images/assembled-gold.png)

It may vary slightly depending on the parameters you used.


Contribute & Future
-------------------

If you would like to add to autogeom or suggest an enhancement, please do! The usual GitHub features (issues & pull requests) should suffice for these purposes. Any questions, contact TJ (email above).

There are a few items already on the to do list:

* Provide interactive tools for getting to a decent initial geometry manually.
* Write unit and integration tests.
* Clean and unify code in `cspad.py`, where redundancies exist from code copied from pyana.


