Autogeom is a module of PyPad that aims to:

* provide a pythonic interface to the CSPAD geometry, including assembling raw images for visualization, access to the pixel positions in a variety of representations (including basis vector representations).
* optimize CSPAD geometries based on calibration standards, correcting e.g. for quad movements at the CXI hutch.
* assess the quality of a CSPAD geometry via known calibration standards
* interact with CSPAD optical metrologies and convert those metrologies into parameters specifying a CSPAD geometry
* provide access to all the above via a clean and intuitive API and CLUI


Help
----

There are three sources of documentation for autogeom:

* This README
* The autogeom API is documented in-place in the code.
* You can also pass -h to the scripts to get some help.


Tutorial
--------
Below, I'll go through a quick example of how to use the autogeom scripts, and demonstrate their function via a relevant example: shifting the CSPAD quads to gain an optimal geometry. For this example we'll use an averaged run from an experiment done on gold nanoparticles at CXI. This file is provided in `tutorial/gold_avg.npz` in numpy zipped format.

Note: before you begin with any geometry optimization, you need a decently good geometry to begin with. Our example of this is contained in the set of parameters in `tutorial/cspad_params`. Providing better tools to get to this initial geometry are in the works, but for many purposes the default values should work well (and these are made available in the scripts).

(1) Generate a filter that will remove experimental noise:

`calibrate-filter tutorial/gold_avg.npz --param-dir tutorial/cspad_params/`

You'll get an interactive window with something that looks like this:

![calibrate-filter screen](https://raw.github.com/tjlane/pypad/master/tutorial/images/calibrate-filter.png)

Follow the instructions that get printed to screen. You'll want to select regions for optimization that contain bright/sharp powder rings. This is demonstrated in the example, where I've selected the brightest ring for these gold nanoparticles.


(2) Optimize the geometry:

`optgeom tutorial/gold_avg.npz filter_params.yaml --param-dir tutorial/cspad_params/`

![assemble output](https://raw.github.com/tjlane/pypad/master/tutorial/images/optgeom.png)



(3) Score your optimized geometry:

`score --energy 9394.363725 --path-len 129.0 --cal-type Au --param-dir  my_cspad_params tutorial/gold_avg.npz`

this should generate a directory `my_cspad_params`, containing your new parameters.

(4) Take a look at the result!

`assemble -param-dir my_cspad_params`

You should get something like this:

![assemble output](https://raw.github.com/tjlane/pypad/master/tutorial/images/assembled-gold.png)

It may vary slightly depending on the parameters you used.