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
Below, I'll go through a quick example of how to use the autogeom scripts, and demonstrate their function via a relevant example: shifting the CSPAD quads to gain an optimal geometry. For this example we'll use an averaged run from an experiment done on gold nanoparticles at CXI. This file is provided in `examples/gold_image.h5` in HDF5 format.

Note: before you begin with any geometry optimization, you need a decently good geometry to begin with. Our example of this is contained in the set of parameters in `examples/cspad_params`. You may have to improve these parameters further for your specific experiment, but for many purposes the default values should work well (and these are made available in the scripts if you run without specifying `--param-dir`).


(1) Generate a filter that will remove experimental noise:

`genfilter examples/gold_image.h5 --param-dir examples/cspad_params/`

You'll get an interactive window with something that looks like this:

![genfilter screen](https://raw.github.com/tjlane/pypad/master/doc/images/genfilter-init.png)

Follow the instructions that get printed to screen. You'll want to change the threshold and filters to obtain a good contrast in the 2D image, then set a rough center and select regions for optimization that contain bright/sharp powder rings. You can select one or several regions, which include one or several powder rings each. This is demonstrated in the example below, where I've selected the two brightest rings (Miller indices (111) and (200)) for these gold nanoparticles:

![genfilter screen](https://raw.github.com/tjlane/pypad/master/doc/images/genfilter-opt.png)


(2) Optimize the geometry:

`optgeom examples/gold_image.h5 filter_params.yaml --param-dir examples/cspad_params/`

![assemble output](https://raw.github.com/tjlane/pypad/master/doc/images/optgeom.png)


(3) Score your optimized geometry:

`score examples/gold_image.h5 --energy 9603.224 --path-len 81.052 --cal-type Au --param-dir my_cspad_params`

This assembles your scattering pattern using the optimized CSPAD alignment parameters specified by `--param-dir` and scores it given the calibration sample, photon energy, and sample-to-detector distance are given as input parameters.


(4) Take a look at the result:

`assemble examples/gold_image.h5 --param-dir my_cspad_params`

You should get something like this:

![assemble output](https://raw.github.com/tjlane/pypad/master/doc/images/assembled-gold.png)

It may vary slightly depending on the parameters you used.


(5) Export the parameters to your front-end of choice:

'geom2cheetah my_cspad_params'

this should generate a Cheetah pixel map called `pixelmap-cheetah-raw.h5`, containing your new parameters. It is now ready to use as input geometry for the Cheetah.

