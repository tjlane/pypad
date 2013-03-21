PyPad
========

PyPad provides a software-independent and intuitive way to:

* Automatically Optimize CSPAD Geometries
* Interactive Mask Generation
* CSPAD Visualization
* CSPAD Analytics (Energy, Pixel Statistics, etc.)

<br>
Principle Authors:
TJ Lane        ||  <tjlane@stanford.edu>         ||  Doniach Group
Jonas Sellberg ||  <sellberg@slac.stanford.edu>  ||  Nilsson Group
<br>
PyPad is currently v0.0.1 alpha.

--------------------------------------------------------------------------------


Install
-------

Installation should be as easy as downloading this package and running:

`python setup.py install`

This will automatically install the API in your python distribution and put the scripts in your path. Autogeom depends on the follow python packages:

* numpy
* scipy
* matplotlib
* pyyaml
* pytables


Functionality
-------------

PyPad is divided into models that each aim to address a different challenge in analyzing CSPad data:

* autogeom: automatically refines a CSPad geometry to
* mask:     an easy and interactive platform for masking out bad or suspect pixels
* analytics: gain access to pixel statistics, energy fluctuations, or other simple analytics

Find out more about each of these modules by checking out the README files at

`pypad/pypad/<module_name>/README.md`

Further, PyPad provides an easy way to visualize your data, and convert CSPad images between different software packages (including psana, Cheetah, myana, ODIN, etc).


Contribute
----------

If you would like to add to PyPad or suggest an enhancement, please do! The usual GitHub features (issues & pull requests) should suffice for these purposes.

There are a few items already on the to do list:

* Provide interactive tools for getting to a decent initial geometry manually.
* Write unit and integration tests.
* Clean and unify code in `cspad.py`, where redundancies exist from code copied from pyana.

There is also a HISTORY file in the master branch that details a version history, along with a nitty-gritty to do list.
