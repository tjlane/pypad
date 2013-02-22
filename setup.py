#coding: utf8

"""
Setup script autogeom
"""

from glob import glob

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(name='autogeom',
      version='0.0.1',
      author="TJ Lane",
      author_email="tjlane@stanford.edu",
      description='CSPad geometry, assembly and optimization',
      packages=["autogeom"],
      package_dir={"autogeom": "autogeom" },
      scripts=[s for s in glob('scripts/*') if not s.endswith('__.py')],
      test_suite="test")