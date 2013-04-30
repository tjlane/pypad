#coding: utf8

"""
Setup script for pypad.
"""

from glob import glob


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(name='pypad',
      version='0.0.1',
      author="TJ Lane",
      author_email="tjlane@stanford.edu",
      description='CSPad geometry, assembly and optimization',
      packages=["pypad", "pypad.autogeom"],
      package_dir={"pypad": "pypad"},
      scripts=[s for s in glob('scripts/*') if not s.endswith('__.py')],
      test_suite="test")