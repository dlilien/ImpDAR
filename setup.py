#! /usr/bin/env python3
# vim:fenc=utf-8
#
# Copyright © 2024 dlilien <dlilien@iu.edu>
#
# Distributed under terms of the GNU-GPL3.0 license.

"""

"""
from setuptools import setup, Extension
import numpy
ext_modules = [Extension("impdar.lib.migrationlib.mig_cython",
                         sources=["src/impdar/lib/migrationlib/_mig_cython.pyx"],
                         include_dirs=[numpy.get_include()],
                         optional=True),
               Extension("impdar.lib.ApresData.coherence",
                         sources=["src/impdar/lib/ApresData/_coherence.pyx"],
                         include_dirs=[numpy.get_include()],
                         optional=True)]

setup(ext_modules=ext_modules)
