#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Created for compilation of fortran code
"""
import setuptools
import numpy as np
from numpy.distutils.core import Extension, setup

try:
    from Cython.Distutils import build_ext
    CYTHON = True
except:
    CYTHON = False

if __name__ == '__main__':
    if CYTHON:
        ext_modules = [Extension("impdar.lib.mig_cython",
                                 sources=["impdar/lib/_mig_cython.pyx", "impdar/lib/mig_cython.c"],
                                 include_dirs=[np.get_include()])]
        cmdclass = {'build_ext': build_ext}
    else:
        ext_modules = None
        cmdclass = None

    console_scripts = ['impdar=impdar.bin.impdarexec:main', 'impproc=impdar.bin.impproc:main', 'imppick=impdar.bin.imppick:main', 'impplot=impdar.bin.impplot:main']
    setuptools.setup(name='impdar',
                     version='0.5a',
                     description='Scripts for impulse radar',
                     cmdclass=cmdclass,
                     url='http://github.com/dlilien/impdar',
                     author='David Lilien',
                     author_email='dal22@uw.edu',
                     license='GNU GPL-3.0',
                     entry_points={'console_scripts': console_scripts},
                     ext_modules=ext_modules,
                     install_requires=['numpy>1.12.0', 'scipy>1.0.0', 'matplotlib>2.0.0', 'h5py'],
                     packages=['impdar', 'impdar.lib', 'impdar.bin', 'impdar.gui', 'impdar.lib.load', 'impdar.lib.RadarData'],
                     test_suite='nose.collector')
