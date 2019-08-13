#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Install ImpDAR. Try to compile C sources for fast mig, else default to pure python.
"""
import setuptools

try:
    import numpy as np
    from numpy.distutils.core import Extension
except ImportError:
    raise ImportError('Numpy is required during build. Install numpy, then retry')

# For now, you have to set this manually to rebuild the c sources
CYTHON = False

if __name__ == '__main__':
    console_scripts = ['impdar=impdar.bin.impdarexec:main',
                       'impproc=impdar.bin.impproc:main',
                       'imppick=impdar.bin.imppick:main',
                       'impplot=impdar.bin.impplot:main']

    ext = '.pyx' if CYTHON else '.c'
    ext_modules = [Extension("impdar.lib.migrationlib.mig_cython",
                             sources=["impdar/lib/migrationlib/_mig_cython" + ext,
                                      "impdar/lib/migrationlib/mig_cython.c"],
                             include_dirs=[np.get_include()])]
    if CYTHON:
        from Cython.Build import cythonize
        ext_modules = cythonize(ext_modules)

    try:
        setuptools.setup(name='impdar',
                         version='0.5a',
                         description='Scripts for impulse radar',
                         url='http://github.com/dlilien/impdar',
                         author='David Lilien',
                         author_email='dal22@uw.edu',
                         license='GNU GPL-3.0',
                         entry_points={'console_scripts': console_scripts},
                         ext_modules=ext_modules,
                         install_requires=['numpy>1.12.0',
                                           'scipy>1.0.0',
                                           'matplotlib>2.0.0',
                                           'h5py'],
                         packages=['impdar',
                                   'impdar.lib',
                                   'impdar.bin',
                                   'impdar.gui',
                                   'impdar.lib.load',
                                   'impdar.lib.RadarData',
                                   'impdar.lib.migrationlib'],
                         test_suite='nose.collector')
    except SystemExit:
        print('Failed to compile c-sources. Using pure python version')
        setuptools.setup(name='impdar',
                         version='0.5a',
                         description='Scripts for impulse radar',
                         url='http://github.com/dlilien/impdar',
                         author='David Lilien',
                         author_email='dal22@uw.edu',
                         license='GNU GPL-3.0',
                         entry_points={'console_scripts': console_scripts},
                         install_requires=['numpy>1.12.0',
                                           'scipy>1.0.0',
                                           'matplotlib>2.0.0',
                                           'h5py'],
                         packages=['impdar',
                                   'impdar.lib',
                                   'impdar.bin',
                                   'impdar.gui',
                                   'impdar.lib.load',
                                   'impdar.lib.RadarData',
                                   'impdar.lib.migrationlib'],
                         test_suite='nose.collector')
