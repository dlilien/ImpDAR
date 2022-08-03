#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2020 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""Install ImpDAR, possibly with C sources."""
import setuptools
import socket

try:
    import numpy as np
    from numpy.distutils.core import Extension
except ImportError:
    raise ImportError('Numpy is required during build.')

# For now, this should just run on my computer and I'll distribute the c code
if socket.gethostname() == 'hozideh':
    CYTHON = True
else:
    CYTHON = False


if __name__ == '__main__':
    console_scripts = ['impdar=impdar.bin.impdarexec:main',
                       'impproc=impdar.bin.impproc:main',
                       'imppick=impdar.bin.imppick:main',
                       'impplot=impdar.bin.impplot:main',
                       'apdar=impdar.bin.apdar:main']

    ext = '.pyx' if CYTHON else '.c'
    ext_modules = [Extension("impdar.lib.migrationlib.mig_cython",
                             sources=["impdar/lib/migrationlib/_mig_cython"
                                      + ext,
                                      "impdar/lib/migrationlib/mig_cython.c"],
                             include_dirs=[np.get_include()])]
    if CYTHON:
        from Cython.Build import cythonize
        ext_modules = cythonize(ext_modules)

    version = '1.1.5'
    packages = ['impdar',
                'impdar.lib',
                'impdar.bin',
                'impdar.gui',
                'impdar.gui.ui',
                'impdar.lib.load',
                'impdar.lib.RadarData',
                'impdar.lib.ApresData',
                'impdar.lib.migrationlib']

    requires = ['numpy>1.12.0',
                'scipy>0.19.0',
                'matplotlib>2.0.0',
                'h5py',
                'segyio']

    try:
        setuptools.setup(name='impdar',
                         version=version,
                         description='Scripts for impulse radar',
                         url='http://github.com/dlilien/impdar',
                         author='David Lilien',
                         author_email='dal22@uw.edu',
                         license='GNU GPL-3.0',
                         entry_points={'console_scripts': console_scripts},
                         ext_modules=ext_modules,
                         install_requires=requires,
                         packages=packages,
                         long_description=open('README.md', 'r').read(),
                         long_description_content_type='text/markdown',
                         test_suite='nose.collector')
    except SystemExit:
        print('Failed to compile c-sources. Using pure python version')
        setuptools.setup(name='impdar',
                         version=version,
                         description='Scripts for impulse radar',
                         url='http://github.com/dlilien/impdar',
                         author='David Lilien',
                         author_email='dal22@uw.edu',
                         license='GNU GPL-3.0',
                         entry_points={'console_scripts': console_scripts},
                         install_requires=requires,
                         packages=packages,
                         long_description=open('README.md', 'r').read(),
                         long_description_content_type='text/markdown',
                         test_suite='nose.collector')
