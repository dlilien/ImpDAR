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
setuptools

if __name__ == '__main__':
    console_scripts = ['impdar=impdar.bin.impdar:main', 'impproc=impdar.bin.impproc:main']
    setuptools.setup(name='impdar',
                     version='0.1b3',
                     description='Scripts for impulse radar',
                     url='http://github.com/dlilien/impdar',
                     author='David Lilien',
                     author_email='dal22@uw.edu',
                     license='GNU GPL-3.0',
                     entry_points={'console_scripts': console_scripts},
                     install_requires=['numpy>1.12.0', 'scipy>1.0.0', 'matplotlib>2.0.0'],
                     packages=['impdar', 'impdar.lib'],
                     test_suite='nose.collector')
