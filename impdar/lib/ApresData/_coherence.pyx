#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
# Worked from minimal example from scipy-examples.org

"""Fast Coherence, just for quadpol calculation"""
import sys
cimport numpy as np
np.import_array()
import numpy as np
import time

# cdefine the signature of our c function
# Need this so that the function is recognized
cdef extern from "coherence.h":
    void coherence2d (double complex * chhvv, double complex * HH, double complex * VV, int nrange, int ntheta, int range_bins, int azimuth_bins)


def coherence2d_loop(np.ndarray[double complex, ndim=2, mode="c"] chhvv not None,
                     np.ndarray[double complex, ndim=2, mode="c"] HH not None,
                     np.ndarray[double complex, ndim=2, mode="c"] VV not None,
                     int nrange,
                     int ntheta,
                     int range_bins,
                     int azimuth_bins
                     ):
    """I am not sure if this wrapper is needed, but I think it gives us type checking so I'm leaving it"""
    coherence2d(<double complex*> np.PyArray_DATA(chhvv),
                <double complex*> np.PyArray_DATA(HH),
                <double complex*> np.PyArray_DATA(VV),
                nrange,
                ntheta,
                range_bins,
                azimuth_bins
                )
    return chhvv.copy()
