#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
# Worked from minimal example from scipy-examples.org

"""Fast Migration. Just Kirchhoff for now"""
import sys
cimport numpy as np
np.import_array()
import numpy as np
import time

# cdefine the signature of our c function
# Need this so that the function is recognized
cdef extern from "mig_cython.h":
    void mig_kirch_loop (double * migdata, int tnum, int snum, double * dist, double * zs, double * zs2, double * tt_sec, double vel, double * gradD, double max_travel_time, int nearfield)


def migrationKirchhoffLoop(np.ndarray[double, ndim=2, mode="c"] migdata not None,
                          int tnum,
                          int snum,
                          np.ndarray[double, ndim=1, mode="c"] dist not None,
                          np.ndarray[double, ndim=1, mode="c"] zs not None,
                          np.ndarray[double, ndim=1, mode="c"] zs2 not None,
                          np.ndarray[double, ndim=1, mode="c"] tt_sec not None,
                          float vel,
                          np.ndarray[double, ndim=2, mode="c"] gradD not None,
                          float max_travel_time,
                          bint nearfield
                          ):
    """I am not sure if this wrapper is needed, but I think it gives us type checking so I'm leaving it"""
    mig_kirch_loop(<double*> np.PyArray_DATA(migdata),
                   tnum,
                   snum,
                   <double*> np.PyArray_DATA(dist),
                   <double*> np.PyArray_DATA(zs),
                   <double*> np.PyArray_DATA(zs2),
                   <double*> np.PyArray_DATA(tt_sec),
                   vel,
                   <double*> np.PyArray_DATA(gradD),
                   max_travel_time,
                   int(nearfield)
                   )


def migrationKirchhoff(dat, vel=1.69e8, nearfield=False):
    """Kirchhoff Migration (Berkhout 1980; Schneider 1978; Berryhill 1979)

    This migration method uses an integral solution to the scalar wave equation Yilmaz (2001) eqn 4.5.
    The algorithm cycles through every sample in each trace, creating a hypothetical diffraciton
    hyperbola for that location,
        t(x)^2 = t(0)^2 + (2x/v)^2
    To migrate, we integrate the power along that hyperbola and assign the solution to the apex point.
    There are two terms in the integral solution, Yilmaz (2001) eqn 4.5, a far-field term and a
    near-field term. Most algorithms ignore the near-field term because it is small. Here there is an option,
    but default is to ignore.

    Parameters
    ---------
    dat: data as a class in the ImpDAR format
    vel: wave velocity, default is for ice
    nearfield: boolean to indicate whether or not to use the nearfield term in summation

    Output
    ---------
    dat: data as a class in the ImpDAR format (with dat.data now being migrated data)

    """

    print('Kirchhoff Migration (diffraction summation) of %.0fx%.0f matrix' % (dat.snum, dat.tnum))
    print('Using compiled cython/c version')
    # check that the arrays are compatible
    _check_data_shape(dat)
    # start the timer
    start = time.time()
    # Calculate the time derivative of the input data
    gradD = np.gradient(np.ascontiguousarray(dat.data, dtype=np.float64), dat.travel_time / 1.0e6, axis=0)
    # Create an empty array to fill with migrated data
    migdata = np.ascontiguousarray(np.zeros_like(dat.data, dtype=np.float64), dtype=np.float64)

    # Try to cache some variables that we need lots
    tt_sec = dat.travel_time / 1.0e6
    max_travel_time = np.max(tt_sec)

    # Cache the depths
    zs = vel * tt_sec / 2.0
    zs2 = zs**2.
    migrationKirchhoffLoop(migdata,
                          dat.tnum,
                          dat.snum,
                          np.ascontiguousarray(dat.dist, dtype=np.float64) * 1.0e3,
                          np.ascontiguousarray(zs, dtype=np.float64),
                          np.ascontiguousarray(zs2, dtype=np.float64),
                          np.ascontiguousarray(tt_sec, dtype=np.float64),
                          vel,
                          np.ascontiguousarray(gradD, dtype=np.float64),
                          max_travel_time,
                          nearfield
                          )
    dat.data = migdata.copy()
    # print the total time
    print('Kirchhoff Migration of %.0fx%.0f matrix complete in %.2f seconds\n'
          % (dat.snum, dat.tnum, time.time() - start))
    return dat


def _check_data_shape(dat):
    if np.size(dat.data, 1) != dat.tnum or np.size(dat.data, 0) != dat.snum:
        raise ValueError('The input array must be of size (tnum,snum)')
