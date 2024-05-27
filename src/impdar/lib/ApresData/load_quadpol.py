#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Load quad-polarized data

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Oct 13 2020
"""

import numpy as np
import glob
import datetime
from scipy.io import loadmat
from . import ApresQuadPol
from .ApresFlags import QuadPolFlags
from .load_apres import load_apres
from ..ImpdarError import ImpdarError


def load_quadpol(fn, ftype='mat', load_single_pol=True, *args, **kwargs):
    """Load processed apres profiles from all four polarizations: hh, hv, vh, vv
    into one data object for a quad polarized acquisition.

    Parameters
    ----------
    fn: either  string of site name where fn = site name + polarization name
        or      list of file names in the correct order - hh, hv, vh, vv

    Returns
    -------
    quadpol_data: class
        quad-polarized apres data object
    """

    if not load_single_pol:
        quadpol_data = ApresQuadPol(fn)
    else:
        # Load each of the individual polarizations as their own ApresData object
        polarizations = ['HH', 'HV', 'VH', 'VV']
        if isinstance(fn, str):
            fns = [glob.glob(fn + '_{:s}.*'.format(pol)) for pol in polarizations]
            for pol, f in zip(polarizations, fns):
                if len(f) != 1:
                    raise FileNotFoundError('Need exactly one file matching each polarization')
            fns = np.squeeze(fns)
        elif len(fn) == 4:
            fns = fn
        else:
            raise ValueError('fn must be a glob for files with _HH, _HV, etc., or a 4-tuple')
        single_acquisitions = [load_apres([f]) for f in fns]

        # Check that the data have gone through the initial processing steps
        # If they haven't do range conversion and stack to one trace
        for i, acquisition in enumerate(single_acquisitions):
            try:
                acquisition.stacking()
                print('Restacked acquisition #{:d} to a 1-d array.'.format(i + 1))
            except ImpdarError:
                print('Acquisition #{:d} is already stacked to shape: {:s}'.format(i + 1, str(np.shape(acquisition.data))))
            if acquisition.flags.range == 0:
                print('Acquisition #',i+1,'has not been converted to range. Range conversion now...')
                acquisition.apres_range(2)

        # Check that all four acquisitions have the same attributes
        from copy import deepcopy
        hh = deepcopy(single_acquisitions[0])
        for xx in single_acquisitions[1:]:
            if hh.snum != xx.snum:
                raise ValueError('Need the same number of vertical samples in each file')
            if not np.all(hh.travel_time == xx.travel_time):
                raise ValueError('Need matching travel time vectors')
            if abs(hh.decday[0] - xx.decday[0]) > 1.:
                # TODO: ask to proceed
                Warning('It looks like these acquisitions were not all taken on the same day.')

        # load into the ApresQuadPol object
        quadpol_data = ApresQuadPol(None)
        quadpol_data.snum = hh.snum
        quadpol_data.shh = hh.data.flatten().astype(np.cdouble)
        quadpol_data.shv = single_acquisitions[1].data.flatten().astype(np.cdouble)
        quadpol_data.svh = single_acquisitions[2].data.flatten().astype(np.cdouble)
        quadpol_data.svv = single_acquisitions[3].data.flatten().astype(np.cdouble)
        quadpol_data.decday = hh.decday
        quadpol_data.range = hh.Rcoarse
        quadpol_data.dt = hh.dt
        quadpol_data.travel_time = hh.travel_time

        # save a 'data' field to make this class play nice with the save scripts etc.
        quadpol_data.data = quadpol_data.shh.copy()
        quadpol_data.data_dtype = quadpol_data.data.dtype

        # take the flags and header from the first data object
        quadpol_data.flags = QuadPolFlags()
        quadpol_data.flags.file_read_code = single_acquisitions[0].flags.file_read_code
        quadpol_data.header = single_acquisitions[0].header

    return quadpol_data


# ------------------------------------------------


def load_quadpol_fujita(model_name):
    """
    Load processed modeled apres profiles from all four polarizations: hh, hv, vh, vv
    into one data object

    Parameters
    ----------
    model : class
        output from the effective medium model

    Returns
    -------
    quadpol_data: class
        quad-polarized apres data object
    """

    if type(model_name) is str:
        class empty:
            pass
        model = empty()
        data = loadmat(model_name)
        for attr in data.keys():
            setattr(model,attr,np.squeeze(data[attr]))
    else:
        model = model_name

    # load into the ApresQuadPol object
    quadpol_data = ApresQuadPol(None)
    quadpol_data.shh = model.shh
    quadpol_data.shv = model.shv
    quadpol_data.svh = model.svh
    quadpol_data.svv = model.svv
    quadpol_data.range = model.range

    now = datetime.datetime.now()
    timezero = datetime.datetime(1, 1, 1, 0, 0, 0)
    offset = now-timezero
    quadpol_data.decday = offset.days + offset.seconds/(3600.*24.) + 366. # Matlab compatable

    quadpol_data.snum = len(model.shh)
    v = model.c/np.sqrt(model.epsr)
    quadpol_data.travel_time = quadpol_data.range/v
    quadpol_data.dt = np.mean(np.gradient(quadpol_data.travel_time))

    quadpol_data.data_dtype = quadpol_data.shh.dtype

    return quadpol_data
