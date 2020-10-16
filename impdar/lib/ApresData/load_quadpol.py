#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Load ApRES data

This code is based on a series of Matlab scripts from Craig Stewart,
Keith Nicholls, and others.
The ApRES (Automated phase-sensitive Radio Echo Sounder) is a self-contained
instrument from BAS.

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 23 2019

"""

import numpy as np
from . import QuadPolData, ApresData
from ._ApresDataProcessing import stacking,apres_range

# -----------------------------------------------------------------------------------------------------

def load_quadpol(fn, *args, **kwargs):
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

    # Load each of the individual polarizations as their own ApresData object
    single_acquisitions = []
    if type(fn) is str:
        single_acquisitions.append(ApresData(fn+'_HH.mat'))
        single_acquisitions.append(ApresData(fn+'_HV.mat'))
        single_acquisitions.append(ApresData(fn+'_VH.mat'))
        single_acquisitions.append(ApresData(fn+'_VV.mat'))
    elif hasattr(fn,'__len__') and len(fn) == 4:
        # TODO: Ask the user to check that the files are correct
        single_acquisitions.append(ApresData(fn[0]))
        single_acquisitions.append(ApresData(fn[1]))
        single_acquisitions.append(ApresData(fn[2]))
        single_acquisitions.append(ApresData(fn[3]))

    # Check that the data have gone through the initial processing steps
    # If they haven't do range conversion and stack to one trace
    for i,xx in enumerate(single_acquisitions):
        print('Restacking acquisition #',i+1,'to a 1-d array...')
        stacking(xx)
        if xx.flags.range[0] == 0:
            print('Acquisition #',i+1,'has not been converted to range. Range conversion now...')
            apres_range(xx,2)

    # Check that all four acquisitions have the same attributes
    from copy import deepcopy
    hh = deepcopy(single_acquisitions[0])
    for xx in single_acquisitions[1:]:
        if hh.snum != xx.snum:
            raise ValueError('Need the same number of vertical samples in each file')
        if not np.all(hh.travel_time == xx.travel_time):
            raise ValueError('Need matching travel time vectors')
        if not np.all(abs(hh.decday - xx.decday)<1.):
            # TODO: ask to proceed
            raise ValueError('It looks like these acquisitions were not all taken on the same day.')

    # load into the QuadPolData object
    quadpol_data = QuadPolData(None)
    quadpol_data.snum = hh.snum
    quadpol_data.shh = hh.data.flatten()
    quadpol_data.shv = single_acquisitions[1].data.flatten()
    quadpol_data.svh = single_acquisitions[2].data.flatten()
    quadpol_data.svv = single_acquisitions[3].data.flatten()
    quadpol_data.decday = hh.decday
    quadpol_data.dt = hh.dt
    quadpol_data.travel_time = hh.travel_time

    quadpol_data.data_dtype = quadpol_data.shh.dtype

    return quadpol_data
