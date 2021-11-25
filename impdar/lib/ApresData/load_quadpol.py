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
from . import QuadPolData
from .load_apres import load_apres

# -----------------------------------------------------------------------------------------------------

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
        quadpol_data = QuadPolData(fn)
    else:
        # Load each of the individual polarizations as their own ApresData object
        single_acquisitions = []
        if type(fn) is str:
            single_acquisitions.append(load_apres(glob.glob(fn+'_HH*')))
            single_acquisitions.append(load_apres(glob.glob(fn+'_HV*')))
            single_acquisitions.append(load_apres(glob.glob(fn+'_VH*')))
            single_acquisitions.append(load_apres(glob.glob(fn+'_VV*')))
        elif hasattr(fn,'__len__') and len(fn) == 4:
            # TODO: Ask the user to check that the files are correct
            single_acquisitions.append(load_apres(fn[0]))
            single_acquisitions.append(load_apres(fn[1]))
            single_acquisitions.append(load_apres(fn[2]))
            single_acquisitions.append(load_apres(fn[3]))

        # Check that the data have gone through the initial processing steps
        # If they haven't do range conversion and stack to one trace
        for i,xx in enumerate(single_acquisitions):
            try:
                xx.stacking()
                print('Restacked acquisition #',i+1,'to a 1-d array.')
            except:
                print('Acquisition #',i+1,'is already stacked to shape:',np.shape(xx.data))
            if xx.flags.range[0] == 0:
                print('Acquisition #',i+1,'has not been converted to range. Range conversion now...')
                xx.apres_range(2)

        # Check that all four acquisitions have the same attributes
        from copy import deepcopy
        hh = deepcopy(single_acquisitions[0])
        for xx in single_acquisitions[1:]:
            if hh.snum != xx.snum:
                raise ValueError('Need the same number of vertical samples in each file')
            if not np.all(hh.travel_time == xx.travel_time):
                raise ValueError('Need matching travel time vectors')
            if abs(hh.decday - xx.decday) > 1.:
                # TODO: ask to proceed
                Warning('It looks like these acquisitions were not all taken on the same day.')

        # load into the QuadPolData object
        quadpol_data = QuadPolData(None)
        quadpol_data.snum = hh.snum
        quadpol_data.shh = hh.data.flatten().astype(np.complex)
        quadpol_data.shv = single_acquisitions[1].data.flatten().astype(np.complex)
        quadpol_data.svh = single_acquisitions[2].data.flatten().astype(np.complex)
        quadpol_data.svv = single_acquisitions[3].data.flatten().astype(np.complex)
        quadpol_data.decday = hh.decday
        quadpol_data.range = hh.Rcoarse
        quadpol_data.dt = hh.dt
        quadpol_data.travel_time = hh.travel_time

        quadpol_data.data_dtype = quadpol_data.shh.dtype

    return quadpol_data
