#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Load ApRES Diff data

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

July 11 2022
"""

import glob
import numpy as np

from .load_apres import load_apres
from . import ApresData, ApresTimeDiff
from .ApresFlags import TimeDiffFlags
from ..ImpdarError import ImpdarError


def load_time_diff(fn, load_single_acquisitions=True, *args, **kwargs):
    """
    Load two Apres acquisitions together into a single object for time differencing.
    """

    if not load_single_acquisitions:
        diff_data = ApresTimeDiff(fn)
    else:
        # Make sure we have the correct fn format to load
        times = ['time1', 'time2']
        if isinstance(fn, str):
            fns = [glob.glob(fn + '_{:s}*'.format(t)) for t in times]
            for t, f in zip(times, fns):
                if len(f) != 1:
                    raise FileNotFoundError('Need exactly one file matching each time acqusition')
        elif len(fn) == 2:
            fns = fn
        else:
            raise ValueError('fn must be a glob for files with _time1, _time2, or a 2-tuple')

        # Load each of the individual acquisitions
        # as their own ApresData object
        if isinstance(fns[0], str):
            single_acquisitions = [load_apres([f]) for f in fns]
        else:
            single_acquisitions = [dat for dat in fns]

        print(single_acquisitions)

        # Check that the data have gone through the initial processing steps
        # If they haven't do range conversion and stack to one trace
        for i, acquisition in enumerate(single_acquisitions):
            try:
                acquisition.stacking()
                print('Restacked acquisition #{:d} to a 1-d array.'.format(i + 1))
            except ImpdarError:
                print('Acquisition #{:d} is already stacked to shape: {:s}'.format(i + 1,
                                          str(np.shape(acquisition.data))))
            if acquisition.flags.range == 0:
                print('Acquisition #', i+1, 'has not been converted to range. Range conversion now...')
                acquisition.apres_range(2)

        # Check that all four acquisitions have the same attributes
        from copy import deepcopy
        dat1 = deepcopy(single_acquisitions[0])
        dat2 = deepcopy(single_acquisitions[1])
        if dat1.snum != dat2.snum:
            raise ValueError('Need the same number of vertical samples in each file')
        if not np.all(dat1.travel_time == dat2.travel_time):
            raise ValueError('Need matching travel time vectors')

        # load into the ApresTimeDiff object
        diff_data = ApresTimeDiff(None)
        diff_data.snum = dat1.snum
        diff_data.data = dat1.data.flatten().astype(complex)
        diff_data.data2 = dat2.data.flatten().astype(complex)
        diff_data.decday = dat1.decday
        diff_data.range = dat1.Rcoarse
        diff_data.dt = dat1.dt
        diff_data.travel_time = dat1.travel_time

        diff_data.fn1 = dat1.header.fn
        diff_data.fn2 = dat2.header.fn
        diff_data.fn = diff_data.fn1+'_diff_'+diff_data.fn2

        # if uncertainties were calculated, bring those in too
        if hasattr(dat1, 'uncertainty'):
            diff_data.unc1 = dat1.uncertainty
        if hasattr(dat2, 'uncertainty'):
            diff_data.unc2 = dat2.uncertainty

        diff_data.data_dtype = diff_data.data.dtype

        # take the flags and header from the first data object
        diff_data.flags = TimeDiffFlags()
        diff_data.flags.file_read_code = dat1.flags.file_read_code
        diff_data.header = dat1.header

    return diff_data
