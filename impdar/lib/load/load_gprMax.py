#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

Load h5 files and convert to the .mat ImpDAR file

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Mar 17 2019

"""

import numpy as np
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags

try:
    import h5py
    H5 = True
except ImportError:
    H5 = False


def load_gprMax(fn_h5, *args, **kwargs):
    """Load a gprMax file, which is really an h5 file, into ImpDAR"""
    if not H5:
        raise ImportError('You need H5 to load gprMax')

    h5_data = RadarData(None)
    h5_data.fn = fn_h5

    # open the h5 file
    with h5py.File(fn_h5, 'r') as f_in:
        h5_data.dt = f_in.attrs['dt']
        h5_data.data = np.array(f_in['/rxs/rx1/Ez'])

    # Remove pretrigger
    trig_threshold = 0.5  # trigger when mean trace gets up to 50% of maximum
    mean_trace = np.nanmean(np.abs(h5_data.data), axis=1)
    idx_threshold = np.argwhere(mean_trace > trig_threshold * np.nanmax(mean_trace))
    idx_trig = np.nanmin(idx_threshold)
    h5_data.data = h5_data.data[idx_trig:]

    # other variables are from the array shape
    h5_data.snum = h5_data.data.shape[0]
    h5_data.tnum = h5_data.data.shape[1]
    h5_data.trace_num = np.arange(h5_data.data.shape[1]) + 1
    h5_data.trig_level = np.zeros((h5_data.tnum,))
    h5_data.pressure = np.zeros((h5_data.tnum,))
    h5_data.flags = RadarFlags()
    h5_data.travel_time = h5_data.dt * 1.0e6 * np.arange(h5_data.snum)
    h5_data.trig = np.zeros((h5_data.tnum,))
    h5_data.lat = np.zeros((h5_data.tnum,))
    h5_data.long = np.zeros((h5_data.tnum,))
    h5_data.x_coord = np.zeros((h5_data.tnum,))
    h5_data.y_coord = np.zeros((h5_data.tnum,))
    h5_data.elev = np.zeros((h5_data.tnum,))
    h5_data.decday = np.arange(h5_data.tnum)
    h5_data.trace_int = np.ones((h5_data.tnum,))

    h5_data.dist = np.arange(h5_data.tnum)
    h5_data.chan = -99.
    h5_data.check_attrs()
    return h5_data
