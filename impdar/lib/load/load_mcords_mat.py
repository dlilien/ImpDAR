#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

Load a MCoRDS data that was downloaded as a .mat file from the CReSIS ftp client.

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

July 17 2019

"""

import datetime
import numpy as np
from ..RadarData import RadarData
from scipy.io import loadmat


def load_mcords_mat(fn_mat):
    """Load MCoRDS data in .mat format downloaded from the CReSIS ftp client

    Parameters
    ----------
    fn: str
        The filename to load
    """
    mcords_data = RadarData(None)

    mat = loadmat(fn_mat)
    if ('Data' not in mat) or ('Longitude' not in mat):
        if ('data' in mat) and ('long' in mat):
            raise KeyError('It appears that this mat file is ImpDAR/StoDeep, not MCoRDS')
        else:
            raise KeyError('ImpDAR cannot read this type of mat file--it does not appear to be MCoRDS')
    mcords_data.data = 10.*np.log10(np.squeeze(mat['Data']))
    mcords_data.long = np.squeeze(mat['Longitude'])
    mcords_data.lat = np.squeeze(mat['Latitude'])
    # time has units of seconds according to documentation, but this seems wrong
    # numbers are way too big. Leaving it since that is how it is documented though?
    partial_days = np.squeeze(mat['GPS_time']) / (24. * 60. * 60.)
    start_day = datetime.datetime(1970,1,1,0,0,0).toordinal() + 366.
    mcords_data.decday = partial_days + start_day
    mcords_data.trace_int = mcords_data.decday[1] - mcords_data.decday[0]
    mcords_data.travel_time = np.squeeze(mat['Time'])*1e6
    mcords_data.dt = np.mean(np.diff(mcords_data.travel_time)) * 1.0e-6
    size = np.shape(mcords_data.data)
    mcords_data.snum, mcords_data.tnum = int(size[0]), int(size[1])
    mcords_data.trace_num = np.arange(mcords_data.tnum) + 1

    mcords_data.chan = 0
    mcords_data.pressure = np.zeros_like(mcords_data.decday)
    mcords_data.trig = np.zeros_like(mcords_data.decday).astype(int)
    mcords_data.trig_level = 0.
    mcords_data.check_attrs()
    return mcords_data
