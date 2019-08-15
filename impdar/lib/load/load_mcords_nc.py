#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com
#
# Distributed under terms of the GNU GPL3.0 License.

"""
Load a NetCDF of MCoRDS data
"""
import datetime
import numpy as np
from ..RadarData import RadarData

try:
    from netCDF4 import Dataset
    NC = True
except ImportError:
    NC = False


def load_mcords_nc(fn_nc):
    """Load MCoRDS data in netcdf downloaded from the NSIDC

    Parameters
    ----------
    fn_nc: str
        The filename to load
    """
    if not NC:
        raise ImportError('Cannot load MCoRDS without netcdf4')
    mcords_data = RadarData(None)
    dst = Dataset(fn_nc, 'r')
    mcords_data.data = dst.variables['amplitude'][:]
    mcords_data.long = dst.variables['lon'][:]
    mcords_data.lat = dst.variables['lat'][:]
    # time has units of seconds according to documentation, but this seems wrong
    # numbers are way too big. Leaving it since that is how it is documented though?
    partial_days = dst.variables['time'][:] / (24. * 60. * 60.)
    start_day = datetime.datetime(int(dst.variables['time'].units[14:18]),
                                  int(dst.variables['time'].units[19:21]),
                                  int(dst.variables['time'].units[22:24])).toordinal() + 366.
    mcords_data.decday = partial_days + start_day
    mcords_data.trace_int = mcords_data.decday[1] - mcords_data.decday[0]
    mcords_data.travel_time = dst.variables['fasttime'][:]
    mcords_data.dt = np.mean(np.diff(mcords_data.travel_time)) * 1.0e-6
    size = dst.variables['amplitude'].matlab_size
    mcords_data.tnum, mcords_data.snum = int(size[0]), int(size[1])
    mcords_data.trace_num = np.arange(mcords_data.tnum) + 1

    mcords_data.chan = 0
    mcords_data.pressure = np.zeros_like(dst.variables['lat'][:])
    mcords_data.trig = np.zeros_like(dst.variables['lat'][:]).astype(int)
    mcords_data.trig_level = 0.
    mcords_data.check_attrs()
    return mcords_data
