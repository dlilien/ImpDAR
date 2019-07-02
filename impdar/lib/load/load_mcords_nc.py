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
    mcords_data.time = dst.variables['time'][:]
    mcords_data.travel_time = dst.variables['fasttime'][:]
    mcords_data.dt = np.mean(np.diff(mcords_data.travel_time)) * 1.0e-6
    size = dst.variables['amplitude'].matlab_size
    mcords_data.snum, mcords_data.tnum = int(size[0]), int(size[1])
    mcords_data.check_attrs()
    return mcords_data
