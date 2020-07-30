#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com
#
# Distributed under terms of the GNU GPL3.0 License.

"""
Load data from the CReSIS radar MCoRDS
"""

import datetime
import numpy as np
from scipy.io import loadmat
from ..RadarData import RadarData

try:
    from netCDF4 import Dataset
    NC = True
except ImportError:
    NC = False

try:
    import h5py
    H5 = True
except ImportError:
    H5 = False



def load_mcords_nc(fn):
    """
    Load MCoRDS data as netcdf downloaded from the NSIDC

    Parameters
    ----------
    fn_nc: str
        The filename to load
    """

    mcords_data = RadarData(None)
    mcords_data.fn = fn

    if not NC:
        raise ImportError('Cannot load MCoRDS without netcdf4')
    dst = Dataset(fn, 'r')
    mcords_data.data = dst.variables['amplitude'][:].T
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
    mcords_data.tnum, mcords_data.snum = int(size[1]), int(size[0])
    mcords_data.trace_num = np.arange(mcords_data.tnum) + 1

    mcords_data.chan = 0
    mcords_data.pressure = np.zeros_like(dst.variables['lat'][:])
    mcords_data.trig = np.zeros_like(dst.variables['lat'][:]).astype(int)
    mcords_data.trig_level = 0.

    mcords_data.check_attrs()

    return mcords_data


def load_mcords_mat(fn_mat):
    """
    Load MCoRDS data as .mat format downloaded from the CReSIS ftp client

    Parameters
    ----------
    fn_mat: str
        The filename to load
    """

    mcords_data = RadarData(None)
    mcords_data.fn = fn_mat

    try:
        mat = loadmat(fn_mat)
    except:
        mat = h5py.File('AR_20140424_03_018.mat', 'r')

    if ('Data' not in mat) or ('Longitude' not in mat):
        if ('data' in mat) and ('long' in mat):
            raise KeyError('It appears that this mat file is ImpDAR/StoDeep, not MCoRDS')
        else:
            raise KeyError('ImpDAR cannot read this type of mat file--it does not appear to be MCoRDS')
    mcords_data.data = 10.*np.log10(np.squeeze(mat['Data']))
    mcords_data.long = np.squeeze(mat['Longitude'])
    mcords_data.lat = np.squeeze(mat['Latitude'])
    # sometimes the mcords data array is transposed, so check and fix
    if np.shape(mcords_data.data)[1] != np.shape(mcords_data.lat) and np.shape(mcords_data.data)[0] == np.shape(mcords_data.lat)[0]:
        mcords_data.data = np.transpose(mcords_data.data)
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

