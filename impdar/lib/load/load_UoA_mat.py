#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""
Load an output file processed from the University of Alabama VHF or UWB radars (2019).

Data are all processed, but some more complex data views are not available.
This allows easy colorbar adjustment, interpolation to geospatial coordinates, and layer picking.
"""

import numpy as np
from scipy.interpolate import interp1d
from ..RadarData import RadarData
from ..gpslib import nmea_info
import h5py


def load_UoA_mat(fn_mat):
    """Load MCoRDS data in .mat format downloaded from the CReSIS ftp client

    Parameters
    ----------
    fn: str
        The filename to load
    """
    UoA_data = RadarData(None)

    with h5py.File(fn_mat) as fin:
        UoA_data.data = 10*np.log10(np.sqrt((fin['Data']['channel1'][:, :]['real'] ** 2.0 + (fin['Data']['channel1'][:, :]['imag']) ** 2.0).T))
        UoA_data.snum, UoA_data.tnum = int(UoA_data.data.shape[0]), int(UoA_data.data.shape[1])
        UoA_data.trace_num = np.arange(UoA_data.tnum) + 1
        UoA_data.travel_time = fin['Data']['time1'][:].flatten() * 1.0e6
        UoA_data.dt = np.mean(np.diff(UoA_data.travel_time)) * 1.0e-6
        nminfo = nmea_info()
        nminfo.time = fin['INS_GPS']['POSIX_time'][:].flatten() / (24. * 60. * 60.)  # Stored in seconds
        nminfo.ppstime = fin['INS_GPS']['pps_cntr'][:].flatten()
        nminfo.lat = fin['INS_GPS']['latitude'][:].flatten()
        nminfo.lon = fin['INS_GPS']['longitude'][:].flatten()
        nminfo.elev = fin['INS_GPS']['altitude_MSL'][:].flatten()
        nminfo.get_utm()
        nminfo.get_dist()

        UoA_data.lat = interp1d(nminfo.ppstime, nminfo.lat)(fin['Data']['PPS_Time'][:].flatten())
        UoA_data.long = interp1d(nminfo.ppstime, nminfo.lon)(fin['Data']['PPS_Time'][:].flatten())
        UoA_data.x_coord = interp1d(nminfo.ppstime, nminfo.x)(fin['Data']['PPS_Time'][:].flatten())
        UoA_data.y_coord = interp1d(nminfo.ppstime, nminfo.y)(fin['Data']['PPS_Time'][:].flatten())
        UoA_data.dist = interp1d(nminfo.ppstime, nminfo.dist)(fin['Data']['PPS_Time'][:].flatten())
        UoA_data.elev = interp1d(nminfo.ppstime, nminfo.elev)(fin['Data']['PPS_Time'][:].flatten())
        UoA_data.decday = interp1d(nminfo.ppstime, nminfo.time)(fin['Data']['PPS_Time'][:].flatten())
        UoA_data.trace_int = UoA_data.decday[1] - UoA_data.decday[0]
        UoA_data.pressure = np.zeros_like(UoA_data.decday)
        UoA_data.trig = np.zeros_like(UoA_data.decday).astype(int)
        UoA_data.trig_level = 0.
        if fn_mat[-10:] == '_files.mat':
            UoA_data.chan = 999
        else:
            if 'hannel' in fn_mat:
                idx = fn_mat.index('hannel')
                UoA_data.chan = int(fn_mat[idx + 6])
            elif 'Ch' in fn_mat:
                idx = fn_mat.index('Ch')
                UoA_data.chan = int(fn_mat[idx + 2])
            else:
                UoA_data.chan = 10
        UoA_data.check_attrs()
        return UoA_data
