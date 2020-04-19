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
from ..RadarData import RadarData, RadarFlags
from ..gpslib import nmea_info
try:
    import h5py
    H5 = True
except ImportError:
    H5 = False


def load_UoA_mat(fn_mat, gps_offset=0.0):
    """Load MCoRDS data in .mat format downloaded from the CReSIS ftp client

    Parameters
    ----------
    fn: str
        The filename to load
    """
    UoA_data = RadarData(None)
    UoA_data.fn = fn_mat

    with h5py.File(fn_mat, 'r') as fin:
        UoA_data.data = fin['Data']['channel'][:, :].T
        # if len(UoA_data.data.dtype) == 2:
        #     UoA_data.data = UoA_data.data['real'] + 1.j * UoA_data.data['imag']

        # complex data
        if len(UoA_data.data.dtype) == 2:
            UoA_data.data = 10 * np.log10(np.sqrt(UoA_data.data['real'] ** 2.0 + UoA_data.data['imag'] ** 2.0))
        else:
           UoA_data.data = 10 * np.log10(UoA_data.data)
        UoA_data.snum, UoA_data.tnum = int(UoA_data.data.shape[0]), int(UoA_data.data.shape[1])
        UoA_data.trace_num = np.arange(UoA_data.tnum) + 1
        UoA_data.travel_time = fin['Data']['fast_time'][:].flatten() * 1.0e6
        UoA_data.dt = np.mean(np.diff(UoA_data.travel_time)) * 1.0e-6
        nminfo = nmea_info()
        nminfo.time = (fin['INS_GPS']['POSIX_time'][:].flatten() + gps_offset) / (24. * 60. * 60.)  # Stored in seconds
        nminfo.ppstime = fin['INS_GPS']['POSIX_time'][:].flatten() + gps_offset
        nminfo.lat = fin['INS_GPS']['latitude'][:].flatten()
        nminfo.lon = fin['INS_GPS']['longitude'][:].flatten()
        nminfo.elev = fin['INS_GPS']['altitude_MSL'][:].flatten()

        UoA_data.lat = interp1d(nminfo.ppstime, nminfo.lat, fill_value='extrapolate')(fin['Data']['POSIX_time'][:].flatten())
        UoA_data.long = interp1d(nminfo.ppstime, nminfo.lon, fill_value='extrapolate')(fin['Data']['POSIX_time'][:].flatten())
        UoA_data.elev = interp1d(nminfo.ppstime, nminfo.elev, fill_value='extrapolate')(fin['Data']['POSIX_time'][:].flatten())
        UoA_data.decday = interp1d(nminfo.ppstime, nminfo.time, fill_value='extrapolate')(fin['Data']['POSIX_time'][:].flatten())

        try:
            UoA_data.get_projected_coords()
        except ImportError:
            pass

        UoA_data.trace_int = UoA_data.decday[1] - UoA_data.decday[0]
        UoA_data.pressure = np.zeros_like(UoA_data.decday)
        UoA_data.trig = np.zeros_like(UoA_data.decday).astype(int)
        UoA_data.trig_level = 0.
        UoA_data.flags = RadarFlags()
        UoA_data.flags.power = False
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
