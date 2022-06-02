#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""
Load an output file processed from a University of Alabama VHF, UWB, or UHF radar.

Data are all processed, but some more complex data views are not available (and some simple things,
like constant trace spacing, are not done).
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
    """Load data processed by the UoA RSC Matlab processor

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


def load_UoA_h5(fn, gps_offset=0.0):
    """Load MCoRDS data in .mat format downloaded from the CReSIS ftp client

    Parameters
    ----------
    fn: str
        The filename to load
    """
    data_list = []
    with h5py.File(fn, 'r') as fin:
        if fin.attrs['Type'] != 'MultiChannel':
            raise ValueError('Can only unpack MultiChannel UoA data')
        if 'processed' in fin:
            for name in fin['processed'].keys():
                for integrator in fin['processed'][name].keys():
                    grp = fin['processed'][name][integrator]
                    UoA_data = RadarData(None)
                    UoA_data.fn = fn[:-3] + name + '_Int' + integrator[-1]
                    UoA_data.chan = 999
                    _load_group(UoA_data, grp, gps_offset)
                    data_list.append(UoA_data)

        else:
            print('No processed data found, reading channels')
            for i in range(8):
                if f'channel_{i}' in fin:
                    for integrator in fin[f'channel_{i}'].keys():
                        grp = fin[f'channel_{i}'][integrator]
                        UoA_data = RadarData(None)
                        UoA_data.fn = fn[:-3] + f'_ch{i}_Int' + integrator[-1]

                        UoA_data.chan = i
                        _load_group(UoA_data, grp, gps_offset)

                        data_list.append(UoA_data)

    return data_list


def _load_group(UoA_data, grp, gps_offset):
    UoA_data.data = grp['Chirps'][()]
    # complex data
    if (len(UoA_data.data.dtype) == 2) or (UoA_data.data.dtype in [np.complex64, np.complex128]):
        UoA_data.data = 10.0 * np.log10(np.sqrt(np.real(UoA_data.data) ** 2.0 + np.imag(UoA_data.data) ** 2.0))
    else:
        UoA_data.data = 10.0 * np.log10(UoA_data.data)
    UoA_data.snum, UoA_data.tnum = int(UoA_data.data.shape[0]), int(UoA_data.data.shape[1])
    UoA_data.trace_num = np.arange(UoA_data.tnum) + 1
    UoA_data.travel_time = grp['_time'][()] * 1.0e6
    UoA_data.dt = np.mean(np.diff(UoA_data.travel_time)) * 1.0e-6
    if 'datetime' in grp:
        nminfo = nmea_info()
        dt = grp['datetime'][()].astype('datetime64[ms]').astype(int) / 1000.0
        nminfo.time = (dt + gps_offset) / (24. * 60. * 60.)  # Stored in seconds
        nminfo.ppstime = dt + gps_offset

        nminfo.lat = grp['lat'][:].flatten()
        nminfo.lon = grp['lon'][:].flatten()
        nminfo.elev = np.zeros_like(nminfo.lat)

        UoA_data.lat = interp1d(nminfo.ppstime, nminfo.lat, fill_value='extrapolate')(dt)
        UoA_data.long = interp1d(nminfo.ppstime, nminfo.lon, fill_value='extrapolate')(dt)
        UoA_data.elev = np.zeros_like(UoA_data.lat)
        UoA_data.elev[:] = np.nan
        UoA_data.decday = interp1d(nminfo.ppstime, nminfo.time, fill_value='extrapolate')(dt)

        if 'x' in grp:
            UoA_data.x_coord = grp['x'][()]
            UoA_data.y_coord = grp['y'][()]
        else:
            try:
                UoA_data.get_projected_coords()
            except ImportError:
                pass
    else:
        print('WARNING: datetime information missing--hopefully this is loopback data???')
        UoA_data.lat = np.zeros((UoA_data.tnum,)) * np.nan
        UoA_data.long = np.zeros((UoA_data.tnum,)) * np.nan
        UoA_data.elev = np.zeros((UoA_data.tnum,)) * np.nan
        UoA_data.decday = np.zeros((UoA_data.tnum,))

    UoA_data.trace_int = UoA_data.decday[1] - UoA_data.decday[0]
    UoA_data.pressure = np.zeros_like(UoA_data.decday)
    UoA_data.trig = np.zeros_like(UoA_data.decday).astype(int)
    UoA_data.trig_level = 0.
    UoA_data.flags = RadarFlags()
    UoA_data.flags.power = False
    UoA_data.check_attrs()
