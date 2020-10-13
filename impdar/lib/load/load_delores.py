#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""
Load files from the BAS system DeLoRES (Deep Looking Radio Echo Sounder).

These are saved as hdf5 files with two channels (A and B).
For consistency with the UW HF system I am using channel names 1 and 2.

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Oct 9 2019

"""

import numpy as np
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags
from ..gpslib import RadarGPS
import datetime

try:
    import h5py
    H5 = True
except ImportError:
    H5 = False


def _get_gps_data(gga, ggis, trace_nums):
    """Read GPS data from NMEA strings.

    Parameters
    ----------
    gga: list of strs
        the gga lines
    ggis: list of strs
        the ggis lines
    trace_nums: list
        The traces that these are associated with

    Returns
    -------
    data: :class:`~impdar.lib.gpslib.nmea_info`
    """
    scans = np.array(list(map(lambda x: int(
        float(x.rstrip('\n\r ').split(' ')[-1])), ggis)))
    data = RadarGPS(gga, scans, trace_nums)
    return data


def load_delores(fn_del, channel=1, *args, **kwargs):
    """Load a DeLoRES file into ImpDAR."""
    if not H5:
        raise ImportError('You need H5 to load gprMax')

    del_data = RadarData(None)
    del_data.fn = fn_del

    # Pull from the file
    with h5py.File(fn_del) as f_in:
        # Get the specific channel
        if channel == 1:
            h5_ch = f_in['Channel_A']
            del_data.chan = 1
        if channel == 2:
            h5_ch = f_in['Channel_B']
            del_data.chan = 2
        # Get the time step
        del_data.dt = h5_ch.attrs['SampleRate'] * 1e-9
        # Get all the traces from this channel
        tr_names = list(h5_ch.keys())
        # Create empty data arrays
        del_data.tnum = len(h5_ch.keys())
        del_data.snum = h5_ch.attrs['NoOfSamples']
        del_data.data = np.empty((del_data.snum, del_data.tnum))
        nmea_strs = np.empty((del_data.tnum)).astype(str)
        decday = np.empty((del_data.tnum))
        # Fill the data arrays
        for i, tr in enumerate(tr_names):
            # Get the trace
            del_data.data[:, i] = h5_ch[tr]
            # Get NMEA strings for each trace
            nmea_strs[i] = h5_ch[tr].attrs['GGA']
            nmea_strs[i] = nmea_strs[i][2:]
            # Get GPS time for each trace
            if not hasattr(h5_ch[tr].attrs['Time'], '__len__'):
                decday[i] = np.nan
            else:
                time_array = h5_ch[tr].attrs['Time']
                date = datetime.date(time_array[0],
                                     time_array[1],
                                     time_array[2])
                time = time_array[3] + (time_array[4] + (
                    time_array[5] + time_array[6] / 1.0e6) / 60.) / 60.
                decday[i] = 366. + datetime.date.toordinal(date) + time / 24.

    """
    Still looking for more info from the delores files
    # GPS information
    del_data.gps_data = _get_gps_data(nmea_strs,gps_times)
    del_data.lat = del_data.gps_data.lat
    del_data.long = del_data.gps_data.lon
    del_data.x_coord = del_data.gps_data.x
    del_data.y_coord = del_data.gps_data.y
    del_data.dist = del_data.gps_data.dist.flatten()
    del_data.elev = del_data.gps_data.z
    del_data.trace_int = np.hstack((np.array(np.nanmean(
        np.diff(del_data.dist))), np.diff(del_data.dist)))
    """

    # all zeros for now
    del_data.lat = np.zeros((del_data.tnum,))
    del_data.long = np.zeros((del_data.tnum,))
    del_data.x_coord = np.zeros((del_data.tnum,))
    del_data.y_coord = np.zeros((del_data.tnum,))
    del_data.dist = np.zeros((del_data.tnum,))
    del_data.elev = np.zeros((del_data.tnum,))
    del_data.trace_int = np.zeros((del_data.tnum,))

    # pretrigger
    del_data.trig = np.zeros((del_data.tnum,))
    del_data.trig_level = np.zeros((del_data.tnum,))

    # other variables are from the array shape
    del_data.decday = decday
    del_data.trace_num = np.arange(del_data.data.shape[1]) + 1
    del_data.pressure = np.zeros((del_data.tnum,))
    del_data.flags = RadarFlags()
    del_data.travel_time = del_data.dt * 1.0e6 * np.arange(del_data.snum)

    del_data.check_attrs()
    return del_data
