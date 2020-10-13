#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Load RAMAC data
"""

import os
import struct
import datetime
import numpy as np
from scipy.interpolate import interp1d
from ..RadarData import RadarData
from ..gpslib import nmea_info, conversions_enabled


def load_ramac(ramac_fn):
    """Return a impdar.lib.RadarData.RadarData object loaded from a Ramac file.

    Parameters
    ----------
    ramac_fn: str
        The file to load. Can include the '.rad' or not.
    """
    ramac_data = RadarData(None)

    # We can take headers, base filenames, or data files
    if len(ramac_fn) <= 4:
        header_fn = ramac_fn + '.rad'
        data_fn = ramac_fn + '.rd3'
        gps_fn = ramac_fn + '.cor'
    elif ramac_fn[-4:] == '.rd3':
        header_fn = ramac_fn[:-3] + 'rad'
        data_fn = ramac_fn[:-3] + 'rd3'
        gps_fn = ramac_fn[:-3] + 'cor'
    elif ramac_fn[-4:] != '.rad':
        header_fn = ramac_fn + '.rad'
        data_fn = ramac_fn + '.rd3'
        gps_fn = ramac_fn + '.cor'
    else:
        header_fn = ramac_fn
        data_fn = ramac_fn[:-3] + 'rd3'
        gps_fn = ramac_fn[:-3] + 'cor'

    ramac_data.fn = data_fn

    with open(header_fn) as f_header:
        header = f_header.readlines()

    ramac_data.chan = ramac_fn[-5]
    ramac_data.snum = int(header[0].rstrip('\n')[8:])
    sampling_freq = float(header[1].rstrip('\n')[10:])
    ramac_data.dt = (1. / sampling_freq) * 1.0e-6
    ramac_data.travel_time = ramac_data.dt * np.arange(ramac_data.snum) * 1.0e6
    # d_x = int(header[2].rstrip('\n')[16:])
    # if d_x == 0:
    #     d_x = 1
    ramac_data.tnum = int(header[22].rstrip('\n')[11:])
    ramac_data.trace_num = np.arange(ramac_data.tnum) + 1
    # FIX TRACE INT
    ramac_data.trace_int = float(header[9].rstrip('\n')[14:]) * np.ones(
        (ramac_data.tnum,))
    ramac_data.trig = np.ones((ramac_data.tnum,)) * 36
    ramac_data.trig_level = 0

    if os.path.exists(gps_fn):
        cor = np.genfromtxt(gps_fn, dtype=[('trace_num', int),
                                           ('date', 'S10'),
                                           ('time', 'S8'),
                                           ('lat', float),
                                           ('north', 'S1'),
                                           ('lon', float),
                                           ('east', 'S1'),
                                           ('elev', float),
                                           ('el_unit', 'S1'),
                                           ('pdop', float)])
        datetimes = np.array([d + b'T' + t for d, t in zip(
            cor['date'], cor['time'])], dtype=np.datetime64)
        decdays = datetimes - np.array(datetime.datetime(1, 1, 1, 0, 0, 0),
                                       dtype=np.datetime64)
        cor['lat'][cor['north'] != b'N'] = -1 * cor[
            'lat'][cor['north'] != b'N']
        cor['lon'][cor['east'] != b'E'] = -1 * cor['lon'][cor['east'] != b'E']

        ramac_data.decday = interp1d(cor['trace_num'],
                                     decdays,
                                     fill_value='extrapolate')(
            ramac_data.trace_num) / (24. * 60. * 60.)
        ramac_data.lat = interp1d(cor['trace_num'],
                                  cor['lat'],
                                  fill_value='extrapolate')(
            ramac_data.trace_num)
        ramac_data.long = interp1d(cor['trace_num'],
                                   cor['lon'],
                                   fill_value='extrapolate')(
            ramac_data.trace_num)
        ramac_data.elev = interp1d(cor['trace_num'],
                                   cor['elev'],
                                   fill_value='extrapolate')(
            ramac_data.trace_num)

        nminfo = nmea_info()
        nminfo.time = ramac_data.decday
        nminfo.lat = ramac_data.lat
        nminfo.lon = ramac_data.long
        nminfo.elev = ramac_data.elev
        if conversions_enabled:
            nminfo.get_utm()
            nminfo.get_dist()
            ramac_data.x_coord = nminfo.x
            ramac_data.y_coord = nminfo.y
            ramac_data.dist = nminfo.dist
        else:
            ramac_data.x_coord = ramac_data.long
            ramac_data.y_coord = ramac_data.lat
            ramac_data.dist = np.sqrt(
                ramac_data.x_coord ** 2.0 + ramac_data.y_coord ** 2.0) / 1000.0

    else:
        ramac_data.decday = np.arange(ramac_data.tnum)
        ramac_data.lat = np.arange(ramac_data.tnum)
        ramac_data.long = np.arange(ramac_data.tnum)
        ramac_data.dist = np.arange(ramac_data.tnum)
        ramac_data.elev = np.arange(ramac_data.tnum)
    ramac_data.pressure = np.zeros_like(ramac_data.dist)

    with open(data_fn, 'rb') as f_data:
        ramac_data.data = np.array(
            struct.unpack('<{:d}h'.format(
                ramac_data.tnum * ramac_data.snum), f_data.read()[:]),
            dtype=np.int16).reshape((ramac_data.snum, ramac_data.tnum),
                                    order='F')
    ramac_data.check_attrs()
    return ramac_data
