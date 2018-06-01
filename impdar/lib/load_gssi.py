#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3 license.

"""

"""
import os.path
import struct
import numpy as np
from .gpslib import RadarGPS
from .RadarData import RadarData, RadarFlags
import datetime


class DZT(RadarData):
    header = None
    samp = None

    def __init__(self, fn):
        rh = RH()
        with open(fn, 'rb') as fid:
            lines = fid.read()
        rh.tag = struct.unpack('<H', lines[0:2])[0]
        rh.data = struct.unpack('<H', lines[2:4])[0]
        rh.nsamp = struct.unpack('<H', lines[4:6])[0]
        rh.bits = struct.unpack('<H', lines[6:8])[0]
        rh.bytes = rh.bits // 8
        if rh.bits == 32:
            rh.us_dattype = 'I'
        elif rh.bits == 16:
            rh.us_dattype = 'H'
        if rh.bits == 32:
            rh.s_dattype = 'i'
        elif rh.bits == 16:
            rh.s_dattype = 'h'
        rh.zero = struct.unpack('<h', lines[8:10])[0]
        rh.sps = struct.unpack('<f', lines[10:14])[0]
        rh.spm = struct.unpack('<f', lines[14:18])[0]
        rh.mpm = struct.unpack('<f', lines[18:22])[0]
        rh.position = struct.unpack('<f', lines[22:26])[0]
        rh.range = struct.unpack('<f', lines[26:30])[0]

        rh.npass = struct.unpack('<h', lines[30:32])[0]

        create_full = struct.unpack('<4s', lines[32:36])[0]
        modify_full = struct.unpack('<4s', lines[36:40])[0]

        rh.Create = _to_date(create_full)
        rh.Modify = _to_date(modify_full)

        rh.rgain = struct.unpack('<H', lines[40:42])[0]
        rh.nrgain = struct.unpack('<H', lines[42:44])[0] + 2
        rh.text = struct.unpack('<H', lines[44:46])[0]
        rh.ntext = struct.unpack('<H', lines[46:48])[0]
        rh.proc = struct.unpack('<H', lines[48:50])[0]
        rh.nproc = struct.unpack('<H', lines[50:52])[0]
        rh.nchan = struct.unpack('<H', lines[52:54])[0]

        rh.epsr = struct.unpack('<f', lines[54:58])[0]
        rh.top = struct.unpack('<f', lines[58:62])[0]
        rh.depth = struct.unpack('<f', lines[62:66])[0]

        rh.reserved = struct.unpack('<31c', lines[66:97])
        rh.dtype = struct.unpack('<c', lines[97:98])[0]
        rh.antname = struct.unpack('<14c', lines[98:112])

        rh.chanmask = struct.unpack('<H', lines[112:114])[0]
        rh.name = struct.unpack('<12c', lines[114:126])
        rh.chksum = struct.unpack('<H', lines[126:128])[0]

        rh.breaks = struct.unpack('<H', lines[rh.rgain:rh.rgain + 2])[0]
        rh.Gainpoints = np.array(struct.unpack('<{:d}i'.format(rh.nrgain), lines[rh.rgain + 2:rh.rgain + 2 + 4 * (rh.nrgain)]))
        rh.Gain = 0
        if rh.ntext != 0:
            rh.comments = struct.unpack('<{:d}s'.format(rh.ntext), lines[130 + 2 * rh.Gain: 130 + rh.bytes * rh.Gain + rh.ntext])[0]
        else:
            rh.comments = ''
        if rh.nproc != 0:
            rh.proccessing = struct.unpack('<{:d}s'.format(rh.nproc), lines[130 + rh.bytes * rh.Gain + rh.ntext:130 + rh.bytes * rh.Gain + rh.ntext + rh.nproc])[0]
        else:
            rh.proc = ''
        d = np.array(struct.unpack('<{:d}'.format((len(lines) - 36 * 4096) // rh.bytes) + rh.us_dattype, lines[36 * 4096:]))
        d = d.reshape((rh.nsamp, -1), order='F')
        d[0, :] = d[2, :]
        d[1, :] = d[2, :]
        d = d + rh.zero

        # legacy from when this was pygssi
        self.rh = rh

        # relevant variables for impdar
        self.data = d
        self.chan = self.rh.nchan
        self.snum = self.rh.nsamp
        self.tnum = self.data.shape[1]
        self.trace_num = np.arange(self.data.shape[1]) + 1
        self.trig_level = np.zeros((self.tnum, ))
        self.pressure = np.zeros((self.tnum, ))
        self.flags = RadarFlags()
        self.dt = self.rh.range / self.rh.nsamp * 1.0e-9
        self.travel_time = np.atleast_2d(np.arange(0, self.rh.range / 1.0e3, self.dt * 1.0e6)).transpose() + self.dt * 1.0e6
        self.trig = self.rh.zero

        # Now deal with the gps info
        self.gps_data = _get_dzg_data(os.path.splitext(fn)[0] + '.DZG', self.trace_num)
        self.lat = self.gps_data.lat
        self.long = self.gps_data.lon
        self.x_coord = self.gps_data.x
        self.y_coord = self.gps_data.y
        self.dist = self.gps_data.dist.flatten()
        self.elev = self.gps_data.z

        timezero = datetime.datetime(2017, 1, 1, 0, 0, 0)
        day_offset = self.rh.Create - timezero
        tmin, tmax = day_offset.days + np.min(self.gps_data.dectime), day_offset.days + np.max(self.gps_data.dectime)
        self.decday = np.linspace(tmin, tmax, self.tnum)
        self.trace_int = np.hstack((np.array(np.nanmean(np.diff(self.dist))), np.diff(self.dist)))

        for attr in ['chan', 'data', 'decday', 'dist', 'dt', 'elev', 'flags', 'lat', 'long', 'pressure', 'snum', 'tnum', 'trace_int', 'trace_num', 'travel_time', 'trig', 'trig_level', 'x_coord', 'y_coord']:
            if getattr(self, attr) is None:
                print(attr + ' is not defined')
                setattr(self, attr, 0)


class _time:
    sec2 = None
    minute = None
    hour = None
    day = None
    month = None
    year = None


class RH:
    tag = None
    data = None
    nsamp = None
    bits = None
    bytes = None
    us_dattype = None
    s_dattype = None
    rgain = None
    nrgain = None
    checksum = None
    antname = None

    def __str__(self):
        return 'rgain: {:d}, nrgain {:d}'.format(self.rgain, self.nrgain)

    def __repr__(self):
        return self.__str__()


def _to_date(bin, le=True):
    def _bit_to_int(bits):
        return sum([(2 ** i) * bit for i, bit in enumerate(bits)])

    def _bits(bytes):
        for b in bytes:
            for i in range(8):
                yield (b >> i) & 1

    a = _time()
    bit = [b for b in _bits(bin)]
    a.sec2 = _bit_to_int(bit[0:5])
    a.minute = _bit_to_int(bit[5:11])
    a.hour = _bit_to_int(bit[11:16])
    a.day = _bit_to_int(bit[16:21])
    a.month = _bit_to_int(bit[21:25])
    a.year = _bit_to_int(bit[25:32])
    if a.year > 0:
        return datetime.datetime(a.year, a.month, a.day, a.hour, a.minute, a.sec2)
    else:
        return None


def _get_dzg_data(fn, trace_nums):
    """Read GPS data associated with a GSSI sir4000 file.

    Parameters
    ----------
    fn: str
        A dzg file with ggis and gga strings.
    trace_nums: np.ndarray
        The traces on which to interpolate the (sparse) GPS input
        
    
    Returns
    -------
    data: :class:`~impdar.lib.gpslib.nmea_info`
    """

    with open(fn) as f:
        lines = f.readlines()
    ggis = lines[::3]
    gga = lines[1::3]
    scans = np.array(list(map(lambda x: int(x.split(',')[1]), ggis)))
    data = RadarGPS(gga, scans, trace_nums)
    return data


def load_gssi(fn, *args, **kwargs):
    return DZT(fn)
