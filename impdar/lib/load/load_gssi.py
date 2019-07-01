#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Load data from SIR3000 or SIR4000
"""
import codecs
import os.path
import struct
import datetime
import numpy as np
from ..gpslib import RadarGPS
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags


class DZT(RadarData):
    """Subclass for GSSI data overriding the __init__ on radardata

    Note that the init here is pretty long and reads everything from the GSSI header.
    I left this in there the in the hopes that it will be useful to somebody,
    but ImpDAR does not use all this information
    """
    header = None
    samp = None

    def __init__(self, fn_dzt):
        """Read in a DZT
        """
        super(DZT, self).__init__(None)
        gssi_header = GSSIHeader()
        with open(fn_dzt, 'rb') as fid:
            lines = fid.read()
        gssi_header.tag = struct.unpack('<H', lines[0:2])[0]
        gssi_header.data = struct.unpack('<H', lines[2:4])[0]
        gssi_header.nsamp = struct.unpack('<H', lines[4:6])[0]
        gssi_header.bits = struct.unpack('<H', lines[6:8])[0]
        gssi_header.bytes = gssi_header.bits // 8
        if gssi_header.bits == 32:
            gssi_header.us_dattype = 'I'
        elif gssi_header.bits == 16:
            gssi_header.us_dattype = 'H'

        if gssi_header.bits == 32:
            gssi_header.s_dattype = 'i'
        elif gssi_header.bits == 16:
            gssi_header.s_dattype = 'h'
        gssi_header.zero = struct.unpack('<h', lines[8:10])[0]
        gssi_header.sps = struct.unpack('<f', lines[10:14])[0]
        gssi_header.spm = struct.unpack('<f', lines[14:18])[0]
        gssi_header.mpm = struct.unpack('<f', lines[18:22])[0]
        gssi_header.position = struct.unpack('<f', lines[22:26])[0]
        gssi_header.range = struct.unpack('<f', lines[26:30])[0]

        gssi_header.npass = struct.unpack('<h', lines[30:32])[0]

        create_full = struct.unpack('<4s', lines[32:36])[0]
        modify_full = struct.unpack('<4s', lines[36:40])[0]

        gssi_header.Create = _to_date(create_full)
        gssi_header.Modify = _to_date(modify_full)

        gssi_header.rgain = struct.unpack('<H', lines[40:42])[0]
        gssi_header.nrgain = struct.unpack('<H', lines[42:44])[0] + 2
        gssi_header.text = struct.unpack('<H', lines[44:46])[0]
        gssi_header.ntext = struct.unpack('<H', lines[46:48])[0]
        gssi_header.proc = struct.unpack('<H', lines[48:50])[0]
        gssi_header.nproc = struct.unpack('<H', lines[50:52])[0]
        gssi_header.nchan = struct.unpack('<H', lines[52:54])[0]

        gssi_header.epsr = struct.unpack('<f', lines[54:58])[0]
        gssi_header.top = struct.unpack('<f', lines[58:62])[0]
        gssi_header.depth = struct.unpack('<f', lines[62:66])[0]

        gssi_header.reserved = struct.unpack('<31c', lines[66:97])
        gssi_header.dtype = struct.unpack('<c', lines[97:98])[0]
        gssi_header.antname = struct.unpack('<14c', lines[98:112])

        gssi_header.chanmask = struct.unpack('<H', lines[112:114])[0]
        gssi_header.name = struct.unpack('<12c', lines[114:126])
        gssi_header.chksum = struct.unpack('<H', lines[126:128])[0]

        gssi_header.breaks = struct.unpack('<H', lines[gssi_header.rgain:gssi_header.rgain + 2])[0]
        gssi_header.Gainpoints = np.array(struct.unpack('<{:d}i'.format(gssi_header.nrgain), lines[gssi_header.rgain + 2:gssi_header.rgain + 2 + 4 * (gssi_header.nrgain)]))
        gssi_header.Gain = 0
        if gssi_header.ntext != 0:
            gssi_header.comments = struct.unpack('<{:d}s'.format(gssi_header.ntext), lines[130 + 2 * gssi_header.Gain: 130 + gssi_header.bytes * gssi_header.Gain + gssi_header.ntext])[0]
        else:
            gssi_header.comments = ''
        if gssi_header.nproc != 0:
            gssi_header.proccessing = struct.unpack('<{:d}s'.format(gssi_header.nproc), lines[130 + gssi_header.bytes * gssi_header.Gain + gssi_header.ntext:130 + gssi_header.bytes * gssi_header.Gain + gssi_header.ntext + gssi_header.nproc])[0]
        else:
            gssi_header.proc = ''
        d = np.array(struct.unpack('<{:d}'.format((len(lines) - 36 * 4096) // gssi_header.bytes) + gssi_header.us_dattype, lines[36 * 4096:]))
        d = d.reshape((gssi_header.nsamp, -1), order='F')
        d[0, :] = d[2, :]
        d[1, :] = d[2, :]
        d = d + gssi_header.zero

        # legacy from when this was pygssi
        self.gssi_header = gssi_header

        # relevant variables for impdar
        self.data = d
        self.chan = self.gssi_header.nchan
        self.snum = self.gssi_header.nsamp
        self.tnum = self.data.shape[1]
        self.trace_num = np.arange(self.data.shape[1]) + 1
        self.trig_level = np.zeros((self.tnum, ))
        self.pressure = np.zeros((self.tnum, ))
        self.flags = RadarFlags()
        self.dt = self.gssi_header.range / self.gssi_header.nsamp * 1.0e-9
        self.travel_time = np.atleast_2d(np.arange(0, self.gssi_header.range / 1.0e3, self.dt * 1.0e6)).transpose() + self.dt * 1.0e6
        self.trig = self.gssi_header.zero

        # Now deal with the gps info
        if os.path.exists(os.path.splitext(fn_dzt)[0] + '.DZG'):
            self.gps_data = _get_dzg_data(os.path.splitext(fn_dzt)[0] + '.DZG', self.trace_num)
            self.lat = self.gps_data.lat
            self.long = self.gps_data.lon
            self.x_coord = self.gps_data.x
            self.y_coord = self.gps_data.y
            self.dist = self.gps_data.dist.flatten()
            self.elev = self.gps_data.z

            timezero = datetime.datetime(1970, 1, 1, 0, 0, 0)
            day_offset = self.gssi_header.Create - timezero
            tmin, tmax = day_offset.days + np.min(self.gps_data.dectime), day_offset.days + np.max(self.gps_data.dectime)
            self.decday = np.linspace(tmin, tmax, self.tnum)
            self.trace_int = np.hstack((np.array(np.nanmean(np.diff(self.dist))), np.diff(self.dist)))

        else:
            self.lat = np.zeros((self.data.shape[1],))
            self.long = np.zeros((self.data.shape[1],))
            self.x_coord = np.zeros((self.data.shape[1],))
            self.y_coord = np.zeros((self.data.shape[1],))
            self.dist = np.zeros((self.data.shape[1],))
            self.elev = np.zeros((self.data.shape[1],))
            self.decday = np.arange(self.data.shape[1])
            self.trace_int = np.ones((self.data.shape[1],))

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


class GSSIHeader:
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


def _to_date(bin, le=True):
    def _bit_to_int(bits):
        return sum([(2 ** i) * bit for i, bit in enumerate(bits)])

    def _bits(bytes):
        for b in bytes:
            for i in range(8):
                yield (int(b) >> i) & 1

    a = _time()
    try:
        bit = [b for b in _bits(bin)]
    except ValueError:
        # No guarantees on this one. I don't think that python2 likes my _bits function.
        return datetime.datetime(2000, 1, 1, 1, 1, 1)
    a.sec2 = _bit_to_int(bit[0:5])
    a.minute = _bit_to_int(bit[5:11])
    a.hour = _bit_to_int(bit[11:16])
    a.day = _bit_to_int(bit[16:21])
    a.month = _bit_to_int(bit[21:25])
    a.year = _bit_to_int(bit[25:32])
    if a.year > 0:
        return datetime.datetime(a.year, a.month, a.day, a.hour, a.minute, a.sec2)
    else:
        return datetime.datetime(2000, 1, 1, 1, 1, 1)


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

    with codecs.open(fn, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    # We have to be careful with this to permit other NMEA strings to have been recorded and to be sure that the indices line up
    gssis_inds = [i for i, line in enumerate(lines) if 'GSSIS' in line]
    gga_inds = [i for i, line in enumerate(lines) if 'GGA' in line]
    # we may have some records without GGA, so check if this is the case; we keep track of the offset if so
    gssis_inds_keep = []
    offset_ind = 0
    for i, j in enumerate(gssis_inds[:-1]):
        if (gga_inds[i + offset_ind] > j and gga_inds[i + offset_ind] < gssis_inds[i + 1]):
            gssis_inds_keep.append(j)
        else:
            offset_ind -= 1
    if gga_inds[-1] > gssis_inds[-1]:
        gssis_inds_keep.append(gssis_inds[-1])

    scans = np.array(list(map(lambda x: int(x.split(',')[1]), [line for i, line in enumerate(lines) if i in gssis_inds_keep])))
    data = RadarGPS([line for i, line in enumerate(lines) if i in gga_inds], scans, trace_nums)
    return data


def load_gssi(fn, *args, **kwargs):
    return DZT(fn)
