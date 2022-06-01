#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
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


class GSSITime:
    """GSSI uses a weird date format that we probably want to read"""
    sec2 = None
    minute = None
    hour = None
    day = None
    month = None
    year = None

    def __init__(self, binary_data, le=True):
        """Read in GSSI binary date data"""
        try:
            bits = [b for b in _bits(binary_data)]
            self.from_bits(bits)
        except ValueError:
            # No guarantees on this one. I don't think that python2 likes my _bits function.
            pass

    def from_bits(self, bits):
        """Define the date from the bits input"""
        self.sec2 = _bit_to_int(bits[0:5])
        self.minute = _bit_to_int(bits[5:11])
        self.hour = _bit_to_int(bits[11:16])
        self.day = _bit_to_int(bits[16:21])
        self.month = _bit_to_int(bits[21:25])
        self.year = _bit_to_int(bits[25:32])

    def to_datetime(self):
        """Convert from the GSSI data to a datetime.datetime output"""
        if self.year > 0:
            # Not 100% sure that 1980 is the correct offset here
            return datetime.datetime(self.year + 1980,
                                     self.month,
                                     self.day,
                                     self.hour,
                                     self.minute,
                                     self.sec2)
        return datetime.datetime(2000, 1, 1, 1, 1, 1)


def _bit_to_int(bits):
    return sum([(2 ** i) * bit for i, bit in enumerate(bits)])


def _bits(bytes_in):
    for b_i in bytes_in:
        for j in range(8):
            yield (int(b_i) >> j) & 1


def _get_dzg_data(fn_dzg, trace_nums):
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

    with codecs.open(fn_dzg, 'r', encoding='utf-8', errors='ignore') as f_in:
        lines = f_in.readlines()
    # We have to be careful with this to permit other NMEA strings to have been recorded
    # and to be sure that the indices line up
    all_gga_inds = [i for i, line in enumerate(lines) if '$GPGGA' == line.split(',')[0]]

    # Get the corresponding GSSI trace numbers
    all_gssis_inds = np.array([i for i, line in enumerate(lines) if line.split(',')[0] == '$GSSIS'])
    gssis_inds = []
    gga_inds = []
    for i, lineind in enumerate(all_gga_inds):
        if i == 0:
            prevind = 0
        else:
            prevind = all_gga_inds[i - 1]
        rel_inds = all_gssis_inds[np.logical_and(all_gssis_inds < lineind, all_gssis_inds > prevind)]
        if len(rel_inds) > 0:
            try:
                # we can still have bad GSSI strings
                if float(lines[np.max(rel_inds)].split(',')[1]).is_integer():
                    gssis_inds.append(np.max(rel_inds))
                    gga_inds.append(lineind)
            except ValueError:
                continue

    # we may have some records without GGA, so check if this is the case;
    # we keep track of the offset if so
    gssis_inds_keep = []
    offset_ind = 0
    for i, j in enumerate(gssis_inds[:-1]):
        if (gga_inds[i + offset_ind] > j and gga_inds[i + offset_ind] < gssis_inds[i + 1]):
            gssis_inds_keep.append(j)
        else:
            offset_ind -= 1
    if gga_inds[-1] > gssis_inds[-1]:
        gssis_inds_keep.append(gssis_inds[-1])

    scans = np.array(list(map(lambda x: int(x.split(',')[1]),
                              [line for i, line in enumerate(lines) if i in gssis_inds_keep])))
    data = RadarGPS([line for i, line in enumerate(lines) if i in gga_inds], scans, trace_nums)
    return data


def load_gssi(fn_dzt, *args, **kwargs):
    """Return a RadarData object with the information from a gssi file

    This reader is has many commented-out liens to read everything from the GSSI header.
    I left this in there the in the hopes that it will be useful to somebody,
    but ImpDAR does not use all this information
    """
    dzt_data = RadarData(None)
    dzt_data.fn = fn_dzt
    with open(fn_dzt, 'rb') as fid:
        lines = fid.read()
    # tag = struct.unpack('<H', lines[0:2])[0]
    # data = struct.unpack('<H', lines[2:4])[0]
    dzt_data.snum = struct.unpack('<H', lines[4:6])[0]
    bits = struct.unpack('<H', lines[6:8])[0]
    n_bytes = bits // 8
    if bits == 32:
        us_dattype = 'I'
        np_dtype = np.int32
    elif bits == 16:
        us_dattype = 'H'
        np_dtype = np.int16
    # if bits == 32:
    #     s_dattype = 'i'
    # elif bits == 16:
    #     s_dattype = 'h'
    # dzt_data.trig = struct.unpack('<h', lines[8:10])[0] * np.ones(dzt_data.tnum)
    # sps = struct.unpack('<f', lines[10:14])[0]
    # spm = struct.unpack('<f', lines[14:18])[0]
    # mpm = struct.unpack('<f', lines[18:22])[0]
    # position = struct.unpack('<f', lines[22:26])[0]
    dzt_data.range = struct.unpack('<f', lines[26:30])[0]

    # npass = struct.unpack('<h', lines[30:32])[0]
    create_full = GSSITime(struct.unpack('<4s', lines[32:36])[0])
    # modify_full = struct.unpack('<4s', lines[36:40])[0]
    dzt_data.create = create_full.to_datetime()
    # Modify = _to_date(modify_full)

    # rgain = struct.unpack('<H', lines[40:42])[0]
    # nrgain = struct.unpack('<H', lines[42:44])[0] + 2
    # text = struct.unpack('<H', lines[44:46])[0]
    # ntext = struct.unpack('<H', lines[46:48])[0]
    # proc = struct.unpack('<H', lines[48:50])[0]
    # nproc = struct.unpack('<H', lines[50:52])[0]
    dzt_data.chan = struct.unpack('<H', lines[52:54])[0]
    # epsr = struct.unpack('<f', lines[54:58])[0]
    # top = struct.unpack('<f', lines[58:62])[0]
    # depth = struct.unpack('<f', lines[62:66])[0]
    # reserved = struct.unpack('<31c', lines[66:97])
    # dtype = struct.unpack('<c', lines[97:98])[0]
    # antname = struct.unpack('<14c', lines[98:112])
    # chanmask = struct.unpack('<H', lines[112:114])[0]
    # name = struct.unpack('<12c', lines[114:126])
    # chksum = struct.unpack('<H', lines[126:128])[0]
    # breaks = struct.unpack('<H', lines[rgain:rgain + 2])[0]
    # Gainpoints = np.array(struct.unpack('<{:d}i'.format(nrgain),
    #                    lines[rgain + 2:rgain + 2 + 4 * (nrgain)]))
    # gain = 0
    # if ntext != 0:
    #     comments = struct.unpack('<{:d}s'.format(ntext),
    #                              lines[130 + 2 * Gain: 130 + bytes * Gain + ntext])[0]
    # else:
    #     comments = ''
    # if nproc != 0:
    #     proccessing = struct.unpack('<{:d}s'.format(nproc),
    #               lines[130 + bytes * Gain + ntext:130 + bytes * Gain + ntext + nproc])[0]
    # else:
    #     processing = ''
    try:
        header_len = 32768*n_bytes # TODO: David originally had this as 36*4096, we still need to figure out when it changes
        try:
            data = np.array(struct.unpack('<{:d}'.format((len(lines) - header_len) // n_bytes) + us_dattype, lines[header_len:]), dtype=np_dtype).reshape((dzt_data.snum, -1), order='F')
        except OverflowError:
            # This is needed on Windows for some reason--I think if OS is 32 bit?
            data = np.array(struct.unpack('<{:d}'.format((len(lines) - header_len) // n_bytes) + us_dattype, lines[header_len:]), dtype=np.int64).reshape((dzt_data.snum, -1), order='F')
    except IndexError:
        header_len = 512*n_bytes
        data = np.array(struct.unpack('<{:d}'.format((len(lines) - header_len) // n_bytes) + us_dattype, lines[header_len:]), dtype=np_dtype).reshape((dzt_data.snum, -1), order='F')
    data[0, :] = data[2, :]
    data[1, :] = data[2, :]
    # data = data + dzt_data.trig
    dzt_data.data = data

    dzt_data.tnum = dzt_data.data.shape[1]
    dzt_data.trace_num = np.arange(dzt_data.data.shape[1]) + 1
    dzt_data.trig_level = 0.
    dzt_data.trig = struct.unpack('<h', lines[8:10])[0] * np.ones(dzt_data.tnum)

    dzt_data.pressure = np.zeros((dzt_data.tnum, ))
    dzt_data.flags = RadarFlags()
    dzt_data.dt = dzt_data.range / dzt_data.snum * 1.0e-9
    dzt_data.travel_time = np.atleast_2d(np.arange(0,
                                                   dzt_data.range / 1.0e3,
                                                   dzt_data.dt * 1.0e6)).transpose()
    dzt_data.travel_time += dzt_data.dt * 1.0e6

    # Now deal with the gps info
    if os.path.exists(os.path.splitext(fn_dzt)[0] + '.DZG'):
        dzt_data.gps_data = _get_dzg_data(os.path.splitext(fn_dzt)[0] + '.DZG', dzt_data.trace_num)
        dzt_data.lat, dzt_data.long = dzt_data.gps_data.lat, dzt_data.gps_data.lon
        dzt_data.x_coord, dzt_data.y_coord = dzt_data.gps_data.x, dzt_data.gps_data.y
        dzt_data.dist = dzt_data.gps_data.dist.flatten()
        dzt_data.elev = dzt_data.gps_data.z

        timezero = datetime.datetime(1, 1, 1, 0, 0, 0)
        day_offset = dzt_data.create - timezero
        tmin = day_offset.days + np.min(dzt_data.gps_data.dectime) + 377.  # matlab compat
        tmax = day_offset.days + np.max(dzt_data.gps_data.dectime) + 377.
        dzt_data.decday = np.linspace(tmin, tmax, dzt_data.tnum)
        dzt_data.trace_int = np.hstack((np.array(np.nanmean(np.diff(dzt_data.dist))),
                                        np.diff(dzt_data.dist)))

    else:
        dzt_data.lat = np.zeros((dzt_data.data.shape[1],))
        dzt_data.long = np.zeros((dzt_data.data.shape[1],))
        dzt_data.x_coord = np.zeros((dzt_data.data.shape[1],))
        dzt_data.y_coord = np.zeros((dzt_data.data.shape[1],))
        dzt_data.dist = np.zeros((dzt_data.data.shape[1],))
        dzt_data.elev = np.zeros((dzt_data.data.shape[1],))
        dzt_data.decday = np.arange(dzt_data.data.shape[1])
        dzt_data.trace_int = np.ones((dzt_data.data.shape[1],))

    dzt_data.check_attrs()
    return dzt_data
