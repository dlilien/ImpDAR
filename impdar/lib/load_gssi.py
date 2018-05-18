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
from .gpslib import nmea_all_info


class DZT:
    header = None
    samp = None

    def __init__(self, header, sample):
        self.header = header
        self.samp = sample


class time:
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


def bits(bytes):
    for b in bytes:
        for i in range(8):
            yield (b >> i) & 1


def to_date(bin, le=True):
    a = time()
    bit = [b for b in bits(bin)]
    a.sec2 = bit_to_int(bit[0:5])
    a.minute = bit_to_int(bit[5:11])
    a.hour = bit_to_int(bit[11:16])
    a.day = bit_to_int(bit[16:21])
    a.month = bit_to_int(bit[21:25])
    a.year = bit_to_int(bit[25:32])
    return a


def bit_to_int(bits):
    return sum([(2 ** i) * bit for i, bit in enumerate(bits)])


def read_dzt(fn, rev=False):
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

    rh.Create = to_date(create_full)
    rh.Modify = to_date(modify_full)

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

    # d = np.array(struct.unpack('<{:d}'.format((len(lines) - 1024) // 2) + dattype, lines[1024:]))
    d = np.array(struct.unpack('<{:d}'.format((len(lines) - 36 * 4096) // rh.bytes) + rh.us_dattype, lines[36 * 4096:]))
    d = d.reshape((rh.nsamp, -1), order='F')
    d[0, :] = d[2, :]
    d[1, :] = d[2, :]
    d = d + rh.zero
    if rev:
        d = np.fliplr(d)
    dat = DZT(rh, d)
    return dat


def get_dzg_data(fn, rev=False):
    """Read GPS data associated with a GSSI sir4000 file.

    Parameters
    ----------
    fn: str
        A dzg file with ggis and gga strings.
    rev: bool, optional
        Reverse the points in this file (used for concatenating radar files). Default False.
    
    Returns
    -------
    data: :class:`~pygssi.lib.gpslib.nmea_info`
    """

    with open(fn) as f:
        lines = f.readlines()
    ggis = lines[::3]
    gga = lines[1::3]
    data = nmea_all_info(gga)
    data.scans = np.array(list(map(lambda x: int(x.split(',')[1]), ggis)))
    if rev:
        data.rev()
    data.get_all()
    return data


def load_gssi(fn, *args, **kwargs):
        gps_data = get_dzg_data(os.path.splitext(fn)[0] + '.DZG')

        # Now find the x coordinates for plotting
        dzt = read_dzt(fn)
        return gps_data, dzt
