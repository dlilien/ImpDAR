#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Read Pulse Ekko data.

# Knut Christianson 6 April 2017; 20 May 2017
# Pythonized by David Lilien, May 2018
# Distributed under the GNU GPL3 license
# Benjamin Hills Added capability for v.1.5.340; Sept 27 2019
"""

import os.path
import struct
import datetime
import numpy as np

from ..gpslib import RadarGPS
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags


class TraceHeaders:
    """Class used internally to handle pulse-ekko headers."""

    def __init__(self, tnum):
        """Create a container for all the trace headers."""
        self.header_index = 0
        self.trace_numbers = np.zeros((1, tnum))
        self.positions = np.zeros((1, tnum))
        self.points_per_trace = np.zeros((1, tnum))
        self.topography = np.zeros((1, tnum))
        self.bytes_per_point = np.zeros((1, tnum))
        self.n_stackes = np.zeros((1, tnum))
        self.time_window = np.zeros((1, tnum))
        self.pos = np.zeros((3, tnum))
        self.receive = np.zeros((3, tnum))
        self.transmit = np.zeros((3, tnum))
        self.tz_adjustment = np.zeros((1, tnum))
        self.zero_flag = np.zeros((1, tnum))
        self.time_of_day = np.zeros((1, tnum))
        self.comment_flag = np.zeros((1, tnum))
        self.comment = ['' for i in range(tnum)]

    def get_header(self, offset, f_lines):
        """Get the header information for a single trace."""
        header = struct.unpack('<25f', f_lines[offset: offset + 25 * 4])
        comment = struct.unpack('<28c',
                                f_lines[offset + 25 * 4: offset + 25 * 4 + 28])
        self.trace_numbers[0, self.header_index] = header[0]
        self.positions[0, self.header_index] = header[1]
        self.points_per_trace[0, self.header_index] = header[2]
        self.topography[0, self.header_index] = header[3]
        self.bytes_per_point[0, self.header_index] = header[5]
        self.n_stackes[0, self.header_index] = header[7]
        self.time_window[0, self.header_index] = header[8]
        self.pos[0, self.header_index] = header[9]
        self.pos[1, self.header_index] = header[11]
        self.pos[2, self.header_index] = header[13]
        self.receive[0, self.header_index] = header[14]
        self.receive[1, self.header_index] = header[15]
        self.receive[2, self.header_index] = header[16]
        self.transmit[0, self.header_index] = header[17]
        self.transmit[1, self.header_index] = header[18]
        self.transmit[2, self.header_index] = header[19]
        self.tz_adjustment[0, self.header_index] = header[20]
        self.zero_flag[0, self.header_index] = header[21]
        self.time_of_day[0, self.header_index] = header[23]
        self.comment_flag[0, self.header_index] = header[24]
        self.comment[self.header_index] = str(comment[0])
        self.header_index += 1


def _get_gps_data(fn_gps, trace_nums):
    """Read GPS data associated with a Pulse Ekko .GPS file.

    Parameters
    ----------
    fn_gps: str
        A dzg file with ggis and gga strings.

    Returns
    -------
    data: :class:`~impdar.lib.gpslib.nmea_info`
    """
    with open(fn_gps) as f_in:
        lines = f_in.readlines()
    ggis = []
    gga = []
    for line in lines:
        if line[:5] == 'Trace':
            ggis.append(line)
        elif line[:6] == '$GPGGA':
            gga.append(line)
        else:
            continue
    if len(gga) == 0:
        raise ValueError('I can only do gga sentences right now')
    scans = np.array(list(map(lambda x: int(float(
        x.rstrip('\n\r ').split(' ')[-1])), ggis)))
    data = RadarGPS(gga, scans, trace_nums)
    return data


def partition_project_file(fn_project):
    """Separate profiles.

    The new pulse ekko dvl writes 'project' files with all the profiles stored
    together. We want to break them out into all the .HD header files and .DT1
    data files.

    Parameters
    ----------
    fn_project: str
        Filename for the .gpz project file
    """
    with open(fn_project, 'rb') as fin:
        f = fin.read()

    profile_num = 1
    while f.find(b'line%d' % profile_num) != -1:
        # Get the header file
        hd_start = f.find(b'line%d.hd' % (profile_num))
        hd_end = f[hd_start:].find(b'PK') + hd_start
        hd_str = str(f[hd_start:hd_end])
        hd_lines = hd_str.split('\\r\\n')
        hd_lines[0] = hd_lines[0][2:]
        hd_lines[-1] = ''

        # Get the 'ini' file
        ini_start = f.find(b'line%d.ini' % (profile_num))
        ini_end = f[ini_start:].find(b'PK') + ini_start
        ini_str = str(f[ini_start:ini_end])
        for i, line in enumerate(ini_str.split('\\r\\n')):
            if i == 0:
                hd_lines.append(line[2:len('line%d.ini' % (profile_num)) + 2])
                hd_lines.append(line[len('line%d.ini' % (profile_num)) + 2:])
            elif i == len(ini_str.split('\\r\\n')) - 1:
                continue
            else:
                hd_lines.append(line)

        # Write to the header file
        with open('LINE' + str(profile_num) + '.HD', 'w') as fout:
            for line in hd_lines:
                fout.write(line + '\n')

        # Get the data file
        dt_start = f.find(b'line%d.dt1' % (profile_num))
        dt_start += len(b'line%d.dt1' % (profile_num))
        dt_end = f[dt_start:].find(b'Lineset') + dt_start
        dt_str = f[dt_start:dt_end]
        # Write to the data file
        with open('LINE' + str(profile_num) + '.DT1', 'wb') as fout:
            fout.write(dt_str)

        profile_num += 1


def load_pe(fn_dt1, *args, **kwargs):
    """Load data from a pulse_ekko file."""
    pe_data = RadarData(None)
    pe_data.fn = fn_dt1
    bn_pe = os.path.splitext(fn_dt1)[0]
    hdname = bn_pe + '.HD'
    true_fn = bn_pe + '.DT1'
    gps_fn = bn_pe + '.GPS'

    try:
        strtypes = (unicode, str)
        openmode_unicode = 'rU'
    except NameError:
        strtypes = (str, )
        openmode_unicode = 'r'

    with open(hdname, openmode_unicode) as fin:
        if fin.read().find('1.5.340') != -1:
            pe_data.version = '1.5.340'
        else:
            pe_data.version = '1.0'
        fin.seek(0)
        for i, line in enumerate(fin):
            if 'TRACES' in line or 'NUMBER OF TRACES' in line:
                pe_data.tnum = int(line.rstrip('\n\r ').split(' ')[-1])
            if 'PTS' in line or 'NUMBER OF PTS/TRC' in line:
                pe_data.snum = int(line.rstrip('\n\r ').split(' ')[-1])
            if ('WINDOW' in line and 'AMPLITUDE' not in line) or 'TOTAL TIME WINDOW' in line:
                window = float(line.rstrip('\n\r ').split(' ')[-1])
            if 'TIMEZERO' in line or 'TIMEZERO AT POINT' in line:
                pe_data.trig = int(float(line.rstrip('\n\r ').split(' ')[-1])
                                   ) * np.ones((pe_data.tnum,))
            if i == 4 and pe_data.version == '1.0':
                try:
                    doy = (int(line[6:10]), int(line[1:2]), int(line[3:5]))
                except ValueError:
                    doy = (int(line[:4]), int(line[5:7]), int(line[8:10]))
            if i == 2 and pe_data.version == '1.5.340':
                doy = (int(line[6:10]), int(line[:2]), int(line[3:5]))

    if pe_data.version == '1.0':
        pe_data.data = np.zeros((pe_data.snum, pe_data.tnum), dtype=np.int16)
    elif pe_data.version == '1.5.340':
        pe_data.data = np.zeros((pe_data.snum, pe_data.tnum), dtype=np.float32)

    pe_data.traceheaders = TraceHeaders(pe_data.tnum)
    with open(true_fn, 'rb') as fin:
        lines = fin.read()

    offset = 0
    for i in range(pe_data.tnum):
        pe_data.traceheaders.get_header(offset, lines)
        offset += 25 * 4 + 28
        if pe_data.version == '1.0':
            trace = struct.unpack('<{:d}h'.format(pe_data.snum),
                                  lines[offset: offset + pe_data.snum * 2])
            offset += pe_data.snum * 2
        elif pe_data.version == '1.5.340':
            fmt = '<%df' % (len(lines[offset: offset + pe_data.snum * 4]) // 4)
            trace = struct.unpack(fmt, lines[offset:offset + pe_data.snum * 4])
            offset += pe_data.snum * 4

        trace -= np.nanmean(trace[:100])
        pe_data.data[:, i] = trace.copy()

    # known vars that are not really set
    pe_data.chan = 1
    pe_data.trace_num = np.arange(pe_data.tnum) + 1
    pe_data.trig_level = 0.
    pe_data.pressure = np.zeros((pe_data.tnum, ))
    pe_data.flags = RadarFlags()

    # Power some more real variables
    pe_data.dt = window / pe_data.snum * 1.0e-9
    pe_data.travel_time = np.atleast_2d(
        np.arange(0, window / 1.e3, pe_data.dt * 1.0e6)).transpose()
    pe_data.travel_time += pe_data.dt * 1.0e6

    # Now deal with the gps info
    if os.path.exists(gps_fn):
        pe_data.gps_data = _get_gps_data(gps_fn, pe_data.trace_num)
        pe_data.lat = pe_data.gps_data.lat
        pe_data.long = pe_data.gps_data.lon
        pe_data.x_coord = pe_data.gps_data.x
        pe_data.y_coord = pe_data.gps_data.y
        pe_data.dist = pe_data.gps_data.dist.flatten()
        pe_data.elev = pe_data.gps_data.z
        day_offset = datetime.datetime(doy[0], doy[1], doy[2], 0, 0, 0)
        tmin = day_offset.toordinal() + np.min(pe_data.gps_data.dectime) + 366.
        tmax = day_offset.toordinal() + np.max(pe_data.gps_data.dectime) + 366.
        # 366 for matlab compat
        pe_data.decday = np.linspace(tmin, tmax, pe_data.tnum)
        pe_data.trace_int = np.hstack((np.array(np.nanmean(
            np.diff(pe_data.dist))), np.diff(pe_data.dist)))
    else:
        print('Warning: Cannot find gps file, %s.' % gps_fn)
        pe_data.lat = np.zeros((pe_data.data.shape[1],))
        pe_data.long = np.zeros((pe_data.data.shape[1],))
        pe_data.x_coord = np.zeros((pe_data.data.shape[1],))
        pe_data.y_coord = np.zeros((pe_data.data.shape[1],))
        pe_data.dist = np.zeros((pe_data.data.shape[1],))
        pe_data.elev = np.zeros((pe_data.data.shape[1],))
        pe_data.decday = np.arange(pe_data.data.shape[1])
        pe_data.trace_int = np.ones((pe_data.data.shape[1],))

    pe_data.check_attrs()
    return pe_data
