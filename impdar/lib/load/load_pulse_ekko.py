#! /usr/bin/env python
# this is a not done function to read pulse ekko data into matlab and
# convert to stro_radar format
# Knut Christianson 6 April 2017; 20 May 2017
# Pythonized by David Lilien, May 2018
# Distributed under the GNU GPL3 license
"""Read Pulse Ekko data"""

import os.path
import struct
import datetime
import numpy as np

from ..gpslib import RadarGPS
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags


class TraceHeaders:
    """Class used internally to handle pulse-ekko headers"""

    def __init__(self, tnum):
        """Create a container for all the trace headers"""
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
        """Get the header information for a single trace"""
        header = struct.unpack('<25f', f_lines[offset: offset + 25 * 4])
        comment = struct.unpack('<28c', f_lines[offset + 25 * 4: offset + 25 * 4 + 28])
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
    """Read GPS data associated with a GSSI sir4000 file.

    Parameters
    ----------
    fn: str
        A dzg file with ggis and gga strings.
    rev: bool, optional
        Reverse the points in this file (used for concatenating radar files). Default False.

    Returns
    -------
    data: :class:`~impdar.lib.gpslib.nmea_info`
    """

    with open(fn_gps) as f_in:
        lines = f_in.readlines()
    ggis = lines[::2]
    gga = lines[1::2]
    scans = np.array(list(map(lambda x: int(float(x.rstrip('\n\r ').split(' ')[-1])), ggis)))
    data = RadarGPS(gga, scans, trace_nums)
    return data


def load_pe(fn_dt1, *args, **kwargs):
    """Load data from a pulse_ekko file"""
    pe_data = RadarData(None)
    bn_pe = os.path.splitext(fn_dt1)[0]
    hdname = bn_pe + '.HD'
    true_fn = bn_pe + '.DT1'
    gps_fn = bn_pe + '.GPS'

    with open(hdname, 'rU') as fin:
        for i, line in enumerate(fin):
            if 'TRACES' in line:
                pe_data.tnum = int(line.rstrip('\n\r ').split(' ')[-1])
            if 'PTS' in line:
                pe_data.snum = int(line.rstrip('\n\r ').split(' ')[-1])
            if 'WINDOW' in line:
                window = float(line.rstrip('\n\r ').split(' ')[-1])
            if 'TIMEZERO' in line:
                pe_data.trig = int(float(line.rstrip('\n\r ').split(' ')[-1]))
            if i == 4:
                doy = (int(line[:4]), int(line[5:7]), int(line[8:10]))

    pe_data.data = np.zeros((pe_data.snum, pe_data.tnum), dtype=np.int16)
    pe_data.traceheaders = TraceHeaders(pe_data.tnum)
    with open(true_fn, 'rb') as fin:
        lines = fin.read()

    offset = 0
    for i in range(pe_data.tnum):
        pe_data.traceheaders.get_header(offset, lines)
        offset += 25 * 4 + 28
        pe_data.data[:, i] = struct.unpack('<{:d}h'.format(pe_data.snum),
                                           lines[offset: offset + pe_data.snum * 2])
        offset += pe_data.snum * 2

    # known vars that are not really set
    pe_data.chan = 1
    pe_data.trace_num = np.arange(pe_data.tnum) + 1
    pe_data.trig_level = np.zeros((pe_data.tnum, ))
    pe_data.pressure = np.zeros((pe_data.tnum, ))
    pe_data.flags = RadarFlags()

    # some more real variables
    pe_data.dt = window / pe_data.snum * 1.0e-9
    pe_data.travel_time = np.atleast_2d(np.arange(0, window / 1.e3, pe_data.dt * 1.0e6)).transpose()
    pe_data.travel_time += pe_data.dt * 1.0e6

    # Now deal with the gps info
    pe_data.gps_data = _get_gps_data(gps_fn, pe_data.trace_num)
    pe_data.lat = pe_data.gps_data.lat
    pe_data.long = pe_data.gps_data.lon
    pe_data.x_coord = pe_data.gps_data.x
    pe_data.y_coord = pe_data.gps_data.y
    pe_data.dist = pe_data.gps_data.dist.flatten()
    pe_data.elev = pe_data.gps_data.z

    day_offset = datetime.datetime(doy[0], doy[1], doy[2], 0, 0, 0)
    tmin = day_offset.toordinal() + np.min(pe_data.gps_data.dectime) + 366.
    tmax = day_offset.toordinal() + np.max(pe_data.gps_data.dectime) + 366.  # 366 for matlab compat
    pe_data.decday = np.linspace(tmin, tmax, pe_data.tnum)
    pe_data.trace_int = np.hstack((np.array(np.nanmean(np.diff(pe_data.dist))),
                                   np.diff(pe_data.dist)))
    pe_data.check_attrs()
    return pe_data
