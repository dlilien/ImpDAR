#! /usr/bin/env python
#function [chan,data,decday,dist,dt,elev,...
#        flags,lat,long,pressure,snum,tnum,trace_int,...
#        trace_num,travel_time,trig,trig_level,x_coord,y_coord] = read_pulse_ekko(filename,sflag)
#
#this is a not done function to read pulse ekko data into matlab and
#convert to stro_radar format
#Knut Christianson 6 April 2017; 20 May 2017
# Pythonized by David Lilien, May 2018
# Distributed under the GNU GPL3 license

import os.path
from .gpslib import RadarGPS
import numpy as np
from .RadarData import RadarData, RadarFlags
import struct
import datetime


class PE(RadarData):

    def __init__(self, fn):
        bn = os.path.splitext(fn)[0]
        hdname = bn + '.HD'
        true_fn = bn + '.DT1'
        gps_fn = bn + '.GPS'

        with open(hdname, 'r') as fin:
            for i, line in enumerate(fin):
                if 'TRACES' in line:
                    self.tnum = int(line.rstrip('\n\r ').split(' ')[-1])
                if 'PTS' in line:
                    self.snum = int(line.rstrip('\n\r ').split(' ')[-1])
                if 'WINDOW' in line:
                    window = float(line.rstrip('\n\r ').split(' ')[-1])
                if 'TIMEZERO' in line:
                    self.trig = int(float(line.rstrip('\n\r ').split(' ')[-1]))
                if i == 4:
                    doy = (int(line[:4]), int(line[5:7]), int(line[8:10]))

        self.data = np.zeros((self.snum, self.tnum), dtype=np.int16)
        self.traceheaders = TH(self.tnum)
        with open(true_fn, 'rb') as fin:
            lines = fin.read()

        offset = 0
        for i in range(self.tnum):
            self.traceheaders.get_header(offset, lines)
            offset += 25 * 4 + 28
            self.data[:, i] = struct.unpack('<{:d}h'.format(self.snum), lines[offset: offset + self.snum * 2])
            offset += self.snum * 2

        # known vars that are not really set
        self.chan = 1
        self.trace_num = np.arange(self.tnum) + 1
        self.trig_level = np.zeros((self.tnum, ))
        self.pressure = np.zeros((self.tnum, ))
        self.flags = RadarFlags()

        # some more real variables
        self.dt = window / self.snum * 1.0e-9
        self.travel_time = np.atleast_2d(np.arange(0, window / 1.0e3, self.dt * 1.0e6)).transpose() + self.dt * 1.0e6

        # Now deal with the gps info
        self.gps_data = _get_gps_data(gps_fn, self.trace_num)
        self.lat = self.gps_data.lat
        self.long = self.gps_data.lon
        self.x_coord = self.gps_data.x
        self.y_coord = self.gps_data.y
        self.dist = self.gps_data.dist.flatten()
        self.elev = self.gps_data.z

        timezero = datetime.datetime(2017, 1, 1, 0, 0, 0)
        day_offset = datetime.datetime(doy[0], doy[1], doy[2], 0, 0, 0) - timezero
        tmin, tmax = day_offset.days + np.min(self.gps_data.dectime), day_offset.days + np.max(self.gps_data.dectime)
        self.decday = np.linspace(tmin, tmax, self.tnum)
        self.trace_int = np.hstack((np.array(np.nanmean(np.diff(self.dist))), np.diff(self.dist)))


class TH:

    def __init__(self, tnum):
        self.header_index = 0
        self.TraceNumbers = np.zeros((1, tnum))
        self.Positions = np.zeros((1, tnum))
        self.NumberOfPointsPerTrace = np.zeros((1, tnum))
        self.Topography = np.zeros((1, tnum))
        self.NumberOfBytesPerPoint = np.zeros((1, tnum))
        self.NumberOfStacks = np.zeros((1, tnum))
        self.TimeWindow = np.zeros((1, tnum))
        self.PosX = np.zeros((1, tnum))
        self.PosY = np.zeros((1, tnum))
        self.PosZ = np.zeros((1, tnum))
        self.RX = np.zeros((1, tnum))
        self.RY = np.zeros((1, tnum))
        self.RZ = np.zeros((1, tnum))
        self.TX = np.zeros((1, tnum))
        self.TY = np.zeros((1, tnum))
        self.TZ = np.zeros((1, tnum))
        self.TimeZeroAdjustment = np.zeros((1, tnum))
        self.ZeroFlag = np.zeros((1, tnum))
        self.TimeOfDay = np.zeros((1, tnum))
        self.CommentFlag = np.zeros((1, tnum))
        self.Comment = ['' for i in range(tnum)]

    def get_header(self, offset, f):
        header = struct.unpack('<25f', f[offset: offset + 25 * 4])
        comment = struct.unpack('<28c', f[offset + 25 * 4: offset + 25 * 4 + 28])
        self.TraceNumbers[0, self.header_index] = header[0]
        self.Positions[0, self.header_index] = header[1]
        self.NumberOfPointsPerTrace[0, self.header_index] = header[2]
        self.Topography[0, self.header_index] = header[3]
        self.NumberOfBytesPerPoint[0, self.header_index] = header[5]
        self.NumberOfStacks[0, self.header_index] = header[7]
        self.TimeWindow[0, self.header_index] = header[8]
        self.PosX[0, self.header_index] = header[9]
        self.PosY[0, self.header_index] = header[11]
        self.PosZ[0, self.header_index] = header[13]
        self.RX[0, self.header_index] = header[14]
        self.RY[0, self.header_index] = header[15]
        self.RZ[0, self.header_index] = header[16]
        self.TX[0, self.header_index] = header[17]
        self.TY[0, self.header_index] = header[18]
        self.TZ[0, self.header_index] = header[19]
        self.TimeZeroAdjustment[0, self.header_index] = header[20]
        self.ZeroFlag[0, self.header_index] = header[21]
        self.TimeOfDay[0, self.header_index] = header[23]
        self.CommentFlag[0, self.header_index] = header[24]
        self.Comment[self.header_index] = str(comment[0])
        self.header_index += 1


def _get_gps_data(fn, trace_nums):
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

    with open(fn) as f:
        lines = f.readlines()
    ggis = lines[::2]
    gga = lines[1::2]
    scans = np.array(list(map(lambda x: int(float(x.rstrip('\n\r ').split(' ')[-1])), ggis)))
    data = RadarGPS(gga, scans, trace_nums)
    return data


def load_pe(fn, *args, **kwargs):
    return PE(fn)
