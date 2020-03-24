#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

Load files from the MALA RAMAC system.
These are saved as:
1) a header file (.rad)
2) a gps file (.cor)
3) a data file (.rd3)

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Nov 1 2019

"""

import numpy as np

import os.path
import struct
import datetime
from scipy.interpolate import interp1d

from ..RadarData import RadarData
from ..RadarFlags import RadarFlags


def load_ramac(fn_ram):
    """Load a MALA RAMAC file into ImpDAR"""

    ram_data = RadarData(None)

    # Read from header file
    if os.path.exists(os.path.splitext(fn_ram)[0]+'.rad'):
        with open(os.path.splitext(fn_ram)[0]+'.rad','r') as fid:
            hd_lines = fid.readlines()
    else:
        raise TypeError('No header file:',os.path.splitext(fn_ram)[0]+'.rad')

    for line in hd_lines:
        # number of samples
        if 'SAMPLES:' in line:
            ram_data.snum = int(line.split(':')[1][:-1])
        # sampling freqency, i.e. time step
        elif 'FREQUENCY:' in line:
            ram_data.dt = 1./float(line.split(':')[1][:-1])*1e-6
        # number of traces
        elif 'LAST TRACE:' in line:
            ram_data.tnum = int(line.split(':')[1][:-1])

    # define the travel time from samples and time step
    ram_data.travel_time = ram_data.dt*1e6*np.arange(ram_data.snum)

    # Read from GPS file
    if os.path.exists(os.path.splitext(fn_ram)[0]+'.cor'):
        with open(os.path.splitext(fn_ram)[0]+'.cor','r') as fid:
            gps_lines = fid.readlines()
    else:
        raise TypeError('No GPS file:',os.path.splitext(fn_ram)[0]+'.rcor')

    """ The gps file is written only for every ~third trace.
    append them all to an array and then do a linear interpolation."""
    lat = []
    lon = []
    elev = []
    decday = []
    trace_num = []
    for i,line in enumerate(gps_lines):
        line = gps_lines[i]
        line_arr = line.split('\t')
        trace_num.append(int(line_arr[0]))
        year = int(line_arr[1][:4])
        month = int(line_arr[1][5:7])
        day = int(line_arr[1][8:10])
        hour = int(line_arr[2][:2])
        minute = int(line_arr[2][3:5])
        second = int(line_arr[2][6:8])
        doy = datetime.datetime(year,month,day).toordinal()
        decday.append(doy + (hour+(minute+second/60.)/60.)/24.)

        hemisphere = line_arr[4]
        if hemisphere == 'N':
            lat.append(float(line_arr[3]))
        elif hemisphere == 'S':
            lat.append(-float(line_arr[3]))
        lon.append(float(line_arr[5]))
        elev.append(float(line_arr[7]))

    # interpolate to get gps information for all traces
    ram_data.trace_num = np.arange(ram_data.tnum) + 1
    decday_interp = interp1d(trace_num,decday,fill_value='extrapolate')
    ram_data.decday = decday_interp(ram_data.trace_num)
    lat_interp = interp1d(trace_num,lat,fill_value='extrapolate')
    ram_data.lat = lat_interp(ram_data.trace_num)
    long_interp = interp1d(trace_num,lon,fill_value='extrapolate')
    ram_data.long = long_interp(ram_data.trace_num)
    elev_interp = interp1d(trace_num,elev,fill_value='extrapolate')
    ram_data.elev = elev_interp(ram_data.trace_num)

    # Read from data file
    with open(fn_ram,'rb') as fid:
        lines = fid.read()
    ram_data.data = np.empty((ram_data.snum,ram_data.tnum))
    # read trace-by-trace from the binary and put them in their appropriate location in the array
    for tr in range(ram_data.tnum):
        trace = struct.unpack('<{:d}h'.format(ram_data.snum),lines[tr*2*ram_data.snum:(tr+1)*2*ram_data.snum])
        ram_data.data[:,tr] = trace

    # fill in the last values that are less important here
    ram_data.chan = 1
    ram_data.pressure = np.zeros_like(ram_data.trace_num)
    ram_data.trace_int = np.zeros_like(ram_data.trace_num)
    ram_data.trig = np.zeros_like(ram_data.trace_num)
    ram_data.trig_level = np.zeros_like(ram_data.trace_num)

    ram_data.flags = RadarFlags()
    ram_data.check_attrs()

    return ram_data
