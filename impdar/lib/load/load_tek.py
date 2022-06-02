#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

Load TEK files (e.g. 01-02 UW Kamb data) and convert to the .mat ImpDAR format

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

December 9 2020

"""

import numpy as np
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags


def fread(fid, nelements, dtype):
    """Equivalent to Matlab fread function"""
    if dtype is np.str:
        dt = np.uint8  # WARNING: assuming 8-bit ASCII for np.str!
    else:
        dt = dtype
    data_array = np.fromfile(fid, dt, nelements)
    data_array.shape = (nelements, 1)
    return data_array


def load_tek(fn_tek, magnets_per_wheel=1, wheel_diameter=0.5, trigger_level=0.1, trigger_sample=None,
            channel=1, *args, **kwargs):
    """Load a TEK file into ImpDAR
    Read file following http://dx.doi.org/10.7265/N5736NTS"""

    # Import the matlab file and setup a RadarData class to fill
    tek_data = RadarData(None)
    tek_data.fn = fn_tek

    # Initialize empty arrays
    tek_data.decday = np.array([]).astype(np.float32)
    wheel_count = np.array([]).astype(np.ushort)
    tek_data.pressure = np.array([]).astype(np.short)
    yinc = np.array([]).astype(np.float32)
    xinc = np.array([]).astype(np.float32)
    averages = np.array([]).astype(np.ushort)
    length = np.array([]).astype(np.ushort)
    tek_data.data = np.empty((0,1000)).astype(np.ushort)

    # Open file to read in the data
    fid = open(fn_tek,'rb')
    more_lines=True
    while more_lines:
        try:
            tek_data.decday = np.append(tek_data.decday,fread(fid, 1, np.float32))
            wheel_count = np.append(wheel_count,fread(fid, 1, np.ushort))
            tek_data.pressure = np.append(tek_data.pressure,fread(fid, 1, np.short))
            yinc = np.append(yinc,fread(fid, 1, np.float32))
            xinc = np.append(xinc,fread(fid, 1, np.float32))
            averages = np.append(averages,fread(fid, 1, np.ushort))
            length = np.append(length,fread(fid, 1, np.ushort))
            tek_data.data = np.append(tek_data.data,
                            np.transpose(fread(fid, length[-1], np.ushort)),axis=0)
        except:
            more_lines=False
    fid.close()

    # Transpose and normalize data around zero
    tek_data.data = np.transpose(tek_data.data)
    tek_data.data.dtype = np.short
    tek_data.data -= 512
    tek_data.snum = tek_data.data.shape[0]
    tek_data.tnum = tek_data.data.shape[1]
    tek_data.trace_num = np.arange(tek_data.tnum)

    # Distance array from wheel count
    tek_data.dist = wheel_count.astype(np.float32)
    tek_data.dist *= np.pi*wheel_diameter/magnets_per_wheel
    tek_data.trace_int = np.gradient(tek_data.dist)

    # Time step
    tek_data.dt = np.median(xinc)

    # Trigger sample
    tek_data.trig_level = trigger_level
    if trigger_sample is None:
        avg_trace = np.mean(tek_data.data,axis=1)
        exceeds_level = np.abs(np.gradient(avg_trace))>tek_data.trig_level*np.max(np.abs(avg_trace))
        trigger_sample = next(x[0] for x in enumerate(exceeds_level) if x[1] > 0.7)
    tek_data.trig = trigger_sample*np.ones(tek_data.tnum)
    tek_data.travel_time = (-trigger_sample+np.arange(tek_data.snum))*tek_data.dt*1e6

    # convert Pressure to an array of differences from 1
    tek_data.pressure -= tek_data.pressure[0]

    # Get variables that are already in the St Olaf .mat file
    tek_data.chan = channel

    # Initialize the flags
    tek_data.flags = RadarFlags()

    # Check attributed in the RadarData class and return
    tek_data.check_attrs()
    return tek_data
