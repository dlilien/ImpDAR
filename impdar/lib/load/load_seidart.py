#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""
Load SeiDarT .csv files and convert to the .mat ImpDAR file.

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 23 2019

"""

import numpy as np
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags


def load_seidart(fn_sd, fn_prj, seismic=False, *args, **kwargs):
    """Load a SeiDarT file into ImpDAR."""
    sd_data = RadarData(None)

    # open the SeiDarT file to get data
    sd_data.data = np.transpose(np.genfromtxt(fn_sd))

    # open the project file to read the time step
    with open(fn_prj, 'r') as fid:
        prj_contents = fid.read()
    if seismic:
        dt_start = prj_contents.find("S,dt,") + 5
    else:
        dt_start = prj_contents.find("E,dt,") + 5
    dt_end = prj_contents[dt_start:].find("\n") + dt_start
    sd_data.dt = prj_contents[dt_start:dt_end]

    # Remove pretrigger
    trig_threshold = 0.5  # trigger when mean trace gets up to 50% of maximum
    mean_trace = np.nanmean(np.abs(sd_data.data), axis=1)
    idx_threshold = np.argwhere(mean_trace > trig_threshold * np.nanmax(
        mean_trace))
    idx_trig = np.nanmin(idx_threshold)
    sd_data.data = sd_data.data[idx_trig:]

    # other variables are from the array shape
    sd_data.snum = sd_data.data.shape[0]
    sd_data.tnum = sd_data.data.shape[1]
    sd_data.trace_num = np.arange(sd_data.data.shape[1]) + 1
    sd_data.trig_level = np.zeros((sd_data.tnum,))
    sd_data.pressure = np.zeros((sd_data.tnum,))
    sd_data.flags = RadarFlags()
    sd_data.travel_time = sd_data.dt * 1.0e6 * np.arange(sd_data.snum)
    sd_data.trig = np.zeros((sd_data.tnum,))
    sd_data.lat = np.zeros((sd_data.tnum,))
    sd_data.long = np.zeros((sd_data.tnum,))
    sd_data.x_coord = np.zeros((sd_data.tnum,))
    sd_data.y_coord = np.zeros((sd_data.tnum,))
    sd_data.elev = np.zeros((sd_data.tnum,))
    sd_data.decday = np.arange(sd_data.tnum)
    sd_data.trace_int = np.ones((sd_data.tnum,))

    sd_data.dist = np.arange(sd_data.tnum)
    sd_data.chan = -99.
    sd_data.fn = fn_sd
    sd_data.check_attrs()
    return sd_data
