#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Read the data from an OSU deep radar file
"""

import datetime
import numpy as np

from ..RadarData import RadarData
from ..RadarFlags import RadarFlags
from .loading_utils import common_start

def load_osu(fns_osu,*args,**kwargs):
    """ load data from OSU deep radar files, written as .txt """
    osu_data = RadarData(None)
    # We want to be able to use this step concatenate a series of files numbered by the controller
    if isinstance(fns_osu, str):
        fns_osu = [fns_osu]
    osu_data.fn = common_start(fns_osu)

    dt_s = []
    osu_data.trace_num = []
    osu_data.lat = []
    osu_data.long = []
    osu_data.decday = []
    osu_data.elev = []
    stacks = []
    for i, fn_i in enumerate(fns_osu):
        osu_data.trace_num.append(fn_i)
        # Read the file into a string
        with open(fn_i, 'r') as fid:
            lines = fid.readlines()

        # Header information
        dt_s.append(float(lines[5]))

        # GPS information
        osu_data.lat = np.append(osu_data.lat,float(lines[2]))
        osu_data.long = np.append(osu_data.long,float(lines[3]))
        osu_data.elev = np.append(osu_data.elev,float(lines[4]))
        month = int(lines[0].split('/')[0])
        day = int(lines[0].split('/')[1])
        year = int(lines[0].split('/')[2][:-1])
        hour = int(lines[1].split(':')[0])
        minute = int(lines[1].split(':')[1])
        second = int(lines[1].split(':')[2][:-1])
        doy = datetime.datetime(year,month,day).toordinal() + 366. # 366 for matlab compat
        decday = doy + (hour+((minute+(second/60.))/60.))/24.
        osu_data.decday = np.append(osu_data.decday,decday)

        # Trace Data
        stacks.append(np.array(lines[6].split('\t')).astype(float))

    # data
    osu_data.data = np.transpose(stacks)
    osu_data.snum = osu_data.data.shape[0]
    osu_data.tnum = osu_data.data.shape[1]
    osu_data.trace_num = np.arange(osu_data.tnum) + 1

    # I don't know if we actually want to do this, but the filenaming scheme is wacky and this
    # will make any logical collection look good
    sort_idx = np.argsort(osu_data.decday)
    osu_data.data = osu_data.data[:,sort_idx]
    osu_data.lat = osu_data.lat[sort_idx]
    osu_data.long = osu_data.long[sort_idx]
    osu_data.elev = osu_data.elev[sort_idx]
    osu_data.decday = osu_data.decday[sort_idx]

    # time and time step
    if all(dt == dt_s[0] for dt in dt_s):
        osu_data.dt = dt_s[0]
        osu_data.travel_time = osu_data.dt*1e6*np.arange(osu_data.snum)
    else:
        raise ValueError('Trace headers have different time steps.')

    # known vars that are not really set
    osu_data.chan = 1
    osu_data.trace_int = np.zeros_like(osu_data.trace_num)
    osu_data.pressure = np.zeros_like(osu_data.trace_num)
    osu_data.trig_level = np.zeros_like(osu_data.trace_num)
    osu_data.trig = np.zeros_like(osu_data.trace_num)
    osu_data.flags = RadarFlags()

    osu_data.check_attrs()
    return osu_data
