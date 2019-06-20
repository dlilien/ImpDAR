#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com
#
# Distributed under terms of the GNU GPL3.0 License.

"""
Load a NetCDF of MCoRDS data
"""
from . import RadarData
import numpy as np

try:
    from netCDF4 import Dataset
    nc = True
except ImportError:
    nc = False


class MCorDSNC(RadarData.RadarData):

    def __init__(self, fn):
        dst = Dataset(fn, 'r')
        self.data = dst.variables['amplitude'][:]
        self.long = dst.variables['lon'][:]
        self.lat = dst.variables['lat'][:]
        self.time = dst.variables['time'][:]
        self.travel_time = dst.variables['fasttime'][:]
        self.dt = np.mean(np.diff(self.travel_time)) * 1.0e-6
        size = dst.variables['amplitude'].matlab_size
        self.snum, self.tnum = int(size[0]), int(size[1])


def load_mcords_nc(fn):
    return MCorDSNC(fn)
