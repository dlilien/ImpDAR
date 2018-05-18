#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the MIT license.

"""
Define a class that just has the necessary attributes for a stodeep file--this should be subclassed per filetype
"""
from scipy.io import loadmat, savemat


class StODeep():
    # First define all the common attributes
    chan = None
    data = None
    decday = None
    dist = None
    dt = None
    elev = None
    flags = None
    lat = None
    long = None
    pressure = None
    snum = None
    tnum = None
    trace_int = None
    trace_num = None
    travel_time = None
    trig = None
    trig_level = None
    x_coord = None
    y_coord = None

    # Now make some load/save methods that will work with the matlab format
    def __init__(self, fn):
        mat = loadmat(fn)
        for attr in ['chan', 'data', 'decday', 'dist', 'dt', 'elev', 'flags', 'lat', 'long', 'pressure', 'snum', 'tnum', 'trace_int', 'trace_num', 'travel_time', 'trig', 'trig_level', 'x_coord', 'y_coord']:
            setattr(self, attr, mat[attr])

    def save(self, fn):
        mat = {}
        for attr in ['chan', 'data', 'decday', 'dist', 'dt', 'elev', 'flags', 'lat', 'long', 'pressure', 'snum', 'tnum', 'trace_int', 'trace_num', 'travel_time', 'trig', 'trig_level', 'x_coord', 'y_coord']:
            mat[attr] = getattr(self, attr)
        savemat(fn, mat)


class StOFlags():

    def __init__(self):
        self.batch = 0
        self.bpass = 0
        self.hfilt = 0
        self.rgain = 0
        self.agc = 0
        self.restack = 0
        self.reverse = 0
        self.crop = 0
        self.nmo = 0
        self.interp = 0
        self.mig = 0
        self.elev = 0
