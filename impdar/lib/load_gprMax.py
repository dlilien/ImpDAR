#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Load h5 files and convert to the .mat ImpDAR file

Author:
Benjamin Hills
benjaminhhills@gmail.com
University of Washington
Earth and Space Sciences

Mar 17 2019

"""

import h5py # TODO: make h5py a dependency?

import numpy as np
from .RadarData import RadarData, RadarFlags


class h5(RadarData):
    def __init__(self, fn):
        # open the h5 file
        with h5py.File(fn) as f:
            self.dt = f.attrs['dt']
            self.data = np.array(f['/rxs/rx1/Ez'])
        # other variables are from the array shape
        self.snum = self.data.shape[0]
        self.tnum = self.data.shape[1]
        self.trace_num = np.arange(self.data.shape[1]) + 1
        self.trig_level = np.zeros((self.tnum, ))
        self.pressure = np.zeros((self.tnum, ))
        self.flags = RadarFlags()
        self.travel_time = self.dt*np.arange(self.snum)
        self.trig = np.zeros((self.tnum,))
        self.lat = np.zeros((self.tnum,))
        self.long = np.zeros((self.tnum,))
        self.x_coord = np.zeros((self.tnum,))
        self.y_coord = np.zeros((self.tnum,))
        self.dist = np.zeros((self.tnum,))
        self.elev = np.zeros((self.tnum,))
        self.decday = np.arange(self.tnum)
        self.trace_int = np.ones((self.tnum,))

        for attr in ['chan', 'data', 'decday', 'dist', 'dt', 'elev', 'flags', 'lat', 'long', 'pressure', 'snum', 'tnum', 'trace_int', 'trace_num', 'travel_time', 'trig', 'trig_level', 'x_coord', 'y_coord']:
            if getattr(self, attr) is None:
                print(attr + ' is not defined')
                setattr(self, attr, 0)

def load_gprMax(fn, *args, **kwargs):
    return h5(fn)
