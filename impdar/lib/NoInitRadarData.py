#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""Define a class with dummy data for use with testing."""
import numpy as np
from .RadarData import RadarData
from .RadarFlags import RadarFlags

DATA_DUMMY = np.ones((500, 400))


class NoInitRadarData(RadarData):
    """Define a very simple radar data class for testing.

    Data size is 2x2 by default.

    Parameters
    ----------
    big: bool
        Use a 10x20 array so more operations are reasonable
    """

    # This only exists so we can do tests on writing without reading

    def __init__(self, big=False):
        super(NoInitRadarData, self).__init__(None)
        if not big:
            self.data = np.array([[2, 2], [1, 1]])
            self.travel_time = 0.001 * np.arange(self.data.shape[0]) + 0.001
        else:
            self.data = np.zeros((10, 20))
            self.travel_time = np.arange(self.data.shape[0])

        self.fn = ''
        self.tnum = self.data.shape[1]
        self.snum = self.data.shape[0]
        self.dist = np.arange(self.tnum,)
        self.elevation = np.zeros((self.tnum,))
        self.long = np.arange(self.tnum) * 3.
        self.lat = np.arange(self.tnum) * 2.
        self.trace_num = np.arange(self.tnum) + 1.
        self.decday = np.arange(self.tnum)
        self.trace_int = 1
        self.dt = 1
        self.trig = np.zeros((self.tnum,))
        self.pressure = np.zeros((self.tnum,))


class NoInitRadarDataFiltering(RadarData):
    """Larger dummy data for testing.

    Data size is 500x400
    """

    # This only exists so we can do tests on writing without reading

    def __init__(self):
        super(NoInitRadarDataFiltering, self).__init__(None)
        self.fn = ''
        self.data = DATA_DUMMY.copy()
        self.dt = 0.1
        self.tnum = self.data.shape[1]
        self.snum = self.data.shape[0]
        self.travel_time = 0.001 * np.arange(self.data.shape[0]) + 0.001
        self.trace_num = np.arange(self.tnum) + 1.
        self.dt = 0.001e-6
        self.trace_int = self.dt * np.ones((self.tnum, ))
        self.flags = RadarFlags()
        self.hfilt_target_output = DATA_DUMMY * np.atleast_2d(1. - np.exp(
            -self.travel_time.flatten() * 0.05) / np.exp(
                -self.travel_time[0] * 0.05)).transpose()
        pexp = np.exp(-self.travel_time.flatten() * 0.05) / np.exp(
            -self.travel_time[0] * 0.05)
        pexp = pexp - pexp[-1]
        pexp = pexp / np.max(pexp)
        self.pexp_target_output = DATA_DUMMY * np.atleast_2d(
            1. - pexp).transpose()
        self.ahfilt_target_output = np.zeros_like(DATA_DUMMY)
        self.long = np.arange(self.tnum) * 3.
        self.lat = np.arange(self.tnum) * 2.
        self.x_coord = np.arange(self.tnum) * 3.
        self.y_coord = np.arange(self.tnum) * 2.
        self.decday = np.arange(self.tnum)
        self.elev = np.arange(self.tnum) * 0.001 + 100
        self.trig = np.zeros_like(self.elev).astype(int)
