#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Define a class with dummy data for use with testing
"""
import numpy as np
from .RadarData import RadarData
from .RadarFlags import RadarFlags

DATA_DUMMY = np.ones((500, 400))


class NoInitRadarData(RadarData):
    """Define a very simple radar data class for testing"""
    # This only exists so we can do tests on writing without reading

    def __init__(self, big=False):
        super(NoInitRadarData, self).__init__(None)
        if not big:
            self.data = np.array([[2, 2], [1, 1]])
            self.travel_time = 0.001 * np.arange(self.data.shape[0]) + 0.001
        else:
            self.data = np.zeros((10, 20))
            self.travel_time = np.arange(self.data.shape[0])

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


class NoInitRadarDataFiltering(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        super(NoInitRadarDataFiltering, self).__init__(None)
        self.data = DATA_DUMMY.copy()
        self.dt = 0.1
        self.tnum = self.data.shape[1]
        self.snum = self.data.shape[0]
        self.travel_time = 0.001 * np.arange(self.data.shape[0]) + 0.001
        self.dt = 0.001e-6
        self.flags = RadarFlags()
        self.hfilt_target_output = DATA_DUMMY * np.atleast_2d(1. - np.exp(-self.travel_time.flatten() * 0.05) / np.exp(-self.travel_time[0] * 0.05)).transpose()
        pexp = np.exp(-self.travel_time.flatten() * 0.05) / np.exp(-self.travel_time[0] * 0.05)
        pexp = pexp - pexp[-1]
        pexp = pexp / np.max(pexp)
        self.pexp_target_output = DATA_DUMMY * np.atleast_2d(1. - pexp).transpose()
        self.ahfilt_target_output = np.zeros_like(DATA_DUMMY)
