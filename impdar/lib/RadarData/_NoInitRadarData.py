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
from . import RadarData


class NoInitRadarData(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        super(NoInitRadarData, self).__init__(None)
        self.data = np.array([[2, 2], [1, 1]])
        # need to set this to avoid divide by zero later
        self.elevation = np.zeros((2,))
        self.lat = np.ones((2,)) * 89.
        self.long = np.ones((2,))
        self.dt = 0.1
        self.tnum = self.data.shape[1]
        self.snum = self.data.shape[0]
        self.travel_time = 0.001 * np.arange(self.data.shape[0]) + 0.001
        self.dt = 1
