#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""

import os
import unittest
import numpy as np
from impdar.lib.RadarData import RadarData
from impdar.lib._RadarDataSaving import conversions_enabled
from impdar.lib.RadarFlags import RadarFlags
from impdar.lib.Picks import Picks

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class NoInitRadarData(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        self.data = np.array([[2, 2], [1, 1]])
        # need to set this to avoid divide by zero later
        self.elevation = np.zeros((2,))
        self.lat = np.zeros((2,))
        self.long = np.zeros((2,))
        self.dt = 0.1
        self.tnum = self.data.shape[1]
        self.snum = self.data.shape[0]
        self.travel_time = 0.001 * np.arange(self.data.shape[0]) + 0.001
        self.dt = 1


class TestRadarData(unittest.TestCase):

    def test_WriteNoFLags(self):
        rd = NoInitRadarData()
        rd.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    def testWriteWithFlags(self):
        rd = NoInitRadarData()
        rd.flags = RadarFlags()
        rd.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    def test_WriteRead(self):
        # We are going to create a really bad file (most info missing) and see if we recover it or get an error
        rd = NoInitRadarData()
        rd.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    @unittest.skipIf(not conversions_enabled, 'No GDAL on this version')
    def test_output_shp_nolayers(self):
        rd = NoInitRadarData()
        rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test_gssi.shp'))

    @unittest.skipIf(not conversions_enabled, 'No GDAL on this version')
    def test_output_shp_picks(self):
        rd = NoInitRadarData()
        rd.Picks = Picks(rd)
        rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test_gssi.shp'))

    def tearDown(self):
        for fn in ['test_out.mat', 'test_gssi.shp', 'test_gssi.shx', 'test_gssi.prj', 'test_gssi.dbf']:
            if os.path.exists(os.path.join(THIS_DIR, 'input_data', fn)):
                os.remove(os.path.join(THIS_DIR, 'input_data', fn))


if __name__ == '__main__':
    unittest.main()
