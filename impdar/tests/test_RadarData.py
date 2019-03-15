#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the basics of RadarData
"""
import os
import unittest
import numpy as np
from impdar.lib.RadarData import RadarData

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class NoInitRadarData(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        self.data = np.array([[2, 2], [1, 1]])
        # need to set this to avoid divide by zero later
        self.dt = 1


class TestRadarDataLoading(unittest.TestCase):

    def test_ReadSucceeds(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.assertEqual(data.data.shape, (20, 40))

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_out.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))


class TestRadarDataMethods(unittest.TestCase):

    def setUp(self):
        self.data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))

    def test_Reverse(self):
        self.data.x_coord = np.array([1, 2])
        data_unrev = self.data.data.copy()
        self.data.reverse()
        self.assertTrue(np.allclose(self.data.data, np.fliplr(data_unrev)))
        self.assertTrue(np.allclose(self.data.x_coord, np.array([2, 1])))
        self.data.reverse()
        self.assertTrue(np.allclose(self.data.data, data_unrev))
        self.assertTrue(np.allclose(self.data.x_coord, np.array([1, 2])))

    def test_NMO_noexcpetion(self):
        # If velocity is 2
        self.data.nmo(0., uice=2.0, uair=2.0)
        self.assertTrue(np.allclose(self.data.travel_time * 1.0e-6, self.data.nmo_depth))
        # shouldn't care about uair if offset=0
        self.data.nmo(0., uice=2.0, uair=200.0)
        self.assertTrue(np.allclose(self.data.travel_time * 1.0e-6, self.data.nmo_depth))

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_out.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))


if __name__ == '__main__':
    unittest.main()
