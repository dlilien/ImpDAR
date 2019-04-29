#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import sys
import os
import unittest
import numpy as np
from impdar.lib.RadarData import RadarData
from impdar.lib import gpslib
if sys.version_info[0] >= 3:
    from unittest.mock import patch, MagicMock
else:
    from mock import patch, MagicMock
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class NoInitRadarData(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        self.data = np.zeros((10, 20))
        # need to set this to avoid divide by zero later
        self.dt = 1
        self.dist = np.arange(20)
        self.tnum = 20
        self.trace_num = np.arange(self.tnum) + 1.
        self.snum = 10
        self.travel_time = np.arange(10)
        self.decday = np.arange(20)


class TestGPS(unittest.TestCase):

    @patch('impdar.lib.gpslib.kinematic_gps_control')    
    def test_kinematic_gps_csv(self, mock_kgc):
        dats = [NoInitRadarData()]
        gpslib.kinematic_gps_csv(dats, os.path.join(THIS_DIR, 'input_data', 'gps_control.csv'))
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), mock_kgc.call_args[0][1]))
        self.assertTrue(np.allclose(np.arange(0, 200, 10), mock_kgc.call_args[0][2]))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), mock_kgc.call_args[0][3]))
        self.assertTrue(np.allclose(np.arange(0, 20, 1), mock_kgc.call_args[0][4]))
        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_csv(dats, os.path.join(THIS_DIR, 'input_data', 'gps_control.csv'), names='dumbanddummer')
        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_csv(dats, os.path.join(THIS_DIR, 'input_data', 'gps_control.csv'), names='dumb,dumb,and,dummer')

    def test_kinematic_gps_control(self):
        dats = [NoInitRadarData()]
        gpslib.kinematic_gps_control(dats, np.arange(0, 2.0, 0.1), np.arange(0, 200, 10), np.arange(0, 2000, 100), np.arange(0, 20, 1))
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), dats[0].lat))
        self.assertTrue(np.allclose(np.arange(0, 200, 10), dats[0].long))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), dats[0].elev))

        dat = NoInitRadarData()
        gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(0, 200, 10), np.arange(0, 2000, 100), np.arange(0, 20, 1))
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), dat.lat))
        self.assertTrue(np.allclose(np.arange(0, 200, 10), dat.long))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), dat.elev))

        dat = NoInitRadarData()
        dat.decday = dat.decday + 10
        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(0, 200, 10), np.arange(0, 2000, 100), np.arange(0, 20, 1))

    @patch('impdar.lib.gpslib.kinematic_gps_control')    
    def test_kinematic_gps_mat(self, mock_kgc):
        dats = [NoInitRadarData()]
        gpslib.kinematic_gps_mat(dats, os.path.join(THIS_DIR, 'input_data', 'gps_control.mat'))
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), mock_kgc.call_args[0][1]))
        self.assertTrue(np.allclose(np.arange(0, 200, 10), mock_kgc.call_args[0][2]))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), mock_kgc.call_args[0][3]))
        self.assertTrue(np.allclose(np.arange(0, 20, 1), mock_kgc.call_args[0][4]))

        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_mat(dats, os.path.join(THIS_DIR, 'input_data', 'gps_control_badfields.mat'))

    @patch('impdar.lib.gpslib.kinematic_gps_mat')    
    @patch('impdar.lib.gpslib.kinematic_gps_csv')    
    def test_interp(self, mock_kgc, mock_kgm):
        dats = [NoInitRadarData()]
        dats[0].constant_space = MagicMock()
        gpslib.interp(dats, 10., fn='dum.csv')
        self.assertTrue(len(mock_kgc.mock_calls) > 0)
        self.assertTrue(len(dats[0].constant_space.mock_calls) > 0)
        dats[0].constant_space = MagicMock()
        gpslib.interp(dats, 10., fn='dum.mat')
        self.assertTrue(len(mock_kgm.mock_calls) > 0)
        self.assertTrue(len(dats[0].constant_space.mock_calls) > 0)
        with self.assertRaises(ValueError):
            gpslib.interp(dats, 10., fn='dum.badext')

        dats[0].constant_space = MagicMock()
        gpslib.interp(dats, 10.)
        self.assertTrue(len(dats[0].constant_space.mock_calls) > 0)


if __name__ == '__main__':
    unittest.main()
