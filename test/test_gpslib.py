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
from impdar.lib.NoInitRadarData import NoInitRadarData
from impdar.lib import gpslib
if sys.version_info[0] >= 3:
    from unittest.mock import patch, MagicMock
else:
    from mock import patch, MagicMock
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestGPS(unittest.TestCase):

    @patch('impdar.lib.gpslib.kinematic_gps_control')
    def test_kinematic_gps_csv(self, mock_kgc):
        dats = [NoInitRadarData(big=True)]
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
        dats = [NoInitRadarData(big=True)]
        gpslib.kinematic_gps_control(dats, np.arange(0, 2.0, 0.1), np.arange(40, 60., 1.), np.arange(0., 2000., 100.), np.arange(0., 20., 1.), guess_offset=False)
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), dats[0].lat))
        self.assertTrue(np.allclose(np.arange(40, 60, 1), dats[0].long))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), dats[0].elev))

        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_control(dats, np.arange(0, 2.0, 0.1), np.arange(0, 200, 10), np.arange(0, 2000, 100), np.arange(0, 20, 1), guess_offset=True)

        dat = NoInitRadarData(big=True)
        gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(40, 60, 1), np.arange(0, 2000, 100), np.arange(0, 20, 1), guess_offset=False)
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), dat.lat))
        self.assertTrue(np.allclose(np.arange(40, 60, 1), dat.long))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), dat.elev))

        # We should be allowed to be off by 360 in longitude
        dat = NoInitRadarData(big=True)
        gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(40, 60, 1) - 360., np.arange(0, 2000, 100), np.arange(0, 20, 1), guess_offset=False)
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), dats[0].lat))
        self.assertTrue(np.allclose(np.arange(40, 60, 1), dats[0].long))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), dats[0].elev))

        # and off the other way
        dat = NoInitRadarData(big=True)
        dat.long = dat.long - 360.
        gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(40, 60, 1), np.arange(0, 2000, 100), np.arange(0, 20, 1), guess_offset=False)
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), dats[0].lat))
        self.assertTrue(np.allclose(np.arange(40, 60, 1), dats[0].long))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), dats[0].elev))

        # but not totally wrong
        dat = NoInitRadarData(big=True)
        dat.long = np.linspace(70, 71, len(dat.long))
        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(40, 60, 1), np.arange(0, 2000, 100), np.arange(0, 20, 1), guess_offset=True)

        dat = NoInitRadarData(big=True)
        gpslib.kinematic_gps_control(dat, np.arange(-1.0, 3.0, 0.1), np.arange(20, 60, 1), np.arange(-1000, 3000, 100), np.arange(-10, 30, 1), guess_offset=True)

        # Multiple inputs
        dats = [NoInitRadarData(big=True), NoInitRadarData(big=True)]
        gpslib.kinematic_gps_control(dats, np.arange(-1.0, 3.0, 0.1), np.arange(40, 80, 1), np.arange(-1000, 3000, 100), np.arange(-10, 30, 1), guess_offset=True)

        # Bad timing
        dat = NoInitRadarData(big=True)
        dat.decday += 10
        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(0, 20, 1), np.arange(0, 2000, 100), np.arange(0, 20, 1))

        dat = NoInitRadarData(big=True)
        dat.decday[10] = np.NaN

        gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(40, 60, 1), np.arange(0, 2000, 100), np.arange(0, 20, 1), guess_offset=False)
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1)[~np.isnan(dat.lat)], dat.lat[~np.isnan(dat.lat)]))
        self.assertTrue(np.allclose(np.arange(40, 60, 1)[~np.isnan(dat.long)], dat.long[~np.isnan(dat.long)]))

        gpslib.kinematic_gps_control(dat, np.arange(0, 2.0, 0.1), np.arange(40, 60, 1), np.arange(0, 2000, 100), np.arange(0, 20, 1), guess_offset=False, old_gps_gaps=True)
        self.assertTrue(np.isnan(dat.lat[10]))
        self.assertTrue(np.isnan(dat.long[10]))
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1)[~np.isnan(dat.lat)], dat.lat[~np.isnan(dat.lat)]))
        self.assertTrue(np.allclose(np.arange(40, 60, 1)[~np.isnan(dat.long)], dat.long[~np.isnan(dat.long)]))

        gpslib.kinematic_gps_control(dat, np.arange(-1.0, 3.0, 0.1), np.arange(20, 60, 1), np.arange(-1000, 3000, 100), np.arange(-10, 30, 1), guess_offset=True, old_gps_gaps=True)
        self.assertTrue(np.isnan(dat.lat[10]))
        self.assertTrue(np.isnan(dat.long[10]))

        # We should overcome a bad timestamp here
        dat.decday[11] += 2000.0
        cachelat = dat.lat[11]
        gpslib.kinematic_gps_control(dat, np.arange(-1.0, 3.0, 0.1), np.arange(20, 60, 1), np.arange(-1000, 3000, 100), np.arange(-10, 30, 1), guess_offset=True, old_gps_gaps=True)
        self.assertEqual(dat.lat[11], cachelat)

        # If we have bad times, we just should fail
        dat.decday += 2000.0
        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_control(dat, np.arange(-1.0, 3.0, 0.1), np.arange(20, 60, 1), np.arange(-1000, 3000, 100), np.arange(-10, 30, 1), guess_offset=True, old_gps_gaps=True)

    @patch('impdar.lib.gpslib.kinematic_gps_control')
    def test_kinematic_gps_mat(self, mock_kgc):
        dats = [NoInitRadarData(big=True)]
        gpslib.kinematic_gps_mat(dats, os.path.join(THIS_DIR, 'input_data', 'gps_control.mat'))
        self.assertTrue(np.allclose(np.arange(0, 2.0, 0.1), mock_kgc.call_args[0][1]))
        self.assertTrue(np.allclose(np.arange(0, 200, 10), mock_kgc.call_args[0][2]))
        self.assertTrue(np.allclose(np.arange(0, 2000, 100), mock_kgc.call_args[0][3]))
        self.assertTrue(np.allclose(np.arange(0, 20, 1), mock_kgc.call_args[0][4]))

        with self.assertRaises(ValueError):
            gpslib.kinematic_gps_mat(dats, os.path.join(THIS_DIR, 'input_data', 'gps_control_badfields.mat'), extrapolate=False)

    @patch('impdar.lib.gpslib.kinematic_gps_mat')
    @patch('impdar.lib.gpslib.kinematic_gps_csv')
    @patch('impdar.lib.gpslib.kinematic_gps_control')
    def test_interp(self, mock_kgctrl, mock_kgc, mock_kgm):
        dats = [NoInitRadarData(big=True)]
        dats[0].constant_space = MagicMock()
        gpslib.interp(dats, 10., fn='dum.csv')
        self.assertTrue(len(mock_kgc.mock_calls) > 0)
        self.assertTrue(len(dats[0].constant_space.mock_calls) > 0)
        dats[0].constant_space = MagicMock()
        gpslib.interp(dats, 10., fn='dum.mat')
        self.assertTrue(len(mock_kgm.mock_calls) > 0)
        self.assertTrue(len(dats[0].constant_space.mock_calls) > 0)
        with self.assertRaises(Exception):
            gpslib.interp(dats, 10., fn='dum.badext')

        dats = [NoInitRadarData(big=True)]
        dats[0].constant_space = MagicMock()
        gpslib.interp(dats, 10.)
        self.assertTrue(len(dats[0].constant_space.mock_calls) > 0)

        # We should be calculating distance if there is none
        dats = [NoInitRadarData(big=True)]
        dats[0].dist = None
        dats[0].constant_space = MagicMock()
        count = len(mock_kgctrl.mock_calls)
        gpslib.interp(dats, 10.)
        self.assertGreater(len(mock_kgctrl.mock_calls), count)

    @unittest.skipIf(not gpslib.conversions_enabled, 'No gdal')
    def test_conversions(self):
        pts = np.array([[-8., 10.], [-9., 11.], [-10., 12.]])
        conv_utm, _ = gpslib.get_utm_conversion(-8.0, 10.0)
        proj_pts = conv_utm(pts)
        self.assertTrue(np.all(~np.isnan(proj_pts)))

        pts = np.array([[-88., 10.], [-89., 11.], [-89.1, 12.]])
        conv_sps, _ = gpslib.get_conversion(t_srs='EPSG:3031')
        proj_pts = conv_sps(pts)
        self.assertTrue(np.all(~np.isnan(proj_pts)))

    @unittest.skipIf(gpslib.conversions_enabled, 'GDAL found, this is a failure test')
    def test_conversions_off(self):
        # we want to be able to import gpslib but later fail
        with self.assertRaises(ImportError):
            conv_utm = gpslib.get_utm_conversion(-8.0, 10.0)
        with self.assertRaises(ImportError):
            conv_sps = gpslib.get_conversion(t_srs='EPSG:3031')


if __name__ == '__main__':
    unittest.main()
