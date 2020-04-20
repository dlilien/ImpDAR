#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
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


class TestRadarDataLoading(unittest.TestCase):

    def test_ReadSucceeds(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.assertEqual(data.data.shape, (20, 40))

    def test_ReadLegacyStodeep(self):
        # This one has data and other attrs
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_otherstodeepattrs.mat'))
        self.assertEqual(data.data.shape, (20, 40))

        # This one has has only other attrs
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_just_otherstodeepattrs.mat'))
        self.assertEqual(data.data.shape, (20, 40))

    def test_badread(self):
        # Data but not other attrs
        with self.assertRaises(KeyError):
            data = RadarData(os.path.join(THIS_DIR, 'input_data', 'nonimpdar_matlab.mat'))

        # All other attrs, no data
        with self.assertRaises(KeyError):
            data = RadarData(os.path.join(THIS_DIR, 'input_data', 'nonimpdar_justmissingdat.mat'))

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_out.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))


class TestRadarDataMethods(unittest.TestCase):

    def setUp(self):
        self.data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.data.x_coord = np.arange(40)
        self.data.nmo_depth = None
        self.data.travel_time = np.arange(0, 0.2, 0.01)
        self.data.dt = 1.0e-8
        self.data.trig = self.data.trig * 0.

    def test_Reverse(self):
        data_unrev = self.data.data.copy()
        self.data.reverse()
        self.assertTrue(np.allclose(self.data.data, np.fliplr(data_unrev)))
        self.assertTrue(np.allclose(self.data.x_coord, np.arange(39, -1, -1)))
        self.data.reverse()
        self.assertTrue(np.allclose(self.data.data, data_unrev))
        self.assertTrue(np.allclose(self.data.x_coord, np.arange(40)))

    def test_CropTWTT(self):
        self.data.crop(0.165, 'bottom', dimension='twtt')
        self.assertTrue(self.data.data.shape == (17, 40))
        self.data.crop(0.055, 'top', dimension='twtt')
        self.assertTrue(self.data.data.shape == (11, 40))

        # do not fail on bad flags
        self.data.flags.crop = False
        self.data.crop(0.055, 'top', dimension='twtt')
        self.assertTrue(self.data.flags.crop.shape == (3,))
        self.assertTrue(self.data.flags.crop[0])

    def test_CropErrors(self):
        with self.assertRaises(ValueError):
            self.data.crop(0.165, 'bottom', dimension='dummy')
        with self.assertRaises(ValueError):
            self.data.crop(0.165, 'dummy', dimension='twtt')

    def test_CropSNUM(self):
        self.data.crop(17, 'bottom', dimension='snum')
        self.assertTrue(self.data.data.shape == (17, 40))
        self.data.crop(6, 'top', dimension='snum')
        self.assertTrue(self.data.data.shape == (11, 40))

    def test_CropTrigInt(self):
        self.data.trig = 2
        with self.assertRaises(ValueError):
            self.data.crop(17, 'bottom', dimension='pretrig')
        self.data.crop(6, 'top', dimension='pretrig')
        self.assertTrue(self.data.data.shape == (18, 40))

    def test_CropTrigMat(self):
        self.data.trig = np.ones((40,), dtype=int)
        self.data.trig[20:] = 2
        self.data.crop(6, 'top', dimension='pretrig')
        self.assertTrue(self.data.data.shape == (19, 40))

    def test_CropDepthOnTheFly(self):
        self.data.crop(0.165, 'bottom', dimension='depth', uice=2.0e6)
        self.assertTrue(self.data.data.shape == (17, 40))
        self.data.crop(0.055, 'top', dimension='depth', uice=2.0e6)
        self.assertTrue(self.data.data.shape == (11, 40))

    def test_CropDepthWithNMO(self):
        self.data.nmo(0., uice=2.0e6, uair=2.0e6)
        self.data.crop(0.165, 'bottom', dimension='depth')
        self.assertTrue(self.data.data.shape == (17, 40))
        self.data.crop(0.055, 'top', dimension='depth')
        self.assertTrue(self.data.data.shape == (11, 40))

    def test_HCropTnum(self):
        self.data.hcrop(2, 'left', dimension='tnum')
        self.assertTrue(self.data.data.shape == (20, 39))
        self.data.hcrop(15, 'right', dimension='tnum')
        self.assertTrue(self.data.data.shape == (20, 14))
        # Make sure we can ditch the last one
        self.data.hcrop(14, 'right', dimension='tnum')
        self.assertTrue(self.data.data.shape == (20, 13))

    def test_HCropInputErrors(self):
        with self.assertRaises(ValueError):
            self.data.hcrop(2, 'left', dimension='dummy')
        with self.assertRaises(ValueError):
            self.data.hcrop(2, 'dummy', dimension='tnum')

    def test_HCropBoundsErrors(self):
        # There are lots of bad inputs for tnum
        with self.assertRaises(ValueError):
            self.data.hcrop(44, 'right', dimension='tnum')
        with self.assertRaises(ValueError):
            self.data.hcrop(-44, 'right', dimension='tnum')
        with self.assertRaises(ValueError):
            self.data.hcrop(0, 'right', dimension='tnum')
        with self.assertRaises(ValueError):
            self.data.hcrop(1, 'right', dimension='tnum')
        with self.assertRaises(ValueError):
            self.data.hcrop(-1, 'right', dimension='tnum')
        with self.assertRaises(ValueError):
            self.data.hcrop(41, 'right', dimension='tnum')

        # Fewer ways to screw up distance
        with self.assertRaises(ValueError):
            self.data.hcrop(1.6, 'right', dimension='dist')
        with self.assertRaises(ValueError):
            self.data.hcrop(0, 'right', dimension='dist')
        with self.assertRaises(ValueError):
            self.data.hcrop(-1, 'right', dimension='dist')

    def test_HCropDist(self):
        self.data.hcrop(0.01, 'left', dimension='dist')
        self.assertTrue(self.data.data.shape == (20, 39))
        self.data.hcrop(1.4, 'right', dimension='dist')
        self.assertTrue(self.data.data.shape == (20, 38))

    def test_agc(self):
        self.data.agc()
        self.assertTrue(self.data.flags.agc)

    def test_rangegain(self):
        self.data.rangegain(1.0)
        self.assertTrue(self.data.flags.rgain)
        self.data.flags.rgain = False

        self.data.trig = np.zeros((self.data.tnum, ))
        self.data.rangegain(1.0)
        self.assertTrue(self.data.flags.rgain)

    def test_NMO_noexcpetion(self):
        # If velocity is 2
        self.data.nmo(0., uice=2.0, uair=2.0)
        self.assertTrue(np.allclose(self.data.travel_time * 1.0e-6, self.data.nmo_depth))
        # shouldn't care about uair if offset=0
        self.data.nmo(0., uice=2.0, uair=200.0)
        self.assertTrue(np.allclose(self.data.travel_time * 1.0e-6, self.data.nmo_depth))

        self.data.flags.nmo = False
        self.data.nmo(0., uice=2.0, uair=200.0)
        self.assertEqual(self.data.flags.nmo.shape, (2,))
        self.assertTrue(self.data.flags.nmo[0])

    def test_restack_odd(self):
        self.data.restack(5)
        self.assertTrue(self.data.data.shape == (20, 8))

    def test_restack_even(self):
        self.data.restack(4)
        self.assertTrue(self.data.data.shape == (20, 8))

    def test_elev_correct(self):
        self.data.elev = np.arange(self.data.data.shape[1]) * 0.002
        with self.assertRaises(ValueError):
            self.data.elev_correct()
        self.data.nmo(0, 2.0e6)
        self.data.elev_correct(v_avg=2.0e6)
        self.assertTrue(self.data.data.shape == (27, 40))

    def test_constant_space(self):
        distlims = (self.data.dist[0], self.data.dist[-1])
        space = 100.
        self.data.constant_space(space)
        self.assertTrue(self.data.data.shape == (20, np.ceil((distlims[-1] - distlims[0]) * 1000. / space)))
        self.assertTrue(self.data.x_coord.shape == (np.ceil((distlims[-1] - distlims[0]) * 1000. / space), ))
        self.assertTrue(self.data.y_coord.shape == (np.ceil((distlims[-1] - distlims[0]) * 1000. / space), ))
        self.assertTrue(self.data.lat.shape == (np.ceil((distlims[-1] - distlims[0]) * 1000. / space), ))
        self.assertTrue(self.data.long.shape == (np.ceil((distlims[-1] - distlims[0]) * 1000. / space), ))
        self.assertTrue(self.data.elev.shape == (np.ceil((distlims[-1] - distlims[0]) * 1000. / space), ))
        self.assertTrue(self.data.decday.shape == (np.ceil((distlims[-1] - distlims[0]) * 1000. / space), ))

        # do not fail because flags structure is weird from matlab
        space = 100.
        self.data.flags.interp = False
        self.data.constant_space(space)
        self.assertTrue(self.data.flags.interp.shape == (2,))
        self.assertTrue(self.data.flags.interp[0])
        self.assertEqual(self.data.flags.interp[1], space)

    def test_constant_sample_depth_spacing(self):
        # first check that it fails if we are not set up
        self.data.nmo_depth = None
        with self.assertRaises(AttributeError):
            self.data.constant_sample_depth_spacing()

        # Spoof variable nmo depths
        self.data.nmo_depth = np.hstack((np.arange(self.data.snum // 2),
                                         self.data.snum // 2 + 2. * np.arange(self.data.snum // 2)))
        self.data.constant_sample_depth_spacing()
        self.assertTrue(np.allclose(np.diff(self.data.nmo_depth), np.ones((self.data.snum - 1,)) * np.diff(self.data.nmo_depth)[0]))

        # So now if we call again, it should do nothing and return 1
        rv = self.data.constant_sample_depth_spacing()
        self.assertEqual(1, rv)

    def test_traveltime_to_depth(self):
        # We are not constant
        depths = self.data.traveltime_to_depth(np.arange(10) - 1., (np.arange(10) + 1) * 91.7)
        self.assertFalse(np.allclose(np.diff(depths), np.ones((len(depths) - 1,)) * (depths[1] - depths[0])))

        # We are constant
        depths = self.data.traveltime_to_depth(np.arange(10) - 1., (np.ones((10,)) * 91.7))
        self.assertTrue(np.allclose(np.diff(depths), np.ones((len(depths) - 1,)) * (depths[1] - depths[0])))

        # we have negative travel times
        self.data.travel_time = self.data.travel_time - 0.01
        depths = self.data.traveltime_to_depth(np.arange(10) - 1., (np.arange(10) + 1) * 91.7)
        self.assertFalse(np.allclose(np.diff(depths), np.ones((len(depths) - 1,)) * (depths[1] - depths[0])))

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_out.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))


if __name__ == '__main__':
    unittest.main()
