#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the machinery of process. This is broken up to match where it would likely fail; tests process wrappers of various methods are with the tests of those methods
"""
import sys
import os
import unittest
import numpy as np
from impdar.lib.NoInitRadarData import NoInitRadarData
from impdar.lib.RadarData import RadarData
from impdar.lib import process
if sys.version_info[0] >= 3:
    from unittest.mock import MagicMock, patch
else:
    from mock import MagicMock, patch

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestConcat(unittest.TestCase):

    def test_concat_nopicks(self):
        dats = process.concat([NoInitRadarData(), NoInitRadarData()])
        self.assertTrue(dats[0].data.shape == (2, 4))

        with self.assertRaises(ValueError):
            d2 = NoInitRadarData()
            d2.snum = 3
            dats = process.concat([NoInitRadarData(), d2])

        with self.assertRaises(ValueError):
            d2 = NoInitRadarData()
            d2.travel_time = np.array((2, 3))
            dats = process.concat([NoInitRadarData(), d2])

    def test_concat_picks(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))

        # no overlapping picks
        data_otherp = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data_otherp.picks.picknums = [pn * 10 - 1 for pn in data_otherp.picks.picknums]

        # one overlapping pick
        data_somepsame = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data_somepsame.picks.picknums = [1, 19]


        dats = process.concat([data, data])
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            self.assertEqual(getattr(dats[0].picks, attr).shape[1], 2 * getattr(data.picks, attr).shape[1])
            self.assertEqual(getattr(dats[0].picks, attr).shape[0], getattr(data.picks, attr).shape[0])

        dats = process.concat([data, data_otherp])
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            self.assertTrue(getattr(dats[0].picks, attr).shape[1] == 2 * data.picks.samp1.shape[1])
            self.assertTrue(getattr(dats[0].picks, attr).shape[0] == 2 * data.picks.samp1.shape[0])
            for pn in data.picks.picknums:
                self.assertTrue(pn in dats[0].picks.picknums)
            for pn in data_otherp.picks.picknums:
                self.assertTrue(pn in dats[0].picks.picknums)

        dats = process.concat([data, data_somepsame])
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            self.assertTrue(getattr(dats[0].picks, attr).shape[1] == 2 * data.picks.samp1.shape[1])
            self.assertTrue(getattr(dats[0].picks, attr).shape[0] == 2 * data.picks.samp1.shape[0] - 1)
            for pn in data.picks.picknums:
                self.assertTrue(pn in dats[0].picks.picknums)
            for pn in data_somepsame.picks.picknums:
                self.assertTrue(pn in dats[0].picks.picknums)

        # no picks
        data_np = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data_np.picks.picknums = 0

        dats = process.concat([data, data_np])
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            self.assertTrue(getattr(dats[0].picks, attr).shape[1] == 2 * data.picks.samp1.shape[1])
            self.assertTrue(np.all(np.isnan(getattr(dats[0].picks, attr)[0, data.picks.samp1.shape[1]:])))
            for pn in data.picks.picknums:
                self.assertTrue(pn in dats[0].picks.picknums)

        data_np.picks.picknums = None
        dats = process.concat([data, data_np])
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            self.assertTrue(getattr(dats[0].picks, attr).shape[1] == 2 * data.picks.samp1.shape[1])
            self.assertTrue(np.all(np.isnan(getattr(dats[0].picks, attr)[0, data.picks.samp1.shape[1]:])))
            for pn in data.picks.picknums:
                self.assertTrue(pn in dats[0].picks.picknums)

        data_np.picks = None
        dats = process.concat([data, data_np])
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            self.assertTrue(getattr(dats[0].picks, attr).shape[1] == 2 * data.picks.samp1.shape[1])
            self.assertTrue(np.all(np.isnan(getattr(dats[0].picks, attr)[0, data.picks.samp1.shape[1]:])))
            for pn in data.picks.picknums:
                self.assertTrue(pn in dats[0].picks.picknums)


class TestProcess_and_exit(unittest.TestCase):

    def test_process_and_exitLOADMAT(self):
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')])

    def test_process_and_exitCAT(self):
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], cat=True)

    def test_process_and_exitPROCESS(self):
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], rev=True)

    def test_process_and_exitOUTNAMING(self):
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'data_raw.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], cat=True)
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'data_cat.mat')))
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], cat=True)
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data_cat.mat')))
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), os.path.join(THIS_DIR, 'input_data', 'data_raw.mat')], rev=True)
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data_proc.mat')))
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT')], filetype='gssi', rev=True)
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_gssi_proc.mat')))
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], cat=True, o=os.path.join(THIS_DIR, 'small_data_cat.mat'))
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'small_data_cat.mat')))
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'data_raw.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], cat=True, o=os.path.join(THIS_DIR, 'data_cat.mat'))
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'data_cat.mat')))
        process.process_and_exit([os.path.join(THIS_DIR, 'input_data', 'data_raw.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], rev=True, o=THIS_DIR)
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'data_proc.mat')))
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'small_data_proc.mat')))

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'small_data_cat.mat')):
            os.remove(os.path.join(THIS_DIR, 'small_data_cat.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data_cat.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'small_data_cat.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data_proc.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'small_data_proc.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'small_data_proc.mat')):
            os.remove(os.path.join(THIS_DIR, 'small_data_proc.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'data_proc.mat')):
            os.remove(os.path.join(THIS_DIR, 'data_proc.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'data_proc.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'data_proc.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'data_cat.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'data_cat.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'data_cat.mat')):
            os.remove(os.path.join(THIS_DIR, 'data_cat.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_gssi_proc.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'test_gssi_proc.mat'))


class TestProcess(unittest.TestCase):

    def setUp(self):
        self.data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.data.x_coord = np.arange(40)
        self.data.nmo_depth = None
        self.data.travel_time = np.arange(0, 0.2, 0.01)
        self.data.dt = 1.0e-8
        self.data.trig = self.data.trig * 0.

    def test_process_Reverse(self):
        self.data.reverse = MagicMock()
        self.assertTrue(process.process([self.data], rev=True))
        self.data.reverse.assert_called_with()

    def test_process_Crop(self):
        with self.assertRaises(TypeError):
            process.process([self.data], crop=True)
        with self.assertRaises(ValueError):
            process.process([self.data], crop=('ugachacka', 'top', 'snum'))
        self.data.crop = MagicMock()
        self.assertTrue(process.process([self.data], crop=(17, 'bottom', 'snum')))
        self.data.crop.assert_called_with(17, 'bottom', 'snum')

    def test_process_Denoise(self):
        self.data.denoise = MagicMock()
        with self.assertRaises(ValueError):
            process.process([self.data], denoise=1)
        with self.assertRaises(ValueError):
            process.process([self.data], denoise='12')
        with self.assertRaises(ValueError):
            process.process([self.data], denoise=(1, ))
        process.process([self.data], denoise=(1, 2))
        self.data.denoise.assert_called_with(1, 2)

    @patch('impdar.lib.process.interpdeep')
    def test_process_Interp(self, mock_interp):
        dl = [self.data]
        with self.assertRaises(ValueError):
            process.process([self.data], interp=('ba', 2))
        with self.assertRaises(ValueError):
            process.process([self.data], interp='ba')
        with self.assertRaises(ValueError):
            process.process([self.data], interp=1)
        with self.assertRaises(ValueError):
            process.process([self.data], interp=(1, ))
        process.process(dl, interp=(1, 2))
        mock_interp.assert_called_with(dl, 1.0, 2)

    def test_process_hcrop(self):
        with self.assertRaises(TypeError):
            process.process([self.data], hcrop=True)
        with self.assertRaises(ValueError):
            process.process([self.data], hcrop=('ugachacka', 'left', 'tnum'))
        self.data.hcrop = MagicMock()
        self.assertTrue(process.process([self.data], hcrop=(17, 'left', 'tnum')))
        self.data.hcrop.assert_called_with(17, 'left', 'tnum')

    def test_process_NMO(self):
        self.data.nmo = MagicMock()
        self.assertTrue(process.process([self.data], nmo=(0., 2.0, 2.0)))
        self.data.nmo.assert_called_with(0., 2.0, 2.0)

        self.data.nmo = MagicMock()
        self.assertTrue(process.process([self.data], nmo=0))
        self.data.nmo.assert_called_with(0, 1.6)

        self.data.nmo = MagicMock()
        self.assertTrue(process.process([self.data], nmo=1.0))
        self.data.nmo.assert_called_with(1.0, 1.6)

    def test_process_restack(self):
        self.data.restack = MagicMock()
        self.assertTrue(process.process([self.data], restack=3))
        self.data.restack.assert_called_with(3)

        self.data.restack = MagicMock()
        self.assertTrue(process.process([self.data], restack=[4., 'dummy']))
        self.data.restack.assert_called_with(4)

    def test_process_vbp(self):
        with self.assertRaises(TypeError):
            process.process([self.data], vbp=3)

        self.data.vertical_band_pass = MagicMock()
        self.assertTrue(process.process([self.data], vbp=(3, 4)))
        self.data.vertical_band_pass.assert_called_with(3, 4)

    def test_migrate(self):
        self.data.migrate = MagicMock()
        self.assertTrue(process.process([self.data], migrate=True))
        self.data.migrate.assert_called_with(mtype='stolt')

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_out.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))


if __name__ == '__main__':
    unittest.main()
