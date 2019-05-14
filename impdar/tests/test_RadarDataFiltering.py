#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""

import os
import unittest
import numpy as np
from impdar.lib.RadarData import RadarData
from impdar.lib.RadarFlags import RadarFlags
from impdar.lib import process

data_dummy = np.ones((500, 400))


class NoInitRadarData(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        self.data = data_dummy.copy()
        self.dt = 0.1
        self.tnum = self.data.shape[1]
        self.snum = self.data.shape[0]
        self.travel_time = 0.001 * np.arange(self.data.shape[0]) + 0.001
        self.dt = 0.001e-6
        self.flags = RadarFlags()
        self.hfilt_target_output = data_dummy * np.atleast_2d(1. - np.exp(-self.travel_time.flatten() * 0.05) / np.exp(-self.travel_time[0] * 0.05)).transpose()
        pexp = np.exp(-self.travel_time.flatten() * 0.05) / np.exp(-self.travel_time[0] * 0.05)
        pexp = pexp - pexp[-1]
        pexp = pexp / np.max(pexp)
        self.pexp_target_output = data_dummy * np.atleast_2d(1. - pexp).transpose()
        self.ahfilt_target_output = np.zeros_like(data_dummy)


class TestAdaptive(unittest.TestCase):

    def test_AdaptiveRun(self):
        radardata = NoInitRadarData()
        radardata.adaptivehfilt()
        # since we subtract average trace and all traces are identical, we should get zeros out
        self.assertTrue(np.all(radardata.data == radardata.ahfilt_target_output))


class TestHfilt(unittest.TestCase):

    def test_HfiltRun(self):
        radardata = NoInitRadarData()
        radardata.horizontalfilt(0, 100)
        # We taper in the hfilt, so this is not just zeros
        self.assertTrue(np.all(radardata.data == radardata.hfilt_target_output))


class TestHighPass(unittest.TestCase):

    def test_HighPass(self):
        radardata = NoInitRadarData()
        radardata.highpass(1000.0, 1)
        # There is no high-frequency variability, so this result should be small
        # We only have residual variability from the quality of the filter
        print(np.abs((radardata.data - radardata.data[0, 0]) / radardata.data[0, 0]))
        self.assertTrue(np.all(np.abs((radardata.data - radardata.data[0, 0]) / radardata.data[0, 0]) < 1.0e-3))

    def test_HighPassBadcutoff(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ValueError):
            # We have a screwed up filter here because of sampling vs. frequency used
            radardata.highpass(1.0e-4, 100)


class TestWinAvgHfilt(unittest.TestCase):

    def test_WinAvgExp(self):
        radardata = NoInitRadarData()
        radardata.winavg_hfilt(11, taper='full')
        self.assertTrue(np.all(radardata.data == radardata.hfilt_target_output))

    def test_WinAvgExpBadwinavg(self):
        # Tests the check on whether win_avg < tnum
        radardata = NoInitRadarData()
        radardata.winavg_hfilt(data_dummy.shape[1] + 10, taper='full')
        self.assertTrue(np.all(radardata.data == radardata.hfilt_target_output))

    def test_WinAvgPexp(self):
        radardata = NoInitRadarData()
        radardata.winavg_hfilt(11, taper='pexp', filtdepth=-1)
        self.assertTrue(np.all(radardata.data == radardata.pexp_target_output))

    def test_WinAvgbadtaper(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ValueError):
            radardata.winavg_hfilt(11, taper='not_a_taper', filtdepth=-1)


class TestVBP(unittest.TestCase):

    def test_vbp_butter(self):
        radardata = NoInitRadarData()
        radardata.vertical_band_pass(0.1, 100., filttype='butter')
        # The filter is not too good, so we have lots of residual
        self.assertTrue(np.all(np.abs(radardata.data) < 1.0e-4))

    def test_vbp_cheb(self):
        radardata = NoInitRadarData()
        radardata.vertical_band_pass(0.1, 100., filttype='cheb')
        # The filter is not too good, so we have lots of residual
        self.assertTrue(np.all(np.abs(radardata.data) < 1.0e-2))

    def test_vbp_bessel(self):
        radardata = NoInitRadarData()
        radardata.vertical_band_pass(0.1, 100., filttype='bessel')
        # The filter is not too good, so we have lots of residual
        self.assertTrue(np.all(np.abs(radardata.data) < 1.0e-1))

    def test_vbp_fir(self):
        radardata = NoInitRadarData()
        radardata.vertical_band_pass(1., 10., filttype='fir', order=100)

        radardata.vertical_band_pass(1., 10., filttype='fir', order=2, fir_window='hanning')

    def test_vbp_badftype(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ValueError):
            radardata.vertical_band_pass(0.1, 100., filttype='dummy')


class TestRadarDataHfiltWrapper(unittest.TestCase):

    def test_AdaptiveRun(self):
        radardata = NoInitRadarData()
        radardata.hfilt('adaptive')
        # since we subtract average trace and all traces are identical, we should get zeros out
        self.assertTrue(np.all(radardata.data == radardata.ahfilt_target_output))

    def test_HfiltRun(self):
        radardata = NoInitRadarData()
        radardata.hfilt('hfilt', (0, 100))
        # We taper in the hfilt, so this is not just zeros
        self.assertTrue(np.all(radardata.data == radardata.hfilt_target_output))


class TestProcessWrapper(unittest.TestCase):
    def test_process_ahfilt(self):
        radardata = NoInitRadarData()
        process.process([radardata], ahfilt=True)
        # We taper in the hfilt, so this is not just zeros
        self.assertTrue(np.all(radardata.data == radardata.ahfilt_target_output))

    def test_process_hfilt(self):
        radardata = NoInitRadarData()
        process.process([radardata], hfilt=(0, 100))
        # We taper in the hfilt, so this is not just zeros
        self.assertTrue(np.all(radardata.data == radardata.hfilt_target_output))

    def test_process_vbp(self):
        radardata = NoInitRadarData()
        process.process([radardata], vbp=(0.1, 100.))
        # The filter is not too good, so we have lots of residual
        self.assertTrue(np.all(np.abs(radardata.data) < 1.0e-4))


if __name__ == '__main__':
    unittest.main()
