#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import sys
import unittest
import numpy as np
from impdar.lib.NoInitRadarData import NoInitRadarDataFiltering as NoInitRadarData
from impdar.lib.RadarData import RadarData
from impdar.lib import process
from impdar.lib.ImpdarError import ImpdarError
if sys.version_info[0] >= 3:
    from unittest.mock import MagicMock, patch
else:
    from mock import MagicMock, patch

data_dummy = np.ones((500, 400))


def Any(cls):
    # to mock data argument in tests
    class Any(cls):
        def __init__(self):
            pass

        def __eq__(self, other):
            return True
    return Any()


class TestAdaptive(unittest.TestCase):

    def test_AdaptiveRun(self):
        radardata = NoInitRadarData()
        radardata.adaptivehfilt(window_size=radardata.tnum // 10)
        self.assertTrue(np.all(radardata.data <= 1.))

        # make sure it works with a big window too
        radardata = NoInitRadarData()
        radardata.adaptivehfilt(window_size=radardata.tnum * 2)
        self.assertTrue(np.all(radardata.data <= 1.))


class TestHfilt(unittest.TestCase):

    def test_horizontalfilt(self):
        radardata = NoInitRadarData()
        radardata.horizontalfilt(0, 100)
        # We taper in the hfilt, so this is not just zeros
        self.assertTrue(np.all(radardata.data == radardata.hfilt_target_output))


class TestHighPass(unittest.TestCase):

    def test_highpass_simple(self):
        radardata = NoInitRadarData()

        # fails without constant-spaced data
        radardata.flags.interp = np.ones((2,))
        radardata.highpass(radardata.tnum * radardata.flags.interp[1] * 0.8)
        # There is no high-frequency variability, so this result should be small
        # We only have residual variability from the quality of the filter
        self.assertTrue(np.all(np.abs((radardata.data - radardata.data[0, 0])) < 1.0e-3))

    def test_highpass_badcutoff(self):
        radardata = NoInitRadarData()

        # fails without constant-spaced data
        radardata.flags.interp = np.ones((2,))

        with self.assertRaises(ValueError):
            radardata.highpass(radardata.flags.interp[1] * 0.5)
        with self.assertRaises(ValueError):
            radardata.highpass(radardata.tnum * radardata.flags.interp[1] * 1.5)

    def test_highpass_errors(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ImpdarError):
            radardata.highpass(100.0)

        # Elevation corrected data should fail
        radardata.flags.interp = np.ones((2,))
        # make sure this throws no error, then
        radardata.highpass(100.0)
        with self.assertRaises(ImpdarError):
            radardata.flags.elev = True
            radardata.highpass(100.0)


class TestHorizontalBandPass(unittest.TestCase):

    def test_hbp_simple(self):
        radardata = NoInitRadarData()

        # fails without constant-spaced data, so pretend it is constant spaced
        radardata.flags.interp = np.ones((2,))
        radardata.horizontal_band_pass(5., radardata.tnum * radardata.flags.interp[1] * 0.9)
        # We cannot really check this since the filter causes some residual variability as an edge effect

    def test_hbp_badcutoff(self):
        radardata = NoInitRadarData()

        # fails without constant-spaced data
        radardata.flags.interp = np.ones((2,))
        with self.assertRaises(ValueError):
            radardata.horizontal_band_pass(0.5, radardata.tnum / 10.)

        with self.assertRaises(ValueError):
            radardata.horizontal_band_pass(radardata.tnum / 10., radardata.tnum * 2.)

    def test_hbp_errors(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ImpdarError):
            # We have a screwed up filter here because of sampling vs. frequency used
            radardata.horizontal_band_pass(1000.0, 2000.0)

        radardata.flags.interp = np.ones((2,))
        # make sure this throws no error, then
        radardata.horizontal_band_pass(radardata.tnum / 10., radardata.tnum / 2.)
        with self.assertRaises(ValueError):
            radardata.horizontal_band_pass(radardata.tnum / 2., radardata.tnum / 10.)

        # Elevation corrected data should fail
        with self.assertRaises(ImpdarError):
            radardata.flags.elev = True
            radardata.horizontal_band_pass(radardata.tnum / 10., radardata.tnum / 2.)


class TestLowPass(unittest.TestCase):

    def test_lowpass_simple(self):
        radardata = NoInitRadarData()

        # fails without constant-spaced data
        radardata.flags.interp = np.ones((2,))
        radardata.lowpass(100.0)
        # There is no high-frequency variability, so this result should be small
        # We only have residual variability from the quality of the filter
        self.assertTrue(np.all(np.abs((radardata.data - radardata.data[0, 0]) / radardata.data[0, 0]) < 1.0e-3))

    def test_lowpass_badcutoff(self):
        # We have a screwed up filter here because of sampling vs. frequency used
        radardata = NoInitRadarData()
        # fails without constant-spaced data
        radardata.flags.interp = np.ones((2,))

        with self.assertRaises(ValueError):
            radardata.lowpass(radardata.flags.interp[1] * 0.5)
        with self.assertRaises(ValueError):
            radardata.lowpass(radardata.tnum * 1.5)

    def test_lowpass_errors(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ImpdarError):
            # We have a screwed up filter here because of sampling vs. frequency used
            radardata.lowpass(100.0)

        # Elevation corrected data should fail
        radardata.flags.interp = np.ones((2,))
        # make sure this throws no error, then
        radardata.lowpass(100.0)
        with self.assertRaises(ImpdarError):
            radardata.flags.elev = True
            radardata.lowpass(100.0)


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
        self.assertTrue(np.all(np.abs(radardata.data) < 1.0e-1))

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


class TestDenoise(unittest.TestCase):

    def test_denoise(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ValueError):
            radardata.denoise()
        radardata.data = np.random.random(radardata.data.shape)
        radardata.denoise()

        radardata = NoInitRadarData()
        radardata.data = np.random.random(radardata.data.shape)
        radardata.denoise(noise=0.1)

    def test_denoise_badftype(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ValueError):
            radardata.denoise(ftype='dummy')


class TestRadarDataHfiltWrapper(unittest.TestCase):

    def test_adaptive(self):
        radardata = NoInitRadarData()
        radardata.adaptivehfilt = MagicMock()
        radardata.hfilt(ftype='adaptive', window_size=1000)
        radardata.adaptivehfilt.assert_called_with(window_size=1000)

    def test_horizontalfilt(self):
        radardata = NoInitRadarData()
        radardata.horizontalfilt = MagicMock()
        radardata.hfilt(ftype='hfilt', bounds=(0, 100))
        radardata.horizontalfilt.assert_called_with(0, 100)

    def test_badfilter(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ValueError):
            radardata.hfilt(ftype='dummy')


class TestProcessWrapper(unittest.TestCase):

    def test_process_ahfilt(self):
        radardata = NoInitRadarData()
        radardata.adaptivehfilt = MagicMock()
        process.process([radardata], ahfilt=1000)
        radardata.adaptivehfilt.assert_called_with(window_size=1000)

    def test_process_hfilt(self):
        radardata = NoInitRadarData()
        radardata.horizontalfilt = MagicMock()
        process.process([radardata], hfilt=(0, 100))
        radardata.horizontalfilt.assert_called_with(0, 100)

    def test_process_vbp(self):
        radardata = NoInitRadarData()
        radardata.vertical_band_pass = MagicMock()
        process.process([radardata], vbp=(0.1, 100.))
        # The filter is not too good, so we have lots of residual
        radardata.vertical_band_pass.assert_called_with(0.1, 100.)


class TestMigrationWrapper(unittest.TestCase):
    """This is only to make sure the calls are setup correctly. Actual tests are separate"""

    @patch('impdar.lib.migrationlib.migrationKirchhoff')
    def test_wrap_kirchhoff(self, patch_ob):
        radardata = NoInitRadarData()
        radardata.migrate(mtype='kirch', vel=10., nearfield=False)
        patch_ob.assert_called_with(Any(RadarData), vel=10., nearfield=False)

    @patch('impdar.lib.migrationlib.migrationStolt')
    def test_wrap_stolt(self, patch_ob):
        radardata = NoInitRadarData()
        radardata.migrate(mtype='stolt', htaper=1, vtaper=2, vel=999.)
        patch_ob.assert_called_with(Any(RadarData), htaper=1, vtaper=2, vel=999.)

    @patch('impdar.lib.migrationlib.migrationPhaseShift')
    def test_wrap_phaseshift(self, patch_ob):
        radardata = NoInitRadarData()
        radardata.migrate(mtype='phsh', vel=1., vel_fn='dummy', htaper=1, vtaper=2)
        patch_ob.assert_called_with(Any(RadarData), vel=1., vel_fn='dummy', htaper=1, vtaper=2)

    @patch('impdar.lib.migrationlib.migrationTimeWavenumber')
    def test_wrap_tk(self, patch_ob):
        radardata = NoInitRadarData()
        radardata.migrate(mtype='tk', vel=1., vel_fn='dummy', htaper=1, vtaper=2)
        patch_ob.assert_called_with(Any(RadarData), vel=1., vel_fn='dummy', htaper=1, vtaper=2)

    @patch('impdar.lib.migrationlib.migrationSeisUnix')
    def test_wrap_seisunix(self, patch_ob):
        radardata = NoInitRadarData()
        radardata.migrate(mtype='su_stolt', vtaper=1, htaper=2, tmig=3, vel_fn=None, vel=1.68e7, nxpad=15, verbose=1)
        patch_ob.assert_called_with(Any(RadarData), vtaper=1, htaper=2, tmig=3, vel_fn=None, vel=1.68e7, nxpad=15, verbose=1, mtype='su_stolt')

    def test_bad_mtype(self):
        radardata = NoInitRadarData()
        with self.assertRaises(ValueError):
            radardata.migrate(mtype='dummy')

        # and spoof the checker for seisunix
        with self.assertRaises(Exception):
            radardata.migrate(mtype='su_dummy')



if __name__ == '__main__':
    unittest.main()
