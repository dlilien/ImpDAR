#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Test the methods on the Pick object
"""

import os
import unittest
import numpy as np
from impdar.lib.RadarData import RadarData

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestPickMods(unittest.TestCase):

    def test_add_pick_loaded(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data.picks.add_pick(2)
        self.assertTrue(data.picks.samp1.shape == (3, data.tnum))
        data.picks.samp1[-1, :] = 0
        data.picks.samp2[-1, :] = 0
        data.picks.samp3[-1, :] = 0
        data.picks.add_pick(10)
        self.assertTrue(data.picks.samp1.shape == (4, data.tnum))

    def test_add_pick_blank(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        self.assertTrue(data.picks.samp1.shape == (1, data.tnum))

    def test_add_pick_badpicknum(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        # need to do this to prevent overwriting
        data.picks.samp1[0, 0] = 1.
        with self.assertRaises(ValueError):
            data.picks.add_pick(1)

    def test_add_pick_overwrite(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        data.picks.add_pick(2)
        self.assertTrue(data.picks.samp1.shape == (1, data.tnum))

    def test_update_pick(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        data.picks.update_pick(1, np.zeros((5, data.tnum)))
        self.assertTrue(np.all(data.picks.samp1 == 0))

    def test_update_pick_badpick_infoshape(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        with self.assertRaises(ValueError):
            data.picks.update_pick(1, np.zeros((4, 2)))

    def test_update_pick_badpicknum(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        with self.assertRaises(ValueError):
            data.picks.update_pick(0, np.zeros((5, data.tnum)))

    def test_smooth(self):
        # first, no NaNs
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        cache_val = data.picks.samp1.copy()
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            val = getattr(data.picks, attr)
            val[np.isnan(val)] = 1
            setattr(data.picks, attr, val)
        data.picks.smooth(4, units='tnum')

        # NaNs ends only
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            val = getattr(data.picks, attr)
            val[:, -1] = np.NaN
            val[:, 0] = np.NaN
            setattr(data.picks, attr, val)
        data.picks.smooth(4, units='tnum')

        # Middle Nans
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            val = getattr(data.picks, attr)
            val[:, -1] = np.NaN
            val[:, 0] = np.NaN
            val[:, 5] = np.NaN
            setattr(data.picks, attr, val)
        data.picks.smooth(4, units='tnum')

        # one row all nans
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            val = getattr(data.picks, attr)
            val[0, :] = np.NaN
            setattr(data.picks, attr, val)
        data.picks.smooth(4, units='tnum')

        # Now with dist
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data.flags.interp = [2, 1]
        for attr in ['samp1', 'samp2', 'samp3', 'power']:
            val = getattr(data.picks, attr)
            val[:, -1] = np.NaN
            val[:, 0] = np.NaN
            val[:, 5] = np.NaN
            setattr(data.picks, attr, val)
        data.picks.smooth(4, units='dist')

        # do not complain if nothing to do
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data.picks.samp1 = None
        data.picks.smooth(4, units='tnum')

        # fail with dist but no interp
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data.flags.interp = None
        with self.assertRaises(Exception):
            data.picks.smooth(4, units='dist')

        # Fail with elevation
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data.flags.elev = True
        with self.assertRaises(Exception):
            data.picks.smooth(4, units='dist')
        with self.assertRaises(Exception):
            data.picks.smooth(4, units='tnum')

        # Fail with bad units
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data.flags.interp = [2, 1]
        with self.assertRaises(ValueError):
            data.picks.smooth(4, 'dum')

        # Now make sure we fail with bad wavelengths--too high or too low for both units
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data.flags.interp = [2, 1]
        with self.assertRaises(ValueError):
            data.picks.smooth(0.5, 'tnum')
        with self.assertRaises(ValueError):
            data.picks.smooth(data.flags.interp[0] / 2, 'dist')
        with self.assertRaises(ValueError):
            data.picks.smooth(data.tnum + 2, 'tnum')
        with self.assertRaises(ValueError):
            data.picks.smooth(data.flags.interp[0] * data.tnum + 2, 'dist')


if __name__ == '__main__':
    unittest.main()
