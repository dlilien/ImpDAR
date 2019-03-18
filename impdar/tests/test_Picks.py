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
from impdar.lib import load

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestLoad(unittest.TestCase):

    def test_load_mat(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        self.assertEqual(data.data.shape, (20, 40))


class TestPickMods(unittest.TestCase):

    def test_add_pick_loaded(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        data.picks.add_pick(1)
        self.assertTrue(data.picks.samp1.shape == (3, data.tnum))

    def test_add_pick_blank(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        self.assertTrue(data.picks.samp1.shape == (1, data.tnum))

    def test_add_pick_badpicknum(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        # need to do this to prevent overwriting
        data.picks.samp1[0, 0] = 1.
        with self.assertRaises(ValueError):
            data.picks.add_pick(1)

    def test_add_pick_overwrite(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        data.picks.add_pick(2)
        self.assertTrue(data.picks.samp1.shape == (1, data.tnum))

    def test_update_pick(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        data.picks.update_pick(1, np.zeros((5, data.tnum)))
        self.assertTrue(np.all(data.picks.samp1 == 0))

    def test_update_pick_badpick_infoshape(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        with self.assertRaises(ValueError):
            data.picks.update_pick(1, np.zeros((4, 2)))

    def test_update_pick_badpicknum(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        data.picks.add_pick(1)
        with self.assertRaises(ValueError):
            data.picks.update_pick(0, np.zeros((5, data.tnum)))


if __name__ == '__main__':
    unittest.main()
