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

import unittest
from impdar.lib import PickParameters, NoInitRadarData


class TestPickParameters(unittest.TestCase):

    def test_init(self):
        rd = NoInitRadarData.NoInitRadarData()
        pick_params = PickParameters.PickParameters(rd)
        for attr in pick_params.attrs:
            self.assertIsNotNone(getattr(pick_params, attr))

    def test_freq_update(self):
        rd = NoInitRadarData.NoInitRadarData()
        pick_params = PickParameters.PickParameters(rd)
        pick_params.freq_update(1000.0)
        self.assertEqual(pick_params.FWW, 1)
        self.assertEqual(pick_params.plength, 3)
        self.assertEqual(pick_params.scst, 1)

        rd = NoInitRadarData.NoInitRadarDataFiltering()
        pick_params = PickParameters.PickParameters(rd)
        pick_params.freq_update(1.0e-8)
        self.assertEqual(pick_params.plength, rd.snum)

    def test_to_struct(self):
        rd = NoInitRadarData.NoInitRadarData()
        pick_params = PickParameters.PickParameters(rd)
        mat = pick_params.to_struct()
        for attr in pick_params.attrs:
            self.assertIsNotNone(mat[attr])

        pick_params.dt = None
        mat = pick_params.to_struct()
        for attr in pick_params.attrs:
            self.assertIsNotNone(mat[attr])


if __name__ == '__main__':
    unittest.main()
