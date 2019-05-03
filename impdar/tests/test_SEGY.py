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
try:
    from impdar.lib.load_segy import load_segy
    segy = True
except ImportError:
    segy = False

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class NoInitRadarData(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        self.data = np.array([[2, 2], [1, 1]])
        # need to set this to avoid divide by zero later
        self.elevation = np.zeros((2,))
        self.lat = np.ones((2,)) * 89.
        self.long = np.ones((2,))
        self.dt = 0.1
        self.tnum = self.data.shape[1]
        self.snum = self.data.shape[0]
        self.travel_time = 0.001 * np.arange(self.data.shape[0]) + 0.001
        self.dt = 1


class TestSEGY(unittest.TestCase):

    @unittest.skipIf(not segy, 'No SEGY on this version')
    def test_ReadSucceeds(self):
        data = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))

    @unittest.skipIf(not segy, 'No SEGY on this version')
    def test_WriteSucceeds(self):
        data = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))
        data.save_as_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))

    @unittest.skipIf(not segy, 'No SEGY on this version')
    def test_ReadWriteRead(self):
        data = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))
        data.save_as_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))
        data2 = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))
        self.assertEqual(data.data.shape, data2.data.shape)

    @unittest.skipIf(segy, 'SEGY on this version, only a graceful failure test')
    def test_SaveFails(self):
        data = NoInitRadarData()
        with self.assertRaises(ImportError):
            data.save_as_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))


if __name__ == '__main__':
    unittest.main()
