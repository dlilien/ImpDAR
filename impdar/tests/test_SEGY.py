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
import pytest
import numpy as np
from impdar.lib.NoInitRadarData import NoInitRadarData
from impdar.lib.load.load_segy import load_segy, SEGY

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestSEGY(unittest.TestCase):

    @unittest.skipIf(not SEGY, 'No SEGY on this version')
    def test_ReadSucceeds(self):
        pytest.importorskip('segyio')
        load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))

    @unittest.skipIf(not SEGY, 'No SEGY on this version')
    def test_WriteSucceeds(self):
        pytest.importorskip('segyio')
        data = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))
        data.save_as_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))

    @unittest.skipIf(not SEGY, 'No SEGY on this version')
    def test_ReadWriteRead(self):
        pytest.importorskip('segyio')
        data = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))
        data.save_as_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))
        data2 = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))
        self.assertEqual(data.data.shape, data2.data.shape)
        self.assertTrue(np.allclose(data.data, data2.data))

    @unittest.skipIf(SEGY, 'SEGY on this version, only a graceful failure test')
    def test_SaveFails(self):
        data = NoInitRadarData()
        with self.assertRaises(ImportError):
            data.save_as_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))


if __name__ == '__main__':
    unittest.main()
