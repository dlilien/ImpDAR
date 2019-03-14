#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the basics of RadarData
"""
import os
import unittest
import numpy as np
try:
    from impdar.lib.load_segy import load_segy
    segy = True
except ImportError:
    segy = False

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


@unittest.skipIf(not segy, 'No SEGY on this version')
class TestSEGY(unittest.TestCase):

    def test_ReadSucceeds(self):
        data = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))

    def test_WriteSucceeds(self):
        data = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))
        data.save_as_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))

    def test_ReadWriteRead(self):
        data = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'))
        data.save_as_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))
        data2 = load_segy(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))
        self.assertEqual(data.data.shape, data2.data.shape)

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy')):
            pass
            # os.remove(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200_resave.segy'))


if __name__ == '__main__':
    unittest.main()
