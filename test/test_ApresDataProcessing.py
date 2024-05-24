#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Benjamin Hills <benjaminhhills@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import os
import sys
import unittest
import numpy as np
from impdar.lib.ApresData import ApresData
from impdar.lib.ApresData.load_apres import load_apres
from impdar.lib.ApresData._ApresDataProcessing import *

data_dummy = np.ones((500, 400))

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestProcessing(unittest.TestCase):

    def test_range_conversion(self):
        apresdata = load_apres([os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT')])
        apresdata.apres_range(2,4000)
        self.assertTrue(np.max(apresdata.Rcoarse)<=4000)
        self.assertTrue(apresdata.data.dtype == 'complex128')

    def test_stacking(self):
        apresdata = load_apres([os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT')])
        apresdata.stacking()
        self.assertTrue(np.shape(apresdata.data)==(1,1,apresdata.snum))

    def test_uncertainty_failrange(self):
        apresdata = load_apres([os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT')])
        with self.assertRaises(TypeError):
            apresdata.phase_uncertainty(3000)

    def test_uncertainty_stackrange(self):
        apresdata = load_apres([os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT')])
        apresdata.apres_range(2,4000)
        with self.assertRaises(IndexError):
            apresdata.phase_uncertainty(3000)

    def test_uncertainty_stackrange(self):
        apresdata = load_apres([os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT')])
        apresdata.apres_range(2,4000)
        apresdata.stacking()
        apresdata.phase_uncertainty(3000)
        self.assertTrue(apresdata.uncertainty.dtype == 'float64')
        self.assertTrue(len(apresdata.uncertainty) == apresdata.snum)

if __name__ == '__main__':
    unittest.main()
