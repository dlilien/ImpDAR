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
from impdar.lib.ApresData import ApresData, ApresTimeDiff
from impdar.lib.ApresData.load_apres import load_apres
from impdar.lib.ApresData.load_time_diff import load_time_diff
from impdar.lib.ApresData._TimeDiffProcessing import *

data_dummy = np.ones((500, 400))

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestProcessing(unittest.TestCase):

    def test_diff_load_save(self):
        apresdata_1 = load_apres([os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT')])
        apresdata_1.apres_range(2)
        apresdata_1.stacking()
        apresdata_1.phase_uncertainty(3000)
        apresdata_2 = load_apres([os.path.join(THIS_DIR, 'input_data', 'apres_2.DAT')])
        apresdata_2.apres_range(2)
        apresdata_2.stacking()
        apresdata_2.phase_uncertainty(3000)
        diffdat = load_time_diff([apresdata_1,apresdata_2])
        self.assertTrue(hasattr(diffdat,'data'))
        self.assertTrue(hasattr(diffdat,'data2'))
        diffdat.save(os.path.join(THIS_DIR, 'input_data', 'diffdat.mat'))

    def test_phase_differencing(self):
        diffdat = ApresTimeDiff(os.path.join(THIS_DIR, 'input_data', 'diffdat.mat'))
        win = 20
        step = 20
        diffdat.phase_diff(win, step)
        self.assertTrue(diffdat.co.dtype == np.complex128)
        self.assertTrue(max(abs(diffdat.co))<=1.)
        self.assertTrue(min(abs(diffdat.co))>=0.)

    def test_unwrap(self):
        diffdat = ApresTimeDiff(os.path.join(THIS_DIR, 'input_data', 'diffdat.mat'))
        diffdat.phase_diff(20, 20)
        diffdat.phase_unwrap()
        self.assertTrue(diffdat.phi.dtype == np.float64)

    def test_range_diff(self):
        diffdat = ApresTimeDiff(os.path.join(THIS_DIR, 'input_data', 'diffdat.mat'))
        diffdat.phase_diff(20, 20)
        diffdat.phase_unwrap()
        diffdat.range_diff()
        diffdat.range_diff(uncertainty='noise_phasor')
        self.assertTrue(max(abs(diffdat.w))<100.)

    def test_strain_rate(self):
        diffdat = ApresTimeDiff(os.path.join(THIS_DIR, 'input_data', 'diffdat.mat'))
        diffdat.phase_diff(20, 20)
        diffdat.phase_unwrap()
        diffdat.range_diff()
        diffdat.strain_rate((200,1000))

    def test_pick_bed(self):
        diffdat = ApresTimeDiff(os.path.join(THIS_DIR, 'input_data', 'diffdat.mat'))
        diffdat.phase_diff(20, 20)
        diffdat.bed_pick()
        self.assertTrue(np.shape(diffdat.bed)==(4,))

if __name__ == '__main__':
    unittest.main()
