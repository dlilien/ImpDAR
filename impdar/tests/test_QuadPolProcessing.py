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
from impdar.lib.ApresData import ApresData, ApresQuadPol
from impdar.lib.ApresData.load_apres import load_apres
from impdar.lib.ApresData.load_quadpol import load_quadpol, load_quadpol_fujita
from impdar.lib.ApresData._QuadPolProcessing import *

data_dummy = np.ones((500, 400))

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestProcessing(unittest.TestCase):

    def test_1_load_save(self):
        qpdat = load_quadpol(os.path.join(THIS_DIR, 'input_data', 'quadpol'))
        self.assertTrue(hasattr(qpdat,'shh'))
        self.assertTrue(hasattr(qpdat,'shv'))
        qpdat.save(os.path.join(THIS_DIR, 'input_data', 'qpdat.mat'))

    def test_2_load_fujita(self):
        data = load_quadpol_fujita(os.path.join(THIS_DIR, 'input_data', 'quadpol_fujita.mat'))
        self.assertTrue(data.data_dtype == np.cdouble().dtype)

    def test_3_rotation(self):
        qpdat = ApresQuadPol(os.path.join(THIS_DIR, 'input_data', 'qpdat.mat'))
        with self.assertRaises(ValueError):
            qpdat.rotational_transform()
        qpdat.rotational_transform(cross_pol_flip='HV')
        self.assertTrue(hasattr(qpdat,'HH'))
        self.assertTrue(hasattr(qpdat,'HV'))

    def test_4_find_cpe(self):
        qpdat = ApresQuadPol(os.path.join(THIS_DIR, 'input_data', 'qpdat.mat'))
        with self.assertRaises(ImpdarError):
            qpdat.find_cpe()
        qpdat.rotational_transform(cross_pol_flip='HV')
        qpdat.find_cpe()
        self.assertTrue(np.shape(qpdat.cpe) == (qpdat.snum,))

    def test_5_coherence(self):
        qpdat = ApresQuadPol(os.path.join(THIS_DIR, 'input_data', 'qpdat.mat'))
        with self.assertRaises(ImpdarError):
            qpdat.coherence2d()
        qpdat.rotational_transform(cross_pol_flip='HV',n_thetas=5)
        qpdat.coherence2d(delta_range=10, delta_theta=60.*np.pi/180.)
        self.assertTrue(qpdat.chhvv.dtype == complex)
        chhvv_hold = qpdat.chhvv.copy()
        qpdat.coherence2d(delta_range=10, delta_theta=60.*np.pi/180., force_python=True)
        self.assertTrue(qpdat.chhvv.dtype == complex)
        self.assertTrue(np.all(abs(qpdat.chhvv - chhvv_hold) < 1e-5))
        qpdat.find_cpe()
        qpdat.save(os.path.join(THIS_DIR, 'input_data', 'qpdat_coh.mat'))

    def test_6_phase_gradient(self):
        qpdat = ApresQuadPol(os.path.join(THIS_DIR, 'input_data', 'qpdat_coh.mat'))
        qpdat.phase_gradient2d(filt='lowpass', Wn=50)
        self.assertTrue(qpdat.dphi_dz.dtype == float)

if __name__ == '__main__':
    unittest.main()
