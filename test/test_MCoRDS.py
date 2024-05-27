#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make sure that we can successfully read gssi input files
"""

import os
import unittest
import numpy as np
from impdar.lib.load import load_mcords

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestMCoRDS_NC(unittest.TestCase):

    @unittest.skipIf(not load_mcords.NC, 'No netcdf on this version')
    def test_loadnc(self):
        dat = load_mcords.load_mcords_nc(os.path.join(THIS_DIR, 'input_data', 'zeros_mcords.nc'))
        self.assertTrue(np.all(dat.data == 0.))

    @unittest.skipIf(load_mcords.NC, 'NETCDF on this version')
    def test_loadnc_failure(self):
        with self.assertRaises(ImportError):
            load_mcords.load_mcords_nc(os.path.join(THIS_DIR, 'input_data', 'zeros_mcords.nc'))


class TestMCoRDS_MAT(unittest.TestCase):

    def test_loadmat(self):
        dat = load_mcords.load_mcords_mat(os.path.join(THIS_DIR,
                                                       'input_data',
                                                       'zeros_mcords_mat.mat'))
        self.assertTrue(np.allclose(dat.data, 0.))

    def test_loadbadmat(self):
        with self.assertRaises(KeyError):
            load_mcords.load_mcords_mat(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))

        with self.assertRaises(KeyError):
            load_mcords.load_mcords_mat(os.path.join(THIS_DIR,
                                                     'input_data',
                                                     'nonimpdar_matlab.mat'))


if __name__ == '__main__':
    unittest.main()
