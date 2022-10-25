#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 bhills <benjaminhhills@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import sys
import os
import unittest
from impdar.lib.ApresData import load_apres

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

try:
    from netCDF4 import Dataset
    nc_load = True
except:
    nc_load = False

class TestLoad(unittest.TestCase):

    def test_loaddat(self):
        data = load_apres.load_apres_single_file(os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT'))
        self.assertEqual(data.data.shape, (data.header.n_subbursts, data.snum))

    def test_load2dat(self):
        data = load_apres.load_apres([os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT'),
                                      os.path.join(THIS_DIR, 'input_data', 'apres_2.DAT')])
        self.assertEqual(data.data.shape, (data.bnum, data.cnum, data.snum))

    def test_bas_nc(self):
        if nc_load:
            data = load_apres.load_BAS_nc(os.path.join(THIS_DIR, 'input_data', 'apres_1.nc'))
            self.assertEqual(data.data.shape, (data.bnum, data.cnum, data.snum))
        else:
            with self.assertRaises(ImportError):
                data = load_apres.load_BAS_nc(os.path.join(THIS_DIR, 'input_data', 'apres_1.nc'))

    def test_loadbad(self):
        with self.assertRaises(ValueError):
            data = load_apres.load_apres_single_file(os.path.join(THIS_DIR, 'input_data', 'small_data.bad'))

if __name__ == '__main__':
    unittest.main()
