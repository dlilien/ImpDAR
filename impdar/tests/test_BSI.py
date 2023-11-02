#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make sure that we can successfully read BSI input files
"""

import os
import unittest
from impdar.lib.load import load_bsi

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestBSI(unittest.TestCase):

    @unittest.skipIf(not load_bsi.H5, 'h5py is not available')
    def test_load_bsi_pre2023(self):
        load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'test_bsi.h5'), XIPR=False)
        load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'test_bsi_xipr.h5'), XIPR=True)

    @unittest.skipIf(not load_bsi.H5, 'h5py is not available')
    def test_load_bsi_2023(self):
        load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'bsi_2023.h5'), XIPR=False)
        load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'bsi_2023.h5'), XIPR=True)
        load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'bsi_2023.h5'), nans="interp")
        load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'bsi_2023.h5'), nans="delete")

    @unittest.skipIf(load_bsi.H5, 'h5py is available')
    def test_load_bsi_noh5py(self):
        with self.assertRaises(ImportError):
            load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'test_bsi.h5'), XIPR=False)

        with self.assertRaises(ImportError):
            load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'test_bsi_xipr.h5'), XIPR=True)


if __name__ == '__main__':
    unittest.main()
