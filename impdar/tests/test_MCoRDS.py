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
from impdar.lib.load import load_mcords_nc

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestMCoRDS(unittest.TestCase):

    @unittest.skipIf(not load_mcords_nc.NC, 'No netcdf on this version')
    def test_loadnc(self):
        dat = load_mcords_nc.load_mcords_nc(os.path.join(THIS_DIR, 'input_data', 'zeros_mcords.nc'))

    @unittest.skipIf(load_mcords_nc.NC, 'NETCDF on this version')
    def test_loadnc_failure(self):
        with self.assertRaises(ImportError):
            load_mcords_nc.load_mcords_nc(os.path.join(THIS_DIR, 'input_data', 'zeros_mcords.nc'))


if __name__ == '__main__':
    unittest.main()
