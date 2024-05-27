#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""

import os
import unittest
import numpy as np
from impdar.lib.ApresData import ApresData
from impdar.lib.ApresData.load_apres import load_apres_single_file

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestRadarDataSaving(unittest.TestCase):

    def test_ReadWrite(self):
        dat = load_apres_single_file(os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT'))
        dat.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        # reload the object and check to make sure they are the same
        dat_reload = ApresData(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    def test_ReadWrite_h5(self):
        dat = load_apres_single_file(os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT'))
        dat.save(os.path.join(THIS_DIR, 'input_data', 'test_out.h5'))
        # reload the object and check to make sure they are the same
        dat_reload = ApresData(os.path.join(THIS_DIR, 'input_data', 'test_out.h5'))


    def tearDown(self):
        for fn in ['test_out.mat', 'test.shp', 'test.shx', 'test.prj', 'test.dbf']:
            if os.path.exists(os.path.join(THIS_DIR, 'input_data', fn)):
                os.remove(os.path.join(THIS_DIR, 'input_data', fn))

if __name__ == '__main__':
    unittest.main()
