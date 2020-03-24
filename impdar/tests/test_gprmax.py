#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make sure that we can successfully read gprMax input files
"""

import os
import unittest
from impdar.lib.load import load_gprMax

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestGPRMax(unittest.TestCase):

    def test_load(self):
        load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', 'rectangle_gprMax_Bscan.h5'))

    def tearDown(self):
        fn = os.path.join(THIS_DIR, 'input_data', 'rectangle_gprMax_Bscan_raw.mat')
        if os.path.exists(fn):
            os.path.remove(fn)


if __name__ == '__main__':
    unittest.main()
