#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import sys
sys.modules['segyio'] = None
sys.modules['_segyio'] = None

import os
import unittest
import numpy as np
from impdar.lib import load

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


@unittest.skipIf(sys.version_info[0] < 3, 'Excluding segyio fails')
class TestLoadNoSEGY(unittest.TestCase):

    def test_loadmat(self):
        # We want normal functionality if not SEGY
        data = load.load('mat', os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.assertEqual(data[0].data.shape, (20, 40))

    def test_segyimporterror(self):
        # Fail if use segy
        with self.assertRaises(ImportError):
            data = load.load('segy', os.path.join(THIS_DIR, 'input_data', 'small_data.segy'))


if __name__ == '__main__':
    unittest.main()
