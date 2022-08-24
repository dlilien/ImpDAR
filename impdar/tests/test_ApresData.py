#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 bhills <benjaminhhills@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the basics of ApresData
"""

import sys
import os
import unittest
from impdar.lib.ApresData import ApresData

if sys.version_info[0] >= 3:
    from unittest.mock import patch
else:
    from mock import patch

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestApresDataLoading(unittest.TestCase):

    def test_ReadSucceeds(self):
        data = ApresData(os.path.join(THIS_DIR, 'input_data', 'apres_1.mat'))
        self.assertEqual(data.data.shape, (data.bnum, data.cnum, data.snum))

    def test_badread(self):
        # Data but not other attrs
        with self.assertRaises(KeyError):
            data = ApresData(os.path.join(THIS_DIR, 'input_data', 'nonimpdar_matlab.mat'))

        # All other attrs, no data
        with self.assertRaises(KeyError):
            data = ApresData(os.path.join(THIS_DIR, 'input_data', 'nonimpdar_justmissingdat.mat'))

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_out.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

if __name__ == '__main__':
    unittest.main()
