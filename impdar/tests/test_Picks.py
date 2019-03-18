#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Test loading of pick data
"""

import os
import unittest
import numpy as np
from impdar.lib import load

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestLoad(unittest.TestCase):

    def test_load_mat(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        self.assertEqual(data.data.shape, (20, 40))


if __name__ == '__main__':
    unittest.main()
