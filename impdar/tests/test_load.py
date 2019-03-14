#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
import os
import unittest
import numpy as np
from impdar.lib import load

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestLoad(unittest.TestCase):

    def test_load_mat(self):
        data = load.load_mat(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.assertEqual(data.data.shape, (20, 40))

    def test_loadmat(self):
        data = load.load('mat', os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.assertEqual(data[0].data.shape, (20, 40))

    def test_loadgssi(self):
        data = load.load('gssi', os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT'))

    def test_loadbad(self):
        with self.assertRaises(ValueError):
            data = load.load('bad', os.path.join(THIS_DIR, 'input_data', 'small_data.bad'))

    def test_load_and_exitmat(self):
        data = load.load_and_exit('mat', os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))


if __name__ == '__main__':
    unittest.main()
