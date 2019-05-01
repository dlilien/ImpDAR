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
from impdar.lib import load_gssi

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Dummy(unittest.TestCase):

    def test_load_withDZG(self):
        load_gssi.DZT(os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT'))

    def test_load_withoutDZG(self):
        load_gssi.DZT(os.path.join(THIS_DIR, 'input_data', 'test_gssi_justdzt.DZT'))

    def test_save_withDZG(self):
        load_gssi.DZT(os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT')).save(os.path.join(THIS_DIR, 'input_data', 'test_gssi_raw.mat'))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_gssi_raw.mat'))

    def test_save_withoutDZG(self):
        load_gssi.DZT(os.path.join(THIS_DIR, 'input_data', 'test_gssi_justdzt.DZT')).save(os.path.join(THIS_DIR, 'input_data', 'test_gssi_justdzt_raw.mat'))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_gssi_justdzt_raw.mat'))


if __name__ == '__main__':
    unittest.main()
