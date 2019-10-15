#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make sure that we can successfully read ramac/mala input files
"""

import os
import unittest
import numpy as np
from impdar.lib.load import load_ramac

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestRAMAC(unittest.TestCase):

    def test_load_ramac_withgps(self):
        a = load_ramac.load_ramac(os.path.join(THIS_DIR, 'input_data', 'ten_col.rd3'))
        b = load_ramac.load_ramac(os.path.join(THIS_DIR, 'input_data', 'ten_col.rad'))
        c = load_ramac.load_ramac(os.path.join(THIS_DIR, 'input_data', 'ten_col'))
        self.assertTrue(np.all(a.data == b.data))
        self.assertTrue(np.all(a.data == c.data))

    def test_load_ramac_nogps(self):
        load_ramac.load_ramac(os.path.join(THIS_DIR, 'input_data', 'ten_col_nogps.rd3'))


if __name__ == '__main__':
    unittest.main()
