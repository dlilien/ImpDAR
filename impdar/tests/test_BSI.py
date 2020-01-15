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

    def test_load_pe(self):
        load_bsi.load_bsi(os.path.join(THIS_DIR, 'input_data', 'test_bsi.h5'))


if __name__ == '__main__':
    unittest.main()
