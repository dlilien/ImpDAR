#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Test the gecko file import for data from the St. Olaf HF radar.

Author:
Benjamin Hills
benjaminhhills@gmail.com
University of Washington
Earth and Space Sciences

Mar 28 2019

"""
import sys
import os
import unittest
from impdar.lib.load import load_olaf

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestLoadGecko(unittest.TestCase):

    @unittest.skipIf(sys.version_info[0] < 3, 'Bytes are weird in 2')
    def test_load_gecko(self):
        load_olaf.load_olaf(os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd'), channel=1)
        load_olaf.load_olaf(os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd'), channel=2)

if __name__ == '__main__':
    unittest.main()
