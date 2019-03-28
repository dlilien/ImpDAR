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

import unittest
import os
from impdar.lib import load_gecko

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestLoadGecko(unittest.TestCase):
    def test_load_gecko(self):
        load_gecko.load_gecko(os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd'))
if __name__ == '__main__':
    unittest.main()
