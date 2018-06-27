#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the basics of RadarData
"""

import unittest
from impdar.lib.RadarData import RadarData, RadarFlags


class TestFlags(unittest.TestCase):

    def setUp(self):
        self.rdf = RadarFlags()

    def test_BoolOutputConversion(self):
        out = self.rdf.to_matlab()
        self.assertFalse(out['reverse'])

    def tearDown(self):
        del self.rdf


if __name__ == '__main__':
    unittest.main()
