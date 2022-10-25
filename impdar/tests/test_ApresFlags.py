#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 bhills <benjaminhhills@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""

import os
import unittest
import numpy as np
from impdar.lib.ApresData import ApresFlags, TimeDiffFlags, QuadPolFlags


class TestFlags(unittest.TestCase):

    def setUp(self):
        self.apf = ApresFlags()
        self.tdf = TimeDiffFlags()
        self.qpf = QuadPolFlags()

    def test_BoolOutputConversion(self):
        # Make sure the value is as expected
        out = self.tdf.to_matlab()
        self.assertFalse(out['unwrap'])

        self.tdf.bed_pick = True
        out = self.tdf.to_matlab()
        self.assertTrue(out['bed_pick'])

        for attr in self.tdf.attrs:
            self.assertTrue(attr in out)

    def test_InputConversion(self):

        in_flags_bad_format = {'file_read_code': None,
                               'range': 0,
                               'stack': 1,
                               'uncertainty': 0}
        in_flags_bad = {'unwrap': False}

        with self.assertRaises(KeyError):
            self.apf.from_matlab(in_flags_bad)
        with self.assertRaises(TypeError):
            self.apf.from_matlab(in_flags_bad_format)

    def tearDown(self):
        del self.apf

if __name__ == '__main__':
    unittest.main()
