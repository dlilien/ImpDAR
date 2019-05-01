#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""

import os
import unittest
import numpy as np
from impdar.lib.RadarFlags import RadarFlags


class TestFlags(unittest.TestCase):

    def setUp(self):
        self.rdf = RadarFlags()

    def test_BoolOutputConversion(self):
        # Make sure the value is as expected
        self.rdf.reverse = False
        out = self.rdf.to_matlab()
        self.assertFalse(out['reverse'])

        self.rdf.rgain = True
        out = self.rdf.to_matlab()
        self.assertTrue(out['rgain'])

        for attr in self.rdf.attrs:
            self.assertTrue(attr in out)

    def test_InputConversion(self):
        
        in_flags_bad_format = {'agc': 0,
                               'batch': 0,
                               'bpass': np.array([0., 0., 0.]),
                               'crop': np.array([0., 0., 0.]),
                               'elev': 0,
                               'hfilt': np.array([0., 0.]),
                               'interp': np.array([0., 0.]),
                               'mig': 0,
                               'nmo': np.array([0., 0.]),
                               'restack': 0,
                               'reverse': 0,
                               'rgain': 0}
        # in_flags_random_arg = {'unknown_val': False}
        in_flags_bad = {'reverse': True}

        with self.assertRaises(KeyError):
            self.rdf.from_matlab(in_flags_bad)
        with self.assertRaises(TypeError):
            self.rdf.from_matlab(in_flags_bad_format)
        # self.rdf.from_matlab(in_flags_random_arg)

    def tearDown(self):
        del self.rdf


if __name__ == '__main__':
    unittest.main()
