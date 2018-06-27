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
import os
import unittest
import numpy as np
from impdar.lib.RadarData import RadarData, RadarFlags

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


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


class TestRadarData(unittest.TestCase):

    def test_ReadSucceeds(self):
        RadarData(os.path.join(THIS_DIR, 'input_data/small_data.mat'))


if __name__ == '__main__':
    unittest.main()
