#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Test converting between filetypes
"""
import os
import pytest
import unittest
from impdar.lib import convert
from impdar.lib.RadarData._RadarDataSaving import CONVERSIONS_ENABLED
from impdar.lib.load.load_segy import SEGY

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestConvert(unittest.TestCase):

    @unittest.skipIf(not CONVERSIONS_ENABLED, 'No GDAL on this version')
    def test_guessload2shp(self):
        convert.convert(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), 'shp')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data.shp')))

        convert.convert([os.path.join(THIS_DIR, 'input_data', 'test_pe.DT1')], 'shp')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_pe.shp')))

        convert.convert([os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT')], 'shp')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_gssi.shp')))

    @unittest.skipIf(not CONVERSIONS_ENABLED, 'No GDAL on this version')
    def test_knownload2shp(self):
        convert.convert(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), 'shp', in_fmt='mat')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data.shp')))

        convert.convert([os.path.join(THIS_DIR, 'input_data', 'test_pe.DT1')], 'shp', in_fmt='pe')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_pe.shp')))

        convert.convert([os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT')], 'shp', in_fmt='gssi')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_gssi.shp')))

    @unittest.skipIf(not CONVERSIONS_ENABLED, 'No GDAL on this version')
    def test_knownload2mat(self):
        convert.convert([os.path.join(THIS_DIR, 'input_data', 'test_pe.DT1')], 'mat', in_fmt='pe')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_pe.mat')))

        convert.convert([os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT')], 'mat', in_fmt='gssi')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_gssi.mat')))

    @unittest.skipIf(SEGY, 'SEGY enabled, this is a failure test')
    def test_nosegy(self):
        with self.assertRaises(ImportError):
            convert.convert([os.path.join(THIS_DIR, 'input_data', 'test_pe.DT1')], 'mat', in_fmt='segy')

        with self.assertRaises(ImportError):
            convert.convert([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], 'sgy', in_fmt='mat')

    @unittest.skipIf(not SEGY, 'SEGY needed for this test')
    def test_segy_save(self):
        pytest.importorskip('segyio', reason='No SEGY on this version')
        convert.convert(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), 'sgy', in_fmt='mat')
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data.sgy')))

        convert.convert(os.path.join(THIS_DIR, 'input_data', 'shots0001_0200.segy'), 'mat', in_fmt='segy')

    def test_badinsout(self):
        with self.assertRaises(ValueError):
            convert.convert([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], 'dummy')
        with self.assertRaises(ValueError):
            convert.convert([os.path.join(THIS_DIR, 'input_data', 'small_data.wtf')], 'shp')

    def tearDown(self):
        for ext in ['shp', 'shx', 'dbf', 'prj', 'sgy']:
            for pref in ['small_data', 'test_gssi', 'test_pe']:
                if os.path.exists(os.path.join(THIS_DIR, 'input_data', pref + '.' + ext)):
                    os.remove(os.path.join(THIS_DIR, 'input_data', pref + '.' + ext))
        for pref in ['test_gssi', 'test_pe', 'shots0001_0200']:
            if os.path.exists(os.path.join(THIS_DIR, 'input_data', pref + '.mat')):
                os.remove(os.path.join(THIS_DIR, 'input_data', pref + '.mat'))


if __name__ == '__main__':
    unittest.main()
