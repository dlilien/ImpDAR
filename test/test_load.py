#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import sys
import os
import unittest
from impdar.lib import load

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestLoad(unittest.TestCase):

    def test_loadmat(self):
        data = load.load('mat', os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.assertEqual(data[0].data.shape, (20, 40))

    def test_loadgssi(self):
        data = load.load('gssi', os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT'))

    def test_loadUoA(self):
        data = load.load('UoA_h5', os.path.join(THIS_DIR, 'input_data', 'UoA_dummy.h5'))

    def test_loadpe(self):
        data = load.load('pe', os.path.join(THIS_DIR, 'input_data', 'test_pe.DT1'))

    @unittest.skipIf(sys.version_info[0] < 3, 'Bytes are weird in 2')
    def test_loadgecko(self):
        data = load.load('gecko', os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd'))
        data = load.load('gecko', [os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd'), os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd')])

    def test_loadbad(self):
        with self.assertRaises(ValueError):
            data = load.load('bad', os.path.join(THIS_DIR, 'input_data', 'small_data.bad'))

    def test_loadtek(self):
        data = load.load('tek', os.path.join(THIS_DIR, 'input_data', 'test_tek.DAT'))

    def test_load_apresprofile(self):
        data = load.load('apres', os.path.join(THIS_DIR, 'input_data', 'apres_1.DAT'))

    def test_load_and_exitmat(self):
        data = load.load_and_exit('mat', os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), o=os.path.join(THIS_DIR, 'input_data', 'small_data_rawrrr.mat'))
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data_rawrrr.mat')))

    @unittest.skipIf(sys.version_info[0] < 3, 'Bytes are weird in 2')
    def test_load_and_exitgecko(self):
        load.load_and_exit('gecko', os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd'))
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_gecko_raw.mat')))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_gecko_raw.mat'))
        load.load_and_exit('gecko', [os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd'), os.path.join(THIS_DIR, 'input_data', 'test_gecko.gtd')])
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_gecko_raw.mat')))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_gecko_raw.mat'))

    def test_load_and_exitcustomfn(self):
        data = load.load_and_exit('mat', os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data_raw.mat')))

    @unittest.skipIf(sys.version_info[0] < 3, 'FileNotFoundError not in 2')
    def test_load_and_exiterror(self):
        # We dont have an output folder
        with self.assertRaises(FileNotFoundError):
            load.load_and_exit('mat', [os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], o='dummy')

    def tearDown(self):
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data_raw.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'small_data_raw.mat'))
        if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'small_data_rawrrr.mat')):
            os.remove(os.path.join(THIS_DIR, 'input_data', 'small_data_rawrrr.mat'))


if __name__ == '__main__':
    unittest.main()
