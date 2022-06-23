#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the machinery of process. This is broken up to match where it would likely fail; tests process wrappers of various methods are with the tests of those methods
"""
import sys
import os
import unittest
from impdar.bin import apdar

if sys.version_info[0] >= 3:
    from unittest.mock import patch, MagicMock
else:
    from mock import patch, MagicMock

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestMain(unittest.TestCase):

    @patch('impdar.bin.apdar.load_apres.load_apres')
    def test_load(self, load_patch):
        apdar.sys.argv = ['dummy', 'load', './input_data/apres_1.DAT']
        apdar.main()
        self.assertTrue(load_patch.called)
        aca, kwca = load_patch.call_args
        self.assertEqual(aca[0], ['./input_data/apres_1.DAT'])

if __name__ == '__main__':
    unittest.main()
