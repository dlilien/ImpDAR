#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Benjamin Hills <bhills@uw.edu>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the apdar command line prompt.
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

class TestApres(unittest.TestCase):

    @patch('impdar.bin.apdar.load_apres.load_apres')
    def test_load(self, load_patch):
        fn = os.path.join(THIS_DIR, 'input_data', 'apres_2.DAT')
        apdar.sys.argv = ['dummy', 'load', fn]
        apdar.main()
        self.assertTrue(load_patch.called)
        aca, kwca = load_patch.call_args
        self.assertEqual(aca[0], [fn])

    @patch('impdar.bin.apdar.full_processing')
    def test_proc(self, proc_patch):
        fn = os.path.join(THIS_DIR, 'input_data', 'apres_1.mat')
        apdar.sys.argv = ['dummy', 'proc', fn]
        apdar.main()
        fn = os.path.join(THIS_DIR, 'input_data', 'apres_2.DAT')
        apdar.sys.argv = ['dummy', 'load', fn]
        apdar.main()
        fn = os.path.join(THIS_DIR, 'input_data', 'apres_2_raw.mat')
        apdar.sys.argv = ['dummy', 'proc', fn]
        apdar.main()
        self.assertTrue(proc_patch.called)

    @patch('impdar.bin.apdar.ApresDiff')
    def test_load_diff(self, load_patch):
        fn1 = os.path.join(THIS_DIR, 'input_data', 'apres_1_proc.mat')
        fn2 = os.path.join(THIS_DIR, 'input_data', 'apres_2_proc.mat')
        apdar.sys.argv = ['dummy', 'diffload', fn1, fn2]
        apdar.main()
        self.assertTrue(load_patch.called)

    @patch('impdar.bin.apdar.full_differencing')
    def test_proc_diff(self, proc_patch):
        fn = os.path.join(THIS_DIR, 'input_data', 'diffdat.mat')
        apdar.sys.argv = ['dummy', 'diffproc', fn]
        apdar.main()
        self.assertTrue(proc_patch.called)

if __name__ == '__main__':
    unittest.main()