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

    @patch('impdar.bin.apdar.load')
    def test_load(self, load_patch):
        apdar.sys.argv = ['dummy', 'load', 'mat', 'fn.mat']
        apdar.main()
        self.assertTrue(load_patch.called)
        aca, kwca = load_patch.call_args
        self.assertEqual(kwca['fns_in'], ['fn.mat'])
        self.assertEqual(kwca['filetype'], 'mat')

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                apdar.sys.argv = ['dummy', 'load', 'notanintype', 'fn.mat']
                apdar.main()

    @patch('impdar.bin.apdar.process')
    def test_process(self, process_patch):
        apdar.sys.argv = ['dummy', 'proc', '-rev', 'fn.mat']
        apdar.main()
        self.assertTrue(process_patch.called)
        aca, kwca = process_patch.call_args
        self.assertEqual(kwca['fn'], ['fn.mat'])
        self.assertEqual(kwca['rev'], True)


if __name__ == '__main__':
    unittest.main()
