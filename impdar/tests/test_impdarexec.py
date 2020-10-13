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
from impdar.bin import impdarexec

if sys.version_info[0] >= 3:
    from unittest.mock import patch, MagicMock
else:
    from mock import patch, MagicMock

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestMain(unittest.TestCase):

    @patch('impdar.bin.impdarexec.load.load_and_exit')
    def test_load(self, load_patch):
        impdarexec.sys.argv = ['dummy', 'load', 'mat', 'fn.mat']
        impdarexec.main()
        self.assertTrue(load_patch.called)
        aca, kwca = load_patch.call_args
        self.assertEqual(kwca['fns_in'], ['fn.mat'])
        self.assertEqual(kwca['filetype'], 'mat')

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impdarexec.sys.argv = ['dummy', 'load', 'notanintype', 'fn.mat']
                impdarexec.main()

    @patch('impdar.bin.impdarexec.process.process_and_exit')
    def test_process(self, process_patch):
        impdarexec.sys.argv = ['dummy', 'proc', '-rev', 'fn.mat']
        impdarexec.main()
        self.assertTrue(process_patch.called)
        aca, kwca = process_patch.call_args
        self.assertEqual(kwca['fn'], ['fn.mat'])
        self.assertEqual(kwca['rev'], True)

    @patch('impdar.bin.impdarexec.plot.plot')
    def test_plot(self, plot_patch):
        impdarexec.sys.argv = ['dummy', 'plot', 'fn.mat']
        impdarexec.main()
        self.assertTrue(plot_patch.called)
        aca, kwca = plot_patch.call_args
        self.assertEqual(kwca['fns'], ['fn.mat'])

    @patch('impdar.bin.impdarexec.convert.convert')
    def test_convert(self, convert_patch):
        impdarexec.sys.argv = ['dummy', 'convert', 'fn.mat', 'shp']
        impdarexec.main()
        self.assertTrue(convert_patch.called)
        aca, kwca = convert_patch.call_args
        self.assertEqual(kwca['fns_in'], ['fn.mat'])
        self.assertEqual(kwca['out_fmt'], 'shp')

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impdarexec.sys.argv = ['dummy', 'convert', 'fn.mat', 'notanoutput']
                impdarexec.main()


if __name__ == '__main__':
    unittest.main()
