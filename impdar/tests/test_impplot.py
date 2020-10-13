#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the machinery of impplot.
This is broken up to match where it would likely fail;
tests process wrappers of various methods are with the tests of those methods
"""
import sys
import unittest
from impdar.bin import impplot

if sys.version_info[0] >= 3:
    from unittest.mock import patch, MagicMock
else:
    from mock import patch, MagicMock


class TestMain(unittest.TestCase):

    # mock so that we have no real processing
    @patch('impdar.bin.impplot.plot.plot')
    def test_badinput(self, plot_patch):
        impplot.sys.argv = ['dummy']

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                impplot.main()

        impplot.sys.argv = ['dummy', 'dummy']
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                impplot.main()

        impplot.sys.argv = ['dummy', 'rg']
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                impplot.main()

        impplot.sys.argv = ['dummy', 'dummy', 'fn']

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impplot.main()

    @patch('impdar.bin.impplot.plot.plot')
    def test_rg(self, plot_patch):
        impplot.sys.argv = ['dummy', 'rg', 'fn']
        impplot.main()
        self.assertTrue(plot_patch.called)
        aca, kwca = plot_patch.call_args
        self.assertEqual(aca[0], ['fn'])

        # we can let these default, but if touched must be None
        if 'power' in kwca:
            self.assertIsNone(kwca['power'])
        if 'tr' in kwca:
            self.assertIsNone(kwca['tr'])

    @patch('impdar.bin.impplot.plot.plot')
    def test_power(self, plot_patch):
        impplot.sys.argv = ['dummy', 'power', 'fn', '16']
        impplot.main()
        self.assertTrue(plot_patch.called)
        aca, kwca = plot_patch.call_args
        self.assertEqual(aca[0], ['fn'])
        self.assertEqual(kwca['power'], 16)
        if 'tr' in kwca:
            self.assertIsNone(kwca['tr'])

        impplot.sys.argv = ['dummy', 'power', 'fn', '16.']

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impplot.main()

    @patch('impdar.bin.impplot.plot.plot')
    def test_traces(self, plot_patch):
        impplot.sys.argv = ['dummy', 'traces', 'fn', '8', '16']
        impplot.main()
        self.assertTrue(plot_patch.called)
        aca, kwca = plot_patch.call_args
        self.assertEqual(aca[0], ['fn'])
        self.assertEqual(kwca['tr'], (8, 16))
        if 'power' in kwca:
            self.assertIsNone(kwca['power'])

        impplot.sys.argv = ['dummy', 'traces', 'fn', '8', '16.']

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impplot.main()

        impplot.sys.argv = ['dummy', 'traces', 'fn', '8.', '16']

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impplot.main()

        impplot.sys.argv = ['dummy', 'traces', 'fn', '8']

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                impplot.main()


if __name__ == '__main__':
    unittest.main()
