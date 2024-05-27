#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the machinery of imppick.
"""
import sys
import unittest
try:
    from impdar.bin import imppick
    QT = True
except ImportError:
    QT = False

if sys.version_info[0] >= 3:
    from unittest.mock import patch, MagicMock
else:
    from mock import patch, MagicMock


class TestMain(unittest.TestCase):

    # mock so that we have no real gui
    @unittest.skipIf(not QT, 'No Qt')
    @patch('impdar.bin.imppick.QtWidgets.QApplication')
    @patch('impdar.bin.imppick.pickgui.InteractivePicker')
    def test_badinput(self, pick_patch, qapppatch):
        imppick.sys.argv = ['dummy']
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                imppick.main()

        imppick.sys.argv = ['dummy', '-xd']
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                imppick.main()

        imppick.sys.argv = ['dummy', '-yd']
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                imppick.main()

        imppick.sys.argv = ['dummy', '-xd', '-yd']
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                imppick.main()

        imppick.sys.argv = ['dummy', 'fn', 'fn2']
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(BaseException):
                imppick.main()

    @unittest.skipIf(not QT, 'No Qt')
    @patch('impdar.bin.imppick.QtWidgets.QApplication')
    @patch('impdar.bin.imppick.pickgui.InteractivePicker')
    @patch('impdar.bin.imppick.load.load')
    def test_pick_tnumsnum(self, load_patch, pick_patch, qapppatch):
        load_patch.return_value = [MagicMock()]
        imppick.sys.argv = ['dummy', 'fn']
        # this is supposed to exit when finished
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                imppick.main()

        self.assertTrue(load_patch.called)
        load_patch.asseert_called_with('mat', ['fn'])
        self.assertTrue(pick_patch.called)
        pick_patch.assert_called_with(load_patch.return_value[0], xdat='tnum', ydat='twtt')

    @unittest.skipIf(not QT, 'No Qt')
    @patch('impdar.bin.imppick.QtWidgets.QApplication')
    @patch('impdar.bin.imppick.pickgui.InteractivePicker')
    @patch('impdar.bin.imppick.load.load')
    def test_pick_tnumdepth(self, load_patch, pick_patch, qapppatch):
        load_patch.return_value = [MagicMock()]
        imppick.sys.argv = ['dummy', 'fn', '-yd']
        # this is supposed to exit when finished
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                imppick.main()

        self.assertTrue(load_patch.called)
        load_patch.asseert_called_with('mat', ['fn'])
        self.assertTrue(pick_patch.called)
        pick_patch.assert_called_with(load_patch.return_value[0], xdat='tnum', ydat='depth')

    @unittest.skipIf(not QT, 'No Qt')
    @patch('impdar.bin.imppick.QtWidgets.QApplication')
    @patch('impdar.bin.imppick.pickgui.InteractivePicker')
    @patch('impdar.bin.imppick.load.load')
    def test_pick_distsnum(self, load_patch, pick_patch, qapppatch):
        load_patch.return_value = [MagicMock()]
        imppick.sys.argv = ['dummy', 'fn', '-xd']
        # this is supposed to exit when finished

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                imppick.main()

        self.assertTrue(load_patch.called)
        load_patch.asseert_called_with('mat', ['fn'])
        self.assertTrue(pick_patch.called)
        pick_patch.assert_called_with(load_patch.return_value[0], xdat='dist', ydat='twtt')


if __name__ == '__main__':
    unittest.main()
