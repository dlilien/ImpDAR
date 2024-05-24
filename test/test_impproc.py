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
from impdar.bin import impproc
from impdar.lib import NoInitRadarData

if sys.version_info[0] >= 3:
    from unittest.mock import patch, MagicMock
else:
    from mock import patch, MagicMock

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestMain(unittest.TestCase):

    # mock so that we have no real processing
    @patch('impdar.bin.impproc.agc')
    @patch('impdar.bin.impproc.load')
    def test_inputfile(self, load_patch, agc_patch):
        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'agc', os.path.join(THIS_DIR, 'input_data', 'small_data.mat')]
        impproc.main()
        self.assertTrue(agc_patch.called)
        self.assertTrue(load_patch.called)
        load_patch.assert_called_with('mat', [os.path.join(THIS_DIR, 'input_data', 'small_data.mat')])

        impproc.sys.argv = ['dummy', 'agc', os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')]
        impproc.main()
        self.assertTrue(agc_patch.called)
        self.assertTrue(load_patch.called)
        load_patch.assert_called_with('mat', [os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')])

    # mock so that we have no real processing
    @patch('impdar.bin.impproc.agc')
    @patch('impdar.bin.impproc.load')
    def test_outputfile(self, load_patch, agc_patch):
        load_patch.return_value = [MagicMock()]
        rd_patch = load_patch.return_value
        rd_patch[0].save = MagicMock()
        impproc.sys.argv = ['dummy', 'agc', '-o', 'dummy', os.path.join(THIS_DIR, 'input_data', 'small_data.mat')]
        impproc.main()
        self.assertTrue(agc_patch.called)
        self.assertTrue(load_patch.called)
        load_patch.assert_called_with('mat', [os.path.join(THIS_DIR, 'input_data', 'small_data.mat')])
        self.assertTrue(rd_patch[0].save.called)
        rd_patch[0].save.assert_called_with('dummy')

    @patch('impdar.bin.impproc.agc')
    @patch('impdar.bin.impproc.load')
    def test_outputraw(self, load_patch, agc_patch):
        load_patch.return_value = [MagicMock()]
        rd_patch = load_patch.return_value
        for p in rd_patch:
            p.save = MagicMock()
        impproc.sys.argv = ['dummy', 'agc', os.path.join(THIS_DIR, 'input_data', 'small_data_raw.mat')]
        impproc.main()
        self.assertTrue(load_patch.called)
        load_patch.assert_called_with('mat', [os.path.join(THIS_DIR, 'input_data', 'small_data_raw.mat')])
        for p in rd_patch:
            self.assertTrue(p.save.called)
            p.save.assert_called_with(os.path.join(THIS_DIR, 'input_data', 'small_data_agc.mat'))

    @patch('impdar.bin.impproc.agc')
    @patch('impdar.bin.impproc.load')
    def test_outputmultiple(self, load_patch, agc_patch):
        load_patch.return_value = [MagicMock(), MagicMock()]
        rd_patch = load_patch.return_value
        for p in rd_patch:
            p.save = MagicMock()
        impproc.sys.argv = ['dummy', 'agc', '-o', 'dummy', os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data.mat')]
        impproc.main()
        for p in rd_patch:
            self.assertTrue(p.save.called)
            p.save.assert_called_with(os.path.join('dummy', 'small_data_agc.mat'))

    @patch('impdar.bin.impproc.agc')
    @patch('impdar.bin.impproc.load')
    def test_outputmultipleraw(self, load_patch, agc_patch):
        load_patch.return_value = [MagicMock(), MagicMock()]
        rd_patch = load_patch.return_value
        for p in rd_patch:
            p.save = MagicMock()
        impproc.sys.argv = ['dummy', 'agc', '-o', 'dummy', os.path.join(THIS_DIR, 'input_data', 'small_data_raw.mat'), os.path.join(THIS_DIR, 'input_data', 'small_data_raw.mat')]
        impproc.main()
        for p in rd_patch:
            self.assertTrue(p.save.called)
            p.save.assert_called_with(os.path.join('dummy', 'small_data_agc.mat'))

    def test_help(self):
        with self.assertRaises(BaseException):
            impproc.sys.argv = ['dummy']
            impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'dummy']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'ahfilt']
                impproc.main()


class TestInputs(unittest.TestCase):

    @patch('impdar.bin.impproc.agc')
    @patch('impdar.bin.impproc.load')
    def test_agc(self, load_patch, agc_patch):
        load_patch.return_value = [MagicMock()]
        window = 10
        impproc.sys.argv = ['dummy', 'agc', 'dummy.mat', '-window', str(window)]
        impproc.main()
        self.assertTrue(agc_patch.called)
        aca, kwca = agc_patch.call_args
        self.assertEqual(kwca['window'], window)

        window = 50
        impproc.sys.argv = ['dummy', 'agc', 'dummy.mat', '-window', str(window)]
        impproc.main()
        self.assertTrue(agc_patch.called)

        aca, kwca = agc_patch.call_args
        self.assertEqual(kwca['window'], window)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'agc', 'dummy.mat', '-window', '10.1']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'agc', 'dummy.mat', '-window', 'badint']
                impproc.main()

    @patch('impdar.bin.impproc.vbp')
    @patch('impdar.bin.impproc.load')
    def test_vbp(self, load_patch, vbp_patch):
        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'vbp', '10', '20', 'dummy.mat']
        impproc.main()
        self.assertTrue(vbp_patch.called)

        aca, kwca = vbp_patch.call_args
        self.assertEqual(kwca['low_MHz'], 10)
        self.assertEqual(kwca['high_MHz'], 20)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'vbp', 'dummy.mat', '10', '20']
                impproc.main()

    @patch('impdar.bin.impproc.rev')
    @patch('impdar.bin.impproc.load')
    def test_rev(self, load_patch, rev_patch):
        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'rev', 'dummy.mat']
        impproc.main()
        self.assertTrue(rev_patch.called)

    @patch('impdar.bin.impproc.ahfilt')
    @patch('impdar.bin.impproc.load')
    def test_ahfilt(self, load_patch, ahfilt_patch):
        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'ahfilt', '1000', 'dummy.mat']
        impproc.main()
        self.assertTrue(ahfilt_patch.called)

    @patch('impdar.bin.impproc.nmo')
    @patch('impdar.bin.impproc.load')
    def test_nmo(self, load_patch, nmo_patch):
        load_patch.return_value = [MagicMock()]
        sep = 123.4
        impproc.sys.argv = ['dummy', 'nmo', str(sep), 'dummy.mat']
        impproc.main()
        self.assertTrue(nmo_patch.called)
        aca, kwca = nmo_patch.call_args
        self.assertEqual(kwca['ant_sep'], sep)

        impproc.sys.argv = ['dummy', 'nmo', '--uice', str(10), str(sep), 'dummy.mat']
        impproc.main()
        aca, kwca = nmo_patch.call_args
        self.assertEqual(kwca['ant_sep'], sep)
        self.assertEqual(kwca['uice'], 10.)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'nmo', '--uice', 'badvel', str(sep), 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'nmo', '--uair', str(10), str(sep), 'dummy.mat']
        impproc.main()

        aca, kwca = nmo_patch.call_args
        self.assertEqual(kwca['ant_sep'], sep)
        self.assertEqual(kwca['uair'], 10.)
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'nmo', '--uair', 'badvel', str(sep), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'nmo', 'dummy.mat', str(sep)]
                impproc.main()

    @patch('impdar.bin.impproc.interp')
    @patch('impdar.bin.impproc.load')
    def test_interp(self, load_patch, interp_patch):
        load_patch.return_value = [MagicMock()]
        spacing = 10.
        impproc.sys.argv = ['dummy', 'interp', str(spacing), 'dummy.mat']
        impproc.main()
        self.assertTrue(interp_patch.called)
        aca, kwca = interp_patch.call_args
        self.assertEqual(kwca['spacing'], spacing)
        self.assertEqual(kwca['extrapolate'], False)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'interp', 'dummy.mat', os.path.join(THIS_DIR, 'input_data', 'small_data.mat')]
                impproc.main()

        impproc.sys.argv = ['dummy', 'interp', '--gps_fn', 'dummy', str(spacing), 'dummy.mat']
        impproc.main()
        aca, kwca = interp_patch.call_args
        self.assertEqual(kwca['spacing'], spacing)
        self.assertEqual(kwca['gps_fn'], 'dummy')

        impproc.sys.argv = ['dummy', 'interp', '--offset', str(10), str(spacing), 'dummy.mat']
        impproc.main()
        aca, kwca = interp_patch.call_args
        self.assertEqual(kwca['spacing'], spacing)
        self.assertEqual(kwca['offset'], 10.)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'interp', '--offset', 'badfloat', str(spacing), 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'interp', '--minmove', str(10), str(spacing), 'dummy.mat']
        impproc.main()
        aca, kwca = interp_patch.call_args
        self.assertEqual(kwca['spacing'], spacing)
        self.assertEqual(kwca['minmove'], 10.)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'interp', '--minmove', 'badfloat', str(spacing), 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'interp', '--extrapolate', str(spacing), 'dummy.mat']
        impproc.main()
        aca, kwca = interp_patch.call_args
        self.assertEqual(kwca['spacing'], spacing)
        self.assertEqual(kwca['extrapolate'], True)

    @patch('impdar.bin.impproc.concat')
    @patch('impdar.bin.impproc.load')
    def test_cat(self, load_patch, cat_patch):

        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'cat', 'dummy.mat', 'dummy.mat']
        impproc.main()
        self.assertTrue(cat_patch.called)

    @patch('impdar.bin.impproc.elev')
    @patch('impdar.bin.impproc.load')
    def test_elev(self, load_patch, elev_patch):
        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'elev', 'dummy.mat']
        impproc.main()
        self.assertTrue(elev_patch.called)

    @patch('impdar.bin.impproc.hfilt')
    @patch('impdar.bin.impproc.load')
    def test_hfilt(self, load_patch, hfilt_patch):
        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'hfilt', '10', '20', 'dummy.mat']
        impproc.main()
        self.assertTrue(hfilt_patch.called)
        aca, kwca = hfilt_patch.call_args
        self.assertEqual(kwca['start_trace'], 10)
        self.assertEqual(kwca['end_trace'], 20)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'hfilt', '10', 'dummy', 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'hfilt', 'dummy', '10', 'dummy.mat']
                impproc.main()

    @patch('impdar.bin.impproc.crop')
    @patch('impdar.bin.impproc.load')
    def test_crop(self, load_patch, crop_patch):
        load_patch.return_value = [MagicMock()]
        lim = 10
        for top_or_bottom in ['top', 'bottom']:
            for dimension in ['snum', 'twtt', 'depth', 'pretrig']:
                impproc.sys.argv = ['dummy', 'crop', top_or_bottom, dimension, str(lim), 'dummy.mat']
                impproc.main()
                self.assertTrue(crop_patch.called)
                aca, kwca = crop_patch.call_args
                self.assertEqual(kwca['top_or_bottom'], top_or_bottom)
                self.assertEqual(kwca['dimension'], dimension)
                self.assertEqual(kwca['lim'], lim)

        # Now bad entries
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'crop', 'top', 'bad', str(lim), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'crop', 'bad', 'snum', str(lim), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'crop', 'top', 'snum', 'notgood', 'dummy.mat']
                impproc.main()

    @patch('impdar.bin.impproc.hcrop')
    @patch('impdar.bin.impproc.load')
    def test_hcrop(self, load_patch, hcrop_patch):
        load_patch.return_value = [MagicMock()]
        lim = 10
        for left_or_right in ['left', 'right']:
            for dimension in ['tnum', 'dist']:
                impproc.sys.argv = ['dummy', 'hcrop', left_or_right, dimension, str(lim), 'dummy.mat']
                impproc.main()
                self.assertTrue(hcrop_patch.called)
                aca, kwca = hcrop_patch.call_args
                self.assertEqual(kwca['left_or_right'], left_or_right)
                self.assertEqual(kwca['dimension'], dimension)
                self.assertEqual(kwca['lim'], lim)

        # Now bad entries
        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'hcrop', 'left', 'bad', str(lim), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'hcrop', 'bad', 'tnum', str(lim), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'hcrop', 'left', 'tnum', 'notgood', 'dummy.mat']
                impproc.main()

    @patch('impdar.bin.impproc.restack')
    @patch('impdar.bin.impproc.load')
    def test_restack(self, load_patch, restack_patch):
        load_patch.return_value = [MagicMock()]
        interval = 3
        impproc.sys.argv = ['dummy', 'restack', str(interval), 'dummy.mat']
        impproc.main()
        self.assertTrue(restack_patch.called)
        aca, kwca = restack_patch.call_args
        self.assertEqual(kwca['traces'], interval)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'restack', 'bad', 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'restack', '0.1', 'dummy.mat']
                impproc.main()

    @patch('impdar.bin.impproc.rgain')
    @patch('impdar.bin.impproc.load')
    def test_rgain(self, load_patch, rgain_patch):
        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'rgain', 'dummy.mat']
        impproc.main()
        self.assertTrue(rgain_patch.called)

        slope = 10
        impproc.sys.argv = ['dummy', 'rgain', '-slope', str(slope), 'dummy.mat']
        impproc.main()
        aca, kwca = rgain_patch.call_args
        self.assertEqual(kwca['slope'], slope)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'rgain', '-slope', 'bad', 'dummy.mat']
                impproc.main()

    @patch('impdar.bin.impproc.mig')
    @patch('impdar.bin.impproc.load')
    def test_migrateTypes(self, load_patch, migrate_patch):
        load_patch.return_value = [MagicMock()]
        impproc.sys.argv = ['dummy', 'migrate', 'dummy.mat']
        impproc.main()
        self.assertTrue(migrate_patch.called)

        # mtype tests
        for mtype in ['stolt', 'kirch', 'phsh', 'tk', 'sustolt', 'sumigtk', 'sumigffd']:
            impproc.sys.argv = ['dummy', 'migrate', '--mtype', mtype, 'dummy.mat']
            impproc.main()
            aca, kwca = migrate_patch.call_args
            self.assertEqual(kwca['mtype'], mtype)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--mtype', 'bad', 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'migrate', '--nearfield', 'dummy.mat']
        impproc.main()
        aca, kwca = migrate_patch.call_args
        self.assertEqual(kwca['nearfield'], True)

        badint = 0.1
        goodint = 10
        worseint = 'hello'

        impproc.sys.argv = ['dummy', 'migrate', '--htaper', str(goodint), 'dummy.mat']
        impproc.main()
        aca, kwca = migrate_patch.call_args
        self.assertEqual(kwca['htaper'], goodint)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--htaper', str(badint), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--htaper', str(worseint), 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'migrate', '--vtaper', str(goodint), 'dummy.mat']
        impproc.main()
        aca, kwca = migrate_patch.call_args
        self.assertEqual(kwca['vtaper'], goodint)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--vtaper', str(badint), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--vtaper', str(worseint), 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'migrate', '--nxpad', str(goodint), 'dummy.mat']
        impproc.main()
        aca, kwca = migrate_patch.call_args
        self.assertEqual(kwca['nxpad'], goodint)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--nxpad', str(badint), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--nxpad', str(worseint), 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'migrate', '--tmig', str(goodint), 'dummy.mat']
        impproc.main()
        aca, kwca = migrate_patch.call_args
        self.assertEqual(kwca['tmig'], goodint)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--tmig', str(badint), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--tmig', str(worseint), 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'migrate', '--verbose', str(goodint), 'dummy.mat']
        impproc.main()
        aca, kwca = migrate_patch.call_args
        self.assertEqual(kwca['verbose'], goodint)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--verbose', str(badint), 'dummy.mat']
                impproc.main()

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--verbose', str(worseint), 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'migrate', '--vel', str(goodint), 'dummy.mat']
        impproc.main()
        aca, kwca = migrate_patch.call_args
        self.assertEqual(kwca['vel'], goodint)

        argparse_mock = MagicMock()
        with patch('argparse.ArgumentParser._print_message', argparse_mock):
            with self.assertRaises(SystemExit):
                impproc.sys.argv = ['dummy', 'migrate', '--vel', 'dummy', 'dummy.mat']
                impproc.main()

        impproc.sys.argv = ['dummy', 'migrate', '--vel_fn', str(goodint), 'dummy.mat']
        impproc.main()
        aca, kwca = migrate_patch.call_args
        self.assertEqual(kwca['vel_fn'], str(goodint))


class TestProc(unittest.TestCase):

    def setUp(self):
        self.data = NoInitRadarData.NoInitRadarDataFiltering()

    def test_hfilt(self):
        impproc.hfilt(self.data)

    def test_ahfilt(self):
        impproc.ahfilt(self.data)

    def test_rev(self):
        impproc.rev(self.data)

    def test_elev(self):
        self.data.nmo_depth = self.data.travel_time * 1.68e8 / 2.
        impproc.elev(self.data)

    def test_vbp(self):
        impproc.vbp(self.data, 0.1, 100.)

    def test_crop(self):
        impproc.crop(self.data, 2)

    def test_hcrop(self):
        impproc.crop(self.data, 2)

    def test_nmo(self):
        impproc.nmo(self.data)

    def test_restack(self):
        impproc.restack(self.data)

    def test_rgain(self):
        impproc.rgain(self.data)

    def test_agc(self):
        impproc.agc(self.data)

    def test_mig(self):
        impproc.mig(self.data)


if __name__ == '__main__':
    unittest.main()
