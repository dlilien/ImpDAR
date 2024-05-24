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
import numpy as np
from impdar.lib.RadarData._RadarDataSaving import CONVERSIONS_ENABLED
from impdar.lib.RadarData import RadarData
from matplotlib import colormaps

try:
    import matplotlib
    matplotlib.use('QT5Agg')
    from PyQt5 import QtWidgets, QtCore
    from impdar.gui.pickgui import InteractivePicker, VBPInputDialog, CropInputDialog, warn, plt
    app = QtWidgets.QApplication(sys.argv)
    qt = True
except ImportError:
    qt = False

if sys.version_info[0] >= 3:
    from unittest.mock import MagicMock, patch
else:
    from mock import MagicMock, patch


THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class DummyEvent:

    def accept(self):
        pass

    def ignore(self):
        pass


@unittest.skipIf(not qt, 'No Qt')
class TestInteractivePicker(unittest.TestCase):

    def setUp(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.ip = InteractivePicker(data)

    def test_other_lims(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        ip = InteractivePicker(data, xdat='dist')
        self.assertEqual(ip.x, 'dist')
        with self.assertRaises(ValueError):
            ip = InteractivePicker(data, xdat='dum')
        ip = InteractivePicker(data, ydat='twtt')
        self.assertEqual(ip.y, 'twtt')
        ip = InteractivePicker(data, ydat='depth')
        self.assertEqual(ip.y, 'depth')
        data.nmo_depth = data.travel_time
        ip = InteractivePicker(data, ydat='depth')
        self.assertEqual(ip.y, 'depth')
        with self.assertRaises(ValueError):
            ip = InteractivePicker(data, ydat='dum')

        with self.assertRaises(ValueError):
            ip = InteractivePicker(data, ydat='elev')
        data.elevation = np.arange(ip.dat.tnum)
        data.flags.elev = True
        ip = InteractivePicker(data, x_range=None)
        self.assertEqual(ip.x_range, (0, ip.dat.tnum))

    def test_PickNum(self):
        self.ip.pickNumberBox.setValue(1)

    def test_update_polarity(self):
        self.assertEqual(self.ip.dat.picks.pickparams.pol, 1)
        self.ip.wbw_radio.setChecked(True)
        self.assertEqual(self.ip.dat.picks.pickparams.pol, -1)

    def test_reverse_color(self):
        if hasattr(plt.cm, "get_cmap"):
            self.assertEqual(self.ip.im.get_cmap(), plt.cm.get_cmap(self.ip.ColorSelector.currentText()))
        else:
            self.assertEqual(self.ip.im.get_cmap(), colormaps[self.ip.ColorSelector.currentText()])
        self.ip._update_color_reversal(QtCore.Qt.Checked)
        if hasattr(plt.cm, "get_cmap"):
            self.assertEqual(self.ip.im.get_cmap(), plt.cm.get_cmap(self.ip.ColorSelector.currentText() + '_r'))
        else:
            self.assertEqual(self.ip.im.get_cmap(), colormaps[self.ip.ColorSelector.currentText() + '_r'])

    def test_select_lines_click(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        self.ip = InteractivePicker(data)
        event = DummyEvent()
        event.artist = self.ip.cline[0]
        self.ip._select_lines_click(event)
        self.assertEqual(self.ip.pickNumberBox.value(), 1)

        event.artist = 'dumdum'
        self.ip._select_lines_click(event)
        self.assertEqual(self.ip.pickNumberBox.value(), 1)

        event.artist = self.ip.cline[1]
        self.ip._select_lines_click(event)
        self.assertEqual(self.ip.pickNumberBox.value(), 5)

    def test_freq_update(self):
        p = self.ip.dat.picks.pickparams.plength
        self.ip._freq_update(678)
        self.assertEqual(self.ip.dat.picks.pickparams.freq, 678)
        self.assertFalse(p == self.ip.dat.picks.pickparams.plength)

    def test_add_pick(self):
        #blank pick
        self.ip._add_pick()
        self.assertTrue(self.ip.dat.picks.samp1.shape[0] == 1)
        
        # Overwrite slot
        self.ip._add_pick()
        self.assertTrue(self.ip.dat.picks.samp1.shape[0] == 1)

        # add to existing
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        self.ip = InteractivePicker(data)
        self.assertTrue(self.ip.dat.picks.samp1.shape[0] == 2)
        self.ip._add_pick()
        self.assertTrue(self.ip.dat.picks.samp1.shape[0] == 3)

        # Check that we can add a non-blank pick
        self.ip._add_pick(snum=10, tnum=2)
        self.assertTrue(self.ip.dat.picks.samp1.shape[0] == 3)
        self.assertTrue(self.ip.current_pick[1, 2] > 5)
        
        # snum but no tnum
        self.ip._add_pick(snum=10, tnum=None)
        self.assertTrue(self.ip.current_pick[1, 0] > 5)

    def test_color_select(self):
        self.ip._color_select('bone')
        self.assertTrue(self.ip.im.get_cmap(), 'bone')
        self.ip._color_select('CEGSIC')
        self.assertTrue(self.ip.im.get_cmap(), 'CEGSIC')

    def test_lim_update(self):
        self.ip._update_lims(-100, 100)
        self.assertEqual(self.ip.im.get_clim(), (-100, 100))
        with self.assertRaises(ValueError):
            self.ip._update_lims(100, -100)
        self.ip.minSpinner.setValue(-999)
        self.ip.maxSpinner.setValue(999)
        self.assertEqual(self.ip.im.get_clim(), (-999, 999))
        self.ip.minSpinner.setValue(1000)
        self.assertEqual(self.ip.im.get_clim(), (1000, 1001))

    def test_mode_update(self):
        self.ip._mode_update()
        self.ip._mode_update()
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data_picks.mat'))
        ip = InteractivePicker(data)
        ip._mode_update()
        ip._mode_update()

    def test_edit_lines_click_existingline(self):
        # First, plain left click
        # event has x and y data
        event = DummyEvent()
        event.xdata = 10.
        event.ydata = 1.0e-1
        event.button = 1

        self.ip._freq_update(800)

        # assume we have a pick
        self.ip._add_pick(snum=10, tnum=1)
        self.ip._add_point_pick = MagicMock()
        self.ip.update_lines = MagicMock()
        self.ip._edit_lines_click(event)
        self.assertTrue(self.ip._add_point_pick.called)
        self.assertTrue(self.ip.update_lines.called)

        # now nanpick
        self.ip._n_pressed = True
        self.ip._add_nanpick = MagicMock()
        self.ip.update_lines = MagicMock()
        # prevent a dialog box
        with patch('impdar.gui.pickgui.warn'):
            self.ip._edit_lines_click(event)
        self.assertTrue(self.ip._add_nanpick.called)
        self.assertTrue(self.ip.update_lines.called)
        self.ip._n_pressed = False

        # now delete pick
        event.button = 3
        self.ip._delete_picks = MagicMock()
        self.ip.update_lines = MagicMock()
        self.ip._edit_lines_click(event)
        self.assertTrue(self.ip._delete_picks.called)
        self.assertTrue(self.ip.update_lines.called)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    def test_edit_lines_click_newline(self):
        # First, plain left click
        # event has x and y data
        event = DummyEvent()
        event.xdata = 10.
        event.ydata = 1.0e-1
        event.button = 1

        # assume we have no picks
        self.ip._add_pick = MagicMock()
        self.ip.update_lines = MagicMock()
        self.ip._edit_lines_click(event)
        self.assertTrue(self.ip._add_pick.called)
        self.assertTrue(self.ip.update_lines.called)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    def test_add_point_pick(self):
        # need to mock a lot to not deal with actually doing any picking
        with patch('impdar.lib.picklib.packet_pick', return_value=np.zeros((5, ))) as mock1:
            with patch('impdar.lib.picklib.pick', return_value=np.zeros((5, self.ip.dat.tnum - 1))) as mock2:
                self.ip._add_pick(0, 0)
                self.ip._add_point_pick(0, self.ip.dat.tnum - 1)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    def test_add_nan_pick(self):
        with patch('impdar.lib.picklib.packet_pick', return_value=np.zeros((5, ))) as mock1:
            self.ip._add_pick(0, 0)
            self.ip._add_nanpick(1, 10)
        self.assertEqual(self.ip.dat.picks.lasttrace.snum[0], 1)
        self.assertEqual(self.ip.dat.picks.lasttrace.tnum[0], 10)


@unittest.skipIf(not qt, 'No Qt')
class TestInteractivePickerLoadingSaving(unittest.TestCase):

    def setUp(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.ip = InteractivePicker(data)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @patch('impdar.gui.pickgui.QMessageBox')
    def test_save_cancel_closeSAVE(self, patchsave):
        patchsave.return_value.exec_.return_value = patchsave.Save
        patchsave.return_value.Save = patchsave.Save
        self.ip.fn = None
        self.ip._save_as = MagicMock(return_value=True)
        event = DummyEvent()
        event.accept = MagicMock()
        self.ip._save_cancel_close(event)
        self.assertTrue(self.ip._save_as.called)
        self.assertTrue(event.accept.called)

        event = DummyEvent()
        event.ignore = MagicMock()
        self.ip._save_as = MagicMock(return_value=False)
        self.ip._save_cancel_close(event)
        self.assertTrue(self.ip._save_as.called)
        self.assertTrue(event.ignore.called)

        self.ip.fn = 'dummy'
        self.ip._save = MagicMock()
        event.accept = MagicMock()
        self.ip._save_cancel_close(event)
        self.assertTrue(self.ip._save.called)
        self.assertTrue(event.accept.called)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @patch('impdar.gui.pickgui.QMessageBox')
    def test_save_cancel_closeCANCEL(self, patchcancel):
        # patchcancel.exec_.return_value = patchcancel.Cancel
        patchcancel.return_value.exec_.return_value = patchcancel.Cancel
        patchcancel.return_value.Cancel = patchcancel.Cancel
        event = DummyEvent()
        event.ignore = MagicMock()
        self.ip._save_cancel_close(event)
        self.assertTrue(event.ignore.called)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @patch('impdar.gui.pickgui.QMessageBox')
    def test_save_cancel_closeClose(self, patchclose):
        # patchcancel.exec_.return_value = patchclose.Close
        patchclose.return_value.exec_.return_value = patchclose.Close
        patchclose.return_value.Close = patchclose.Close
        event = DummyEvent()
        event.accept = MagicMock()
        self.ip._save_cancel_close(event)
        self.assertTrue(event.accept.called)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @patch('impdar.gui.pickgui.QFileDialog')
    def test_load_cp(self, patchqfd):
        patchqfd.getOpenFileName.return_value = ('not_a_file', True)
        with self.assertRaises(IOError):
            self.ip._load_cp(DummyEvent())

        patchqfd.getOpenFileName.return_value = (os.path.join(THIS_DIR, 'input_data', 'small_data.mat'), True)
        with patch('impdar.gui.pickgui.warn') as patchwarn:
            self.ip._load_cp(DummyEvent())
            self.assertTrue(patchwarn.called)

        patchqfd.getOpenFileName.return_value = (os.path.join(THIS_DIR, 'input_data', 'cross_picked.mat'), True)
        self.ip._load_cp(DummyEvent())

        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.ip = InteractivePicker(data, ydat='depth')
        self.ip._load_cp(DummyEvent())

    def test_save(self):
        self.ip.fn = None
        with self.assertRaises(AttributeError):
            self.ip._save(DummyEvent())

        self.ip.fn = os.path.join(THIS_DIR, 'input_data', 'test_out.mat')
        self.ip._save(DummyEvent())
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_out.mat')))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @patch('impdar.gui.pickgui.QFileDialog')
    def test_save_as(self, patchqfd):
        patchqfd.getSaveFileName.return_value = (os.path.join(THIS_DIR, 'input_data', 'test_out.mat'), True)
        self.ip._save_as(DummyEvent())
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test_out.mat')))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @patch('impdar.gui.pickgui.QFileDialog')
    def test_export_csv(self, patchqfd):
        patchqfd.getSaveFileName.return_value = (os.path.join(THIS_DIR, 'input_data', 'test.csv'), True)
        self.ip._export_csv(DummyEvent())
        self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test.csv')))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test.csv'))

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @unittest.skipIf(not CONVERSIONS_ENABLED, 'No GDAL on this version')
    @patch('impdar.gui.pickgui.QFileDialog')
    def test_export_shp(self, patchqfd):
        patchqfd.getSaveFileName.return_value = (os.path.join(THIS_DIR, 'input_data', 'test.shp'), True)
        self.ip._export_shp(DummyEvent())
        for ext in ['shp', 'shx', 'prj', 'dbf']:
            self.assertTrue(os.path.exists(os.path.join(THIS_DIR, 'input_data', 'test.' + ext)))
            os.remove(os.path.join(THIS_DIR, 'input_data', 'test.' + ext))


@unittest.skipIf(not qt, 'No Qt')
class TestInteractivePickerProcessing(unittest.TestCase):

    def setUp(self):
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        self.ip = InteractivePicker(data)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    def test_ahfilt(self):
        self.ip.dat.adaptivehfilt = MagicMock()

        # takes a dummy event arg
        self.ip._ahfilt(DummyEvent())
        self.assertTrue(self.ip.dat.adaptivehfilt.called)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @patch('impdar.gui.pickgui.VBPInputDialog', exec_=lambda x: None, lims=(100, 200))
    def test_vbp(self, vbpmock):
        self.ip.dat.vertical_band_pass = MagicMock()
        self.ip._vbp(DummyEvent())
        self.assertTrue(self.ip.dat.vertical_band_pass.called)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    @patch('impdar.gui.pickgui.CropInputDialog', exec_=lambda x: None, top_or_bottom='top', inputtype='twtt')
    def test_crop(self, cropinputmock):
        self.ip.dat.crop = MagicMock()
        self.ip._crop(DummyEvent())
        self.assertTrue(self.ip.dat.crop.called)

    @unittest.skipIf(sys.version_info[0] < 3, 'Mock is only on 3+')
    def test_reverse(self):
        self.ip.dat.reverse = MagicMock()
        self.ip._reverse(DummyEvent())
        self.assertTrue(self.ip.dat.reverse.called)


@unittest.skipIf(not qt, 'No Qt')
class TestVBP(unittest.TestCase):

    def test_VBPInputDialog(self):
        vbp = VBPInputDialog()
        vbp._click_ok()
        self.assertTrue(vbp.lims == (50, 250))
        self.assertTrue(vbp.accepted)

        vbp = VBPInputDialog()
        vbp.minspin.setValue(2)
        vbp.maxspin.setValue(298)
        vbp._click_ok()
        self.assertTrue(vbp.lims == (2, 298))
        self.assertTrue(vbp.accepted)

        vbp = VBPInputDialog()
        vbp.minspin.setValue(299)
        vbp.maxspin.setValue(298)
        # click OK twice since we have bad lims
        vbp._click_ok()
        vbp._click_ok()
        self.assertTrue(vbp.lims == (297, 298))
        self.assertTrue(vbp.accepted)


@unittest.skipIf(not qt, 'No Qt')
class TestCrop(unittest.TestCase):

    def test_CropInputDialog(self):
        cid = CropInputDialog()
        cid._click_ok()
        self.assertTrue(cid.accepted)

        cid = CropInputDialog()
        cid.inputtype.setCurrentText('snum')
        self.assertTrue(cid.spinnerlabel.text() == 'Cutoff (sample num):')

        cid.inputtype.setCurrentText('twtt')
        self.assertTrue(cid.spinnerlabel.text() == 'Cutoff in TWTT (usec):')

        cid.inputtype.setCurrentText('depth')
        self.assertTrue(cid.spinnerlabel.text() == 'Cutoff in depth (m):')

        cid._click_ok()
        self.assertTrue(cid.accepted)


if __name__ == '__main__':
    unittest.main()
