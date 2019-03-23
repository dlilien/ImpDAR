#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
import sys
import os
import unittest
import numpy as np
try:
    import matplotlib
    matplotlib.use('QT5Agg')
    from PyQt5 import QtCore, QtWidgets, QtGui
    from impdar.gui.pickgui import InteractivePicker, VBPInputDialog, CropInputDialog
    from PyQt5.QtTest import QTest
    app = QtWidgets.QApplication(sys.argv)
    qt = True
except ImportError:
    qt = False
from impdar.lib.RadarData import RadarData


THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class DummyEvent:
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
        self.assertEqual(ip.y, 'nmo_depth')
        with self.assertRaises(ValueError):
            ip = InteractivePicker(data, ydat='dum')

        with self.assertRaises(ValueError):
            ip = InteractivePicker(data, ydat='elev')
        data.elevation = np.arange(ip.dat.tnum)
        data.flags.elev = True
        ip = InteractivePicker(data, ydat='elev')
        self.assertEqual(ip.y, 'elev')

        ip = InteractivePicker(data, x_range=None)
        self.assertEqual(ip.x_range, (0, ip.dat.tnum))

    def test_PickNum(self):
        self.ip.pickNumberBox.setValue(1)

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
        self.assertTrue(self.im.get_cmap(), 'bone')
        self.ip._color_select('CEGSIC')
        self.assertTrue(self.im.get_cmap(), 'CEGSIC')


@unittest.skipIf(not qt, 'No Qt')
class TestVBP(unittest.TestCase):

    def test_VBPInputDialog(self):
        vbp = VBPInputDialog()
        vbp.clickOK()
        self.assertTrue(vbp.lims == (50, 250))
        self.assertTrue(vbp.accepted)

        vbp = VBPInputDialog()
        vbp.minspin.setValue(2)
        vbp.maxspin.setValue(298)
        vbp.clickOK()
        self.assertTrue(vbp.lims == (2, 298))
        self.assertTrue(vbp.accepted)

        vbp = VBPInputDialog()
        vbp.minspin.setValue(299)
        vbp.maxspin.setValue(298)
        # click OK twice since we have bad lims
        vbp.clickOK()
        vbp.clickOK()
        self.assertTrue(vbp.lims == (297, 298))
        self.assertTrue(vbp.accepted)


@unittest.skipIf(not qt, 'No Qt')
class TestCrop(unittest.TestCase):

    def test_CropInputDialog(self):
        cid = CropInputDialog()
        cid.clickOK()
        self.assertTrue(cid.accepted)

        cid = CropInputDialog()
        cid.inputtype.setCurrentText('snum')
        self.assertTrue(cid.spinnerlabel.text() == 'Cutoff (sample num):')

        cid.inputtype.setCurrentText('twtt')
        self.assertTrue(cid.spinnerlabel.text() == 'Cutoff in TWTT (usec):')

        cid.inputtype.setCurrentText('depth')
        self.assertTrue(cid.spinnerlabel.text() == 'Cutoff in depth (m):')

        cid.clickOK()
        self.assertTrue(cid.accepted)


if __name__ == '__main__':
    unittest.main()
