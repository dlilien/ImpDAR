#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Test the machinery of plotting. We will not try the "show" lines.
"""
import sys
import os
import unittest
import numpy as np
from impdar.lib.RadarData import RadarData
from impdar.lib import plot
if sys.version_info[0] >= 3:
    from unittest.mock import patch, MagicMock
else:
    from mock import patch, MagicMock

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class NoInitRadarData(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        self.data = np.zeros((10, 20))
        # need to set this to avoid divide by zero later
        self.dt = 1
        self.dist = np.arange(20)
        self.tnum = 20
        self.trace_num = np.arange(self.tnum) + 1.
        self.snum = 10
        self.travel_time = np.arange(10)


class DummyFig:
    # to mock saving
    def __init__(self):
        sfcalled = False

    def savefig(self, fn, dpi=None, ftype=None):
        sfcalled = True


def Any(cls):
    # to mock data argument in tests
    class Any(cls):
        def __init__(self):
            pass

        def __eq__(self, other):
            return True
    return Any()


class TestPlot(unittest.TestCase):
    
    def test_plot_traces(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData()
        fig, ax = plot.plot_traces(dat, 0)
        fig, ax = plot.plot_traces(dat, [1, 1])
        fig, ax = plot.plot_traces(dat, [1, 18])
        with self.assertRaises(TypeError):
            fig, ax = plot.plot_traces(dat, np.arange(10))
        with self.assertRaises(IndexError):
            fig, ax = plot.plot_traces(dat, 999)

        # no nmo
        fig, ax = plot.plot_traces(dat, 0, ydat='depth')

        # with nmo
        dat.nmo_depth = np.arange(10)
        fig, ax = plot.plot_traces(dat, 0, ydat='depth')
        with self.assertRaises(ValueError):
            fig, ax = plot.plot_traces(dat, 0, ydat='dum')

    @patch('impdar.lib.plot.plot_radargram', returns=[DummyFig(), None])
    def test_plotPLOTARGS(self, mock_plot_rad):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')])
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], xd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='dist', ydat='twtt', x_range=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], yd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='depth', x_range=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], xd=True, yd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='dist', ydat='depth', x_range=None)
        mock_plot_rad.reset_called()

    @patch('impdar.lib.plot.plot_traces', returns=[DummyFig(), None])
    def test_plotPLOTTRACES(self, mock_plot_tr):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], tr=0)
        mock_plot_tr.assert_called_with(Any(RadarData), 0, ydat='twtt')

    @patch('impdar.lib.plot.plot_radargram', returns=[DummyFig(), None])
    def test_plotLOADGSSI(self, mock_plot_rad):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT')], gssi=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None)

    @patch('impdar.lib.plot.plot_radargram', returns=[DummyFig(), None])
    def test_plotLOADPE(self, mock_plot_rad):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'test_pe.DT1')], pe=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None)

    def test_plot(self):
        pass


if __name__ == '__main__':
    unittest.main()
