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
from impdar.lib.NoInitRadarData import NoInitRadarData
from impdar.lib.Picks import Picks
from impdar.lib import plot
import matplotlib.pyplot as plt
if sys.version_info[0] >= 3:
    from unittest.mock import patch
else:
    from mock import patch

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


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
    @patch('impdar.lib.plot.plot_radargram', returns=[DummyFig(), None])
    def test_plotPLOTARGS(self, mock_plot_rad):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')])
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None, pick_colors=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], xd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='dist', ydat='twtt', x_range=None, pick_colors=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], yd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='depth', x_range=None, pick_colors=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], xd=True, yd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='dist', ydat='depth', x_range=None, pick_colors=None)
        mock_plot_rad.reset_called()

        # Check that we can save
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], xd=True, yd=True, s=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='dist', ydat='depth', x_range=None, pick_colors=None)
        mock_plot_rad.reset_called()

    @patch('impdar.lib.plot.plot_traces', returns=[DummyFig(), None])
    def test_plotPLOTTRACES(self, mock_plot_tr):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], tr=0)
        mock_plot_tr.assert_called_with(Any(RadarData), 0, ydat='twtt')

    @patch('impdar.lib.plot.plot_power', returns=[DummyFig(), None])
    def test_plotPLOTPOWER(self, mock_plot_power):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], power=0)
        mock_plot_power.assert_called_with(Any(RadarData), 0)

    @patch('impdar.lib.plot.plot_radargram', returns=[DummyFig(), None])
    def test_plotLOADGSSI(self, mock_plot_rad):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT')], filetype='gssi')
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None, pick_colors=None)

    @patch('impdar.lib.plot.plot_radargram', returns=[DummyFig(), None])
    def test_plotLOADPE(self, mock_plot_rad):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'test_pe.DT1')], filetype='pe')
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None, pick_colors=None)


    def test_plotBADINPUT(self):
        with self.assertRaises(ValueError):
            plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], tr=0, power=1)


class TestPlotTraces(unittest.TestCase):

    def test_plot_traces(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)
        fig, ax = plot.plot_traces(dat, 0)
        fig, ax = plot.plot_traces(dat, [1, 1])
        fig, ax = plot.plot_traces(dat, [1, 18])
        with self.assertRaises(ValueError):
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

        # Make sure we handle axes rescaling ok
        dat.data[:, 0] = 10
        dat.data[:, 1] = -10
        fig, ax = plot.plot_traces(dat, (0, 2))


class TestPlotPower(unittest.TestCase):

    def test_plot_power(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)
        with self.assertRaises(TypeError):
            fig, ax = plot.plot_power(dat, [12, 14])
        with self.assertRaises(ValueError):
            fig, ax = plot.plot_power(dat, 0)

        dat.picks = Picks(dat)
        dat.picks.add_pick(10)
        dat.picks.power[:] = 10
        # works with constant power
        fig, ax = plot.plot_power(dat, 10)

        with self.assertRaises(ValueError):
            fig, ax = plot.plot_power(dat, 0)

        # gets ok lims with variable power?
        dat.picks.power[:, 0] = 1
        fig, ax = plot.plot_power(dat, 10)


class TestPlotRadargram(unittest.TestCase):

    def test_plot_radargram_figaxin(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)
        fig, ax = plot.plot_radargram(dat)

        fig, ax = plt.subplots()
        fig, ax = plot.plot_radargram(dat, fig=fig, ax=ax)
        fig, ax = plt.subplots()
        fig, ax = plot.plot_radargram(dat, fig=fig)


class TestPlotPicks(unittest.TestCase):

    def test_plot_picks(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)

        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time)

        dat.picks = Picks(dat)
        dat.picks.samp1 = np.ones((2, len(dat.lat)))
        dat.picks.samp2 = np.ones((2, len(dat.lat)))
        dat.picks.samp3 = np.ones((2, len(dat.lat)))

        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time)
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors='g')
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors='gmm')
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=['c', 'g'])
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=['cmy', 'brb'])
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=['cm', 'br'])
        with self.assertRaises(ValueError):
            fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=['c', 'm', 'b'])


class TestPlotSpectral(unittest.TestCase):
    def test_plot_spectrogram(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)

        dat.picks = Picks(dat)
        dat.picks.samp1 = np.ones((2, len(dat.lat)))
        dat.picks.samp2 = np.ones((2, len(dat.lat)))
        dat.picks.samp3 = np.ones((2, len(dat.lat)))

        fig, ax = plot.plot_spectrogram(dat, (0.,5.0))
        plot.plot_spectrogram(dat, (0.,5.0), fig=fig)
        plot.plot_spectrogram(dat, (0.,5.0), fig=fig, ax=ax)
        plot.plot_spectrogram(dat, (0.,5.0), window='hamming')
        plot.plot_spectrogram(dat, (0.,5.0), scaling='density')

        # freq too low
        with self.assertRaises(ValueError):
            plot.plot_spectrogram(dat, (0.,-100))

        with self.assertRaises(ValueError):
            plot.plot_spectrogram(dat, (0.,5), scaling='dummy')

        with self.assertRaises(ValueError):
            plot.plot_spectrogram(dat, (0.,5), window='dummy')


    @unittest.skipIf(sys.version_info[0] < 3, 'Att error on 2')
    def test_failure_3(self):
        dat = NoInitRadarData(big=True)
        dat.picks = Picks(dat)
        dat.picks.samp1 = np.ones((2, len(dat.lat)))
        dat.picks.samp2 = np.ones((2, len(dat.lat)))
        dat.picks.samp3 = np.ones((2, len(dat.lat)))
        with self.assertRaises(TypeError):
            plot.plot_specdense(dat, 'bad')

    @unittest.skipIf(sys.version_info[0] >= 3, 'Type error on 3')
    def test_failure_3(self):
        dat = NoInitRadarData(big=True)
        dat.picks = Picks(dat)
        dat.picks.samp1 = np.ones((2, len(dat.lat)))
        dat.picks.samp2 = np.ones((2, len(dat.lat)))
        dat.picks.samp3 = np.ones((2, len(dat.lat)))
        with self.assertRaises(AttributeError):
            plot.plot_specdense(dat, 'bad')


if __name__ == '__main__':
    unittest.main()
