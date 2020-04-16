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
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None, pick_colors=None, clims=None, cmap=Any(object), flatten_layer=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], xd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='dist', ydat='twtt', x_range=None, pick_colors=None, clims=None, cmap=Any(object), flatten_layer=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], yd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='depth', x_range=None, pick_colors=None, clims=None, cmap=Any(object), flatten_layer=None)
        mock_plot_rad.reset_called()
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], xd=True, yd=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='dist', ydat='depth', x_range=None, pick_colors=None, clims=None, cmap=Any(object), flatten_layer=None)
        mock_plot_rad.reset_called()

        # Check that we can save
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], xd=True, yd=True, s=True)
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='dist', ydat='depth', x_range=None, pick_colors=None, clims=None, cmap=Any(object), flatten_layer=None)
        mock_plot_rad.reset_called()

    @patch('impdar.lib.plot.plot_traces', returns=[DummyFig(), None])
    def test_plotPLOTTRACES(self, mock_plot_tr):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], tr=0)
        mock_plot_tr.assert_called_with(Any(RadarData), 0, ydat='twtt')

    @patch('impdar.lib.plot.plot_spectrogram', returns=[DummyFig(), None])
    def test_plotPLOTSPECDENSE(self, mock_plot_specdense):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], spectra=(0, 1), window=0, scaling=1)
        mock_plot_specdense.assert_called_with(Any(RadarData), (0, 1), window=0, scaling=1)

    @patch('impdar.lib.plot.plot_ft', returns=[DummyFig(), None])
    def test_plotFT(self, mock_plot_ft):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], ft=True)
        mock_plot_ft.assert_called_with(Any(RadarData))

    @patch('impdar.lib.plot.plot_hft', returns=[DummyFig(), None])
    def test_plotHFT(self, mock_plot_hft):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], hft=True)
        mock_plot_hft.assert_called_with(Any(RadarData))

    @patch('impdar.lib.plot.plot_power', returns=[DummyFig(), None])
    def test_plotPLOTPOWER(self, mock_plot_power):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], power=0)
        mock_plot_power.assert_called_with(Any(RadarData), 0)

    @patch('impdar.lib.plot.plot_radargram', returns=[DummyFig(), None])
    def test_plotLOADGSSI(self, mock_plot_rad):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'test_gssi.DZT')], filetype='gssi')
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None, pick_colors=None, clims=None, cmap=Any(object), flatten_layer=None)

    @patch('impdar.lib.plot.plot_radargram', returns=[DummyFig(), None])
    def test_plotLOADPE(self, mock_plot_rad):
        plot.plot([os.path.join(THIS_DIR, 'input_data', 'test_pe.DT1')], filetype='pe')
        mock_plot_rad.assert_called_with(Any(RadarData), xdat='tnum', ydat='twtt', x_range=None, pick_colors=None, clims=None, cmap=Any(object), flatten_layer=None)


    def test_plotBADINPUT(self):
        with self.assertRaises(ValueError):
            plot.plot([os.path.join(THIS_DIR, 'input_data', 'small_data.mat')], tr=0, power=1)

    def tearDown(self):
        plt.close('all')


class TestPlotTraces(unittest.TestCase):

    def test_plot_traces(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)
        fig, ax = plot.plot_traces(dat, 0)
        fig, ax = plt.subplots()
        plot.plot_traces(dat, 0, fig=fig)
        plot.plot_traces(dat, 0, fig=fig, ax=ax)
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

    def tearDown(self):
        plt.close('all')


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
        dat.picks.power[:] = 10.5
        # works with constant power
        fig, ax = plot.plot_power(dat, 10)
        
        # works with various inputs
        fig, ax = plt.subplots()
        plot.plot_power(dat, 10, fig=fig)
        plot.plot_power(dat, 10, fig=fig, ax=ax)
        plot.plot_power(dat, 10, clims=(-100, 100))

        # works with multiple inputs
        fig, ax = plot.plot_power([dat, dat], 10)

        # works with projected coordinates
        dat.x_coord = np.arange(dat.data.shape[1])
        dat.y_coord = np.arange(dat.data.shape[1])
        fig, ax = plot.plot_power(dat, 10)
        fig, ax = plot.plot_power([dat, dat], 10)

        with self.assertRaises(ValueError):
            fig, ax = plot.plot_power(dat, 0)

        # gets ok lims with variable power?
        dat.picks.power[:, 0] = 1
        fig, ax = plot.plot_power(dat, 10)

    def tearDown(self):
        plt.close('all')


class TestPlotRadargram(unittest.TestCase):

    def test_plot_radargram_figaxin(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)
        fig, ax = plot.plot_radargram(dat)

        fig, ax = plt.subplots()
        fig, ax = plot.plot_radargram(dat, fig=fig, ax=ax)
        fig, ax = plt.subplots()
        fig, ax = plot.plot_radargram(dat, fig=fig)

        # Varying xdata
        fig, ax = plot.plot_radargram(dat, x_range=None)
        fig, ax = plot.plot_radargram(dat, xdat='dist')
        with self.assertRaises(ValueError):
            fig, ax = plot.plot_radargram(dat, xdat='dummy')

        fig, ax = plot.plot_radargram(dat, y_range=None)
        fig, ax = plot.plot_radargram(dat, ydat='depth')
        with self.assertRaises(ValueError):
            fig, ax = plot.plot_radargram(dat, ydat='dummy')

        # Elevation offsets
        with self.assertRaises(ValueError):
            plot.plot_radargram(dat, ydat='elev')
        dat.flags.elev = True
        dat.elev = np.zeros(dat.data.shape[1])
        dat.elev[1:] = 1
        plot.plot_radargram(dat, ydat='elev')

    def test_plot_radargram_flattenlayer(self):
        dat = NoInitRadarData(big=True)
        dat.picks = Picks(dat)
        dat.picks.add_pick(10)
        dat.picks.power[:] = 10
        dat.picks.samp1[:] = 0
        dat.picks.samp2[:] = 1  # make sure no bugs if this is actually constant
        dat.picks.samp3[:] = 3
        # works with constant power
        fig, ax = plot.plot_radargram(dat, flatten_layer=10)

        # make sure we can actually follow a variable layer
        dat.picks.samp2[:, 1:] = 2
        dat.picks.samp2[:, -1] = 4
        # works with constant power
        fig, ax = plot.plot_radargram(dat, flatten_layer=10)

        dat.picks.samp2[:] = 0  # make sure no bugs if this is at the top
        fig, ax = plot.plot_radargram(dat, flatten_layer=10)

        dat.picks.samp2[:] = dat.data.shape[0] - 1  # make sure no bugs if this is at the bottom
        fig, ax = plot.plot_radargram(dat, flatten_layer=10)

        dat.picks.samp2[:, 1] = np.NaN  # make sure no bugs if this is at the bottom
        fig, ax = plot.plot_radargram(dat, flatten_layer=10)

        with self.assertRaises(ValueError):
            fig, ax = plot.plot_radargram(dat, flatten_layer=1)

    def tearDown(self):
        plt.close('all')


class TestPlotFT(unittest.TestCase):
    
    def test_plot_ft(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)
        fig, ax = plot.plot_ft(dat)
        fig, ax = plt.subplots()
        fig, ax = plot.plot_ft(dat, fig=fig, ax=ax)
        fig, ax = plt.subplots()
        fig, ax = plot.plot_ft(dat, fig=fig)

    def tearDown(self):
        plt.close('all')


class TestPlotHFT(unittest.TestCase):
    
    def test_plot_hft(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)
        fig, ax = plot.plot_hft(dat)
        fig, ax = plt.subplots()
        fig, ax = plot.plot_hft(dat, fig=fig, ax=ax)
        fig, ax = plt.subplots()
        fig, ax = plot.plot_hft(dat, fig=fig)

    def tearDown(self):
        plt.close('all')


class TestPlotPicks(unittest.TestCase):
    
    def test_plot_picks_via_radargram(self):
        """We want to be able to call this via plot_radargram"""
        dat = NoInitRadarData(big=True)
        dat.picks = Picks(dat)
        dat.picks.samp1 = np.ones((2, len(dat.lat)))
        dat.picks.samp2 = np.ones((2, len(dat.lat)))
        dat.picks.samp3 = np.ones((2, len(dat.lat)))

        fig, ax = plot.plot_radargram(dat, pick_colors='mgm')

    def test_plot_picks(self):
        # Only checking that these do not throw errors
        dat = NoInitRadarData(big=True)

        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time)

        dat.picks = Picks(dat)
        dat.picks.samp1 = np.ones((2, len(dat.lat)))
        dat.picks.samp2 = np.ones((2, len(dat.lat)))
        dat.picks.samp3 = np.ones((2, len(dat.lat)))

        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time)
        fig, ax = plt.subplots()
        plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, fig=fig)
        plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, fig=fig, ax=ax)
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors='g')
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors='gmm')
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=['c', 'g'])
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=['cmy', 'brb'])
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=['cm', 'br'])
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=True)
        fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=False)
        with self.assertRaises(ValueError):
            fig, ax = plot.plot_picks(dat, np.arange(int(dat.tnum)), dat.travel_time, colors=['c', 'm', 'b'])

    def tearDown(self):
        plt.close('all')


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

        # no error if freq high
        plot.plot_spectrogram(dat, 100)

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

    def tearDown(self):
        plt.close('all')


if __name__ == '__main__':
    unittest.main()
