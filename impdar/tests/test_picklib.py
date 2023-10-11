#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Test the methods on the Pick object
"""

import os
import unittest
import numpy as np
from impdar.lib.NoInitRadarData import NoInitRadarData
from impdar.lib import picklib, Picks, RadarData

traces = np.random.random((300, 200))

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class BareRadarData(NoInitRadarData):

    def __init__(self):
        super(BareRadarData, self).__init__()
        self.dt = 1.0e-7
        self.data = traces
        self.snum, self.tnum = self.data.shape
        self.picks = Picks.Picks(self)


class TestPickLib(unittest.TestCase):

    def test_midpoint(self):
        # Fully test that we find midpoints as expected
        self.assertTrue(np.allclose(picklib._midpoint(200, 100, 100), np.ones((200,)) * 100.))
        self.assertTrue(np.allclose(picklib._midpoint(200, -9999, 100), np.ones((200,)) * 100.))
        self.assertTrue(np.allclose(picklib._midpoint(200, 0, 200), np.arange(200)))

    def test_packet_power(self):
        with self.assertRaises(ValueError):
            picklib.packet_power(traces, 2, 100)
        self.assertTrue(len(picklib.packet_power(traces[:, 0], 10, 100)[0]) == 10)
        self.assertTrue(len(picklib.packet_power(traces[:, 0], 11, 100)[0]) == 11)

    def test_packet_pick(self):
        easy_pick_trace = np.zeros((traces.shape[1], ))
        cpeak = 100
        bpeak = -200
        tpeak = -100
        easy_pick_trace[101] = cpeak
        easy_pick_trace[107] = bpeak
        easy_pick_trace[95] = tpeak
        data = BareRadarData()
        # do something ill-advised where we now have mismatched plength, scst, and FWW
        data.picks.pickparams.scst = 200
        data.picks.pickparams.FWW = 200
        with self.assertRaises(ValueError):
            picklib.packet_pick(traces[:, 0], data.picks.pickparams, 100)

        # This should also be an error due to mismatched plength, scst, and FWW
        data.picks.pickparams.scst = 2
        data.picks.pickparams.FWW = 0
        with self.assertRaises(ValueError):
            picklib.packet_pick(traces[:, 0], data.picks.pickparams, 100)

        # We should be able to pick within 3
        # for a variety of frequency this is how the discrete rounding works out
        data = BareRadarData()
        for freq in [0.85, 0.9, 0.95]:
            data.picks.pickparams.freq_update(freq)
            for pick in [98, 101, 104]:
                pickout = picklib.packet_pick(easy_pick_trace, data.picks.pickparams, pick)
                self.assertEqual(pickout[0], 95)
                self.assertEqual(pickout[1], 101)
                self.assertEqual(pickout[2], 107)

        # and now we should be able to pick slightly wider
        data.picks.pickparams.freq_update(0.8)
        for pick in [97, 101, 105]:
            pickout = picklib.packet_pick(easy_pick_trace, data.picks.pickparams, pick)
            self.assertEqual(pickout[0], 95)
            self.assertEqual(pickout[1], 101)
            self.assertEqual(pickout[2], 107)

        # if our plength is really short we should still hit the middle
        # sides are undefined
        data.picks.pickparams.freq_update(4.0)
        pickout = picklib.packet_pick(easy_pick_trace, data.picks.pickparams, 101)
        self.assertEqual(pickout[1], 101)
        pickout = picklib.packet_pick(easy_pick_trace, data.picks.pickparams, 102)
        self.assertEqual(pickout[1], 101)

    def test_pick(self):
        easy_pick_traces = np.zeros_like(traces)
        cpeak = 100
        bpeak = -200
        tpeak = -100
        easy_pick_traces[101, :] = cpeak
        easy_pick_traces[107, :] = bpeak
        easy_pick_traces[95, :] = tpeak

        data = BareRadarData()
        data.picks.pickparams.freq_update(1.0)

        # first just do a line across the middle guessing correctly
        picks = picklib.pick(easy_pick_traces, 101, 101, data.picks.pickparams)
        self.assertTrue(np.all(picks[0, :] == 95))
        self.assertTrue(np.all(picks[1, :] == 101))
        self.assertTrue(np.all(picks[2, :] == 107))

        # now do a line across the middle guessing slanty
        picks = picklib.pick(easy_pick_traces, 99, 103, data.picks.pickparams)
        self.assertTrue(np.all(picks[0, :] == 95))
        self.assertTrue(np.all(picks[1, :] == 101))

    def test_intersection(self):
        thisdata = RadarData.RadarData(os.path.join(THIS_DIR, 'input_data', 'along_picked.mat'))
        thatdata = RadarData.RadarData(os.path.join(THIS_DIR, 'input_data', 'cross_picked.mat'))

        nopickdata = RadarData.RadarData(os.path.join(THIS_DIR, 'input_data', 'small_data.mat'))
        tnum, sn = picklib.get_intersection(thisdata, thatdata, multiple_int=True, return_nans=False)
        self.assertTrue(len(sn) == len(thatdata.picks.picknums))
        tnum, sn = picklib.get_intersection(thatdata, thisdata, multiple_int=False, return_nans=False)
        self.assertTrue(len(sn) == len(thisdata.picks.picknums))
        tnum, sn = picklib.get_intersection(thatdata, thisdata, multiple_int=True, return_nans=True)
        self.assertTrue(len(sn) == len(thisdata.picks.picknums))
        tnum, sn = picklib.get_intersection(thatdata, thisdata, multiple_int=False, return_nans=True)
        self.assertTrue(len(sn) == len(thisdata.picks.picknums))
        with self.assertRaises(AttributeError):
            tnum, sn = picklib.get_intersection(thisdata, nopickdata)

        thatdata.picks.samp1[:, :] = np.NaN
        thatdata.picks.samp2[:, :] = np.NaN
        thatdata.picks.samp3[:, :] = np.NaN
        tnum, sn = picklib.get_intersection(thisdata, thatdata, multiple_int=False, return_nans=True)
        self.assertTrue(len(sn) == len(thisdata.picks.picknums))
        tnum, sn = picklib.get_intersection(thisdata, thatdata, multiple_int=False, return_nans=False)
        self.assertTrue(len(sn) == len(thisdata.picks.picknums))

    def test_autopick(self):
        easy_pick_traces = np.zeros_like(traces)
        cpeak = 100
        bpeak = -200
        tpeak = -100
        easy_pick_traces[101, :] = cpeak
        easy_pick_traces[107, :] = bpeak
        easy_pick_traces[95, :] = tpeak

        data = BareRadarData()
        data.picks.pickparams.freq_update(1.0)
        data.data = easy_pick_traces

        # first just do a line across the middle guessing correctly
        picks = picklib.auto_pick(data, [101], [9])
        self.assertTrue(np.all(picks[0, 0, :] == 95))
        self.assertTrue(np.all(picks[0, 1, :] == 101))
        self.assertTrue(np.all(picks[0, 2, :] == 107))

        with self.assertRaises(ValueError):
            picks = picklib.auto_pick(data, [101, 10], [9])

        with self.assertRaises(TypeError):
            picks = picklib.auto_pick(data, 101, 9)


if __name__ == '__main__':
    unittest.main()
