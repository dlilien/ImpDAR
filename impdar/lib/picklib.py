#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
This file contains functions that are a for the mechanics of picking, not for the display (e.g. the wavelet stuff)
"""
import numpy as np


def pick(traces, snum_start, snum_end, pickparams):
    # This is similar to stp_pickloop
    picks_out = np.zeros((5, traces.shape[1]))
    dmid = midpoint(traces.shape[1], snum_start, snum_end)
    for i in range(traces.shape[1]):
        pickpacket = packet_pick(traces[:, i], pickparams, dmid[i])
        picks_out[:, i] = pickpacket
    return picks_out


def midpoint(len_tnums, snum_start, snum_end):
    if snum_start == -9999:
        snum_start = snum_end
    return np.round(np.arange(len_tnums) * (snum_end - snum_start) / len_tnums) + snum_start


def packet_power(trace, plength, midpoint):
    # Calculate top snum from plength and midpoint
    topsnum = int(midpoint - ((plength / 2.) - 1.))

    # Calculate bottom snum from plenght and midpoint
    bottom = int(midpoint + (plength / 2.))

    # read off power values from top snum to bottom snum
    powerpacket = trace[topsnum:bottom]
    return powerpacket, topsnum


def packet_pick(trace, pickparams, midpoint):
    powerpacket, topsnum = packet_power(trace, pickparams.plength, midpoint)
    
    # Find the center peak
    cpeak = int(np.argmax(powerpacket[pickparams.scst: pickparams.scst + pickparams.FWW] * pickparams.pol) + pickparams.scst - 1)

    # Find a peak with opposite polarity higher up
    if cpeak > pickparams.FWW:
        tpeak = int(np.argmin(powerpacket[cpeak - pickparams.FWW:cpeak] * pickparams.pol)) + (cpeak - pickparams.FWW - 1)
    else:
        tpeak = int(np.argmin(powerpacket[0:cpeak] * pickparams.pol))

    # Find a peak with opposite polarity lower down
    if cpeak + pickparams.FWW < pickparams.plength:
        bpeak = int(np.argmin(powerpacket[cpeak:cpeak + pickparams.FWW] * pickparams.pol)) + cpeak - 1
    else:
        bpeak = int(np.argmin(powerpacket[cpeak:pickparams.plength] * pickparams.pol)) + cpeak - 1
    power = np.sum(powerpacket[tpeak:bpeak] ** 2. / (bpeak - tpeak + 1))
    return [tpeak + topsnum - 1, cpeak + topsnum - 1, bpeak + topsnum - 1, np.nan, power]
