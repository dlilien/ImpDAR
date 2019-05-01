#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
This file contains functions that are a for the mechanics of picking, not for the display (e.g. the wavelet stuff)
"""
import numpy as np
from scipy.spatial import cKDTree as KDTree


def pick(traces, snum_start, snum_end, pickparams):
    """Pick a reflector in some traces. Uses a line between the starting and ending picks to guide it.

    Parameters
    ----------
    traces: np.ndarray
        The chunk of data we are picking. Should be a slice of the overall data with shape snum x (size to pick)
    snum_start: int
        The index of the pick in the leftmost trace. We would normally get this from the last pick.
    snum_end: int
        The index of the pick in the rightmost trace.
    pickparams: `impdar.lib.PickParameters.PickParameters`
        Use for polarity, frequency, plength.
    """
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

    if len(trace.shape) > 1:
        raise ValueError('Need a single, flat trace')
    # Calculate top snum from plength and midpoint
    topsnum = int(midpoint - ((plength / 2.)))

    # Calculate bottom snum from plenght and midpoint
    bottom = int(midpoint + (plength / 2.))

    # read off power values from top snum to bottom snum
    powerpacket = trace[topsnum:bottom]
    return powerpacket, topsnum


def packet_pick(trace, pickparams, midpoint):
    powerpacket, topsnum = packet_power(trace, pickparams.plength, midpoint)

    # Check if we are taking a bad slice
    if len(powerpacket[pickparams.scst: pickparams.scst + pickparams.FWW]) == 0:
        raise ValueError('Your choice of frequency (too low) is causing the pick window to be too large')

    # Find the center peak
    cpeak = int(np.argmax(powerpacket[pickparams.scst: pickparams.scst + pickparams.FWW] * pickparams.pol) + pickparams.scst)

    # Find a peak with opposite polarity higher up
    if cpeak > pickparams.FWW:
        tpeak = int(np.argmin(powerpacket[cpeak - pickparams.FWW:cpeak] * pickparams.pol)) + (cpeak - pickparams.FWW)
    else:
        tpeak = int(np.argmin(powerpacket[:cpeak] * pickparams.pol))

    # Find a peak with opposite polarity lower down
    if cpeak + pickparams.FWW < pickparams.plength:
        bpeak = int(np.argmin(powerpacket[cpeak:cpeak + pickparams.FWW] * pickparams.pol)) + cpeak
    else:
        # I can't seem to hit this line in tests. Might be due to offset between matlab and python
        bpeak = int(np.argmin(powerpacket[cpeak:] * pickparams.pol)) + cpeak
    power = np.sum(powerpacket[tpeak:bpeak] ** 2. / (bpeak - tpeak + 1))

    return [tpeak + topsnum, cpeak + topsnum, bpeak + topsnum, np.nan, power]


def get_intersection(data_main, data_cross):
    if data_cross.picks is None or data_cross.picks.picknums is None or len(data_cross.picks.picknums) == 0 or data_cross.picks.samp1 is None:
        raise AttributeError('We do not have viable cross picks')

    out_tnums = np.zeros_like(data_cross.picks.picknums, dtype=int)
    out_sns = np.zeros_like(data_cross.picks.picknums, dtype=int)

    tree = KDTree(np.vstack((data_main.x_coord.flatten(), data_main.y_coord.flatten())).transpose())
    for i in range(len(out_tnums)):
        mask_pick_not_nan = ~np.isnan(data_cross.picks.samp1[i])
        closest_dist, closest_inds = tree.query(np.vstack((data_cross.x_coord[mask_pick_not_nan].flatten(), data_cross.y_coord[mask_pick_not_nan].flatten())).transpose())

        # need the spot in the cross profile that is closest
        ind_dat_cross = np.argmin(closest_dist)

        # Where to plot this on the main profile
        out_tnums[i] = closest_inds[ind_dat_cross]

        out_sns[i] = data_cross.picks.samp1[i, :][mask_pick_not_nan][ind_dat_cross].astype(int)

    return out_tnums, out_sns
