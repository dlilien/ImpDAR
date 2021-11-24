#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""Functions that are a for the mechanics of picking, not for the display."""

import numpy as np
from scipy.spatial import cKDTree as KDTree



def pick(traces, snum_start, snum_end, pickparams):
    """Pick a reflector in some traces.

    Uses a line between the starting and ending picks to guide picking the
    maximum (or minimum) return, and the surrounding peaks with opposite
    polarity.

    Parameters
    ----------
    traces: numpy.ndarray
        The chunk of data we are picking. Should be a slice of the overall
        data with shape snum x (size to pick)
    snum_start: int
        The index of the pick in the leftmost trace. We would normally get
        this from the last pick.
    snum_end: int
        The index of the pick in the rightmost trace.
    pickparams: `impdar.lib.PickParameters.PickParameters`
        Use for polarity, frequency, plength.

    Returns
    -------
    numpy.ndarray
        The picks selected. Rows are: top of packet, center pick, bottom of
        packet, time (deprecated, all nans), and power. Size 5xtnum
    """
    # This is similar to stp_pickloop
    picks_out = np.zeros((5, traces.shape[1]))
    dmid = _midpoint(traces.shape[1], snum_start, snum_end)
    for i in range(traces.shape[1]):
        pickpacket = packet_pick(traces[:, i], pickparams, dmid[i])
        picks_out[:, i] = pickpacket
    return picks_out


def auto_pick(dat,snums,tnums):
    """Automatically pick any number of reflectors.

    Parameters
    ----------
    dat: object class
        data object
    indices: numpy.ndarray
        indices on the left side of the image where the picks will start from
        These are the centerpoint of the wavelet

    Returns
    -------
    numpy.ndarray
        The picks selected. Rows are: top of packet, center pick, bottom of
        packet, time (deprecated, all nans), and power. Size 5xtnum
    """

    picks_out = np.empty((len(snums),5,dat.tnum))

    for i in range(len(snums)):
        j = int(tnums[i])
        t_start = int(tnums[i])
        dmid = snums[i]
        for n in range(dat.tnum):
            pp = packet_pick(dat.data[:,j],dat.picks.pickparams,dmid)
            picks_out[i,:,j] = pp
            if j <= t_start and j > 0:
                dmid = (pp[0]+pp[2])//2
                j -= 1
            elif j == 0:
                dmid = (picks_out[i,0,t_start]+picks_out[i,2,t_start])//2
                j = t_start + 1
            elif j > t_start:
                dmid = (pp[0]+pp[2])//2
                j += 1

    return picks_out


def _midpoint(len_tnums, snum_start, snum_end):
    if snum_start == -9999:
        snum_start = snum_end
    return np.round(np.arange(len_tnums) * (
        snum_end - snum_start) / len_tnums) + snum_start


def packet_power(trace, plength, midpoint):
    """Return power of a packet.

    This function is pretty boring.
    It just finds a window around a point in a trace.

    Parameters
    ----------
    trace: numpy.ndarray
        (snum,) The trace in which to find the window
    plength: float
        The size of the packet (in samples)
    midpoint: int
        The central sample

    Returns
    -------
    numpy.ndarray
        The packet of length plength
    int
        The index of the top sample (desirable for calculating overall indices
        in functions that call this one).
    """
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
    """Really do the picking.

    This is where we look for the highest amplitude return, and the opposite
    polarity returns surrounding it

    Parameters
    ----------
    trace: 1d numpy.ndarray
        The trace in which we are looking for our return
    pickparams: impdar.lib.PickParameters.PickParameters
        The information about picking that we need for determining window
        size and polarity
    midpoint: int
        The guess at the index where a pick is, used for searching out the
        highest amplitude return

    Returns
    -------
    list
        len=5. Top of packet, middle of packet, bottom of packet, nan, power
    """
    powerpacket, topsnum = packet_power(trace, pickparams.plength, midpoint)

    # Check if we are taking a bad slice
    if len(powerpacket) < pickparams.scst + pickparams.FWW:
        raise ValueError('Your choice of frequency is too high, \
                         making the pick window sub-pixel in size')
    if len(powerpacket[pickparams.scst:
                       pickparams.scst + pickparams.FWW]) == 0:
        raise ValueError('Your choice of frequency (too low) is causing the\
                         pick window to be too large')

    # Find the center peak
    cpeak = int(np.argmax(powerpacket[pickparams.scst + 1:
                                      pickparams.scst + pickparams.FWW + 1]
                          * pickparams.pol) + pickparams.scst + 1)

    # Find a peak with opposite polarity higher up
    if cpeak > pickparams.FWW:
        tpeak = int(np.argmin(powerpacket[cpeak - pickparams.FWW:
                                          cpeak] * pickparams.pol)) + (
                                              cpeak - pickparams.FWW)
    elif cpeak <= 1:
        tpeak = 0
    else:
        tpeak = int(np.argmin(powerpacket[:cpeak] * pickparams.pol))

    # Find a peak with opposite polarity lower down
    if cpeak + pickparams.FWW < pickparams.plength:
        bpeak = int(np.argmin(powerpacket[cpeak + 1:
                                          cpeak + pickparams.FWW + 1]
                              * pickparams.pol)) + cpeak + 1
    elif cpeak >= pickparams.plength - 1:
        bpeak = pickparams.plength - 1
    else:
        bpeak = int(np.argmin(powerpacket[cpeak + 1:] * pickparams.pol)
                    ) + cpeak + 1
    power = np.sum(powerpacket[tpeak:bpeak + 1] ** 2.) / (bpeak - tpeak + 1)

    return [tpeak + topsnum, cpeak + topsnum, bpeak + topsnum, np.nan, power]


def get_intersection(data_main, data_cross, multiple_int=True, return_nans=False,
                     cutoff=10.0):
    """Find the intersection of two radar datasets.

    Used for plotting up where pick depths at places where two profiles cross.
    This is a pretty simple function as implemented,
    so it is not going to find you multiple intersections.
    Rather, you are just going to get the single point where the traces are
    closest. Note that this will work picking sequential profiles too.
    If there are nans at the intersection, we look for the closest non-nan.

    Parameters
    ----------
    data_main: impdar.lib.RadarData.RadarData
        The RadarData object upon which you want the intersection in reference
        to (i.e. then one you will plot up)
    data_cross: impdar.lib.RadarData.RadarData
        The crossprofile from which you are going to plot up layers.
        Must have some picks in it.
    return_nans: bool, optional
        Return the closest sample, even if it was nanpicked.
        Default is false (find closest non-nan value)
    cutoff: float, optional
        The maximum distance for multiple intersections.

    Returns
    -------
    np.ndarray (npicks,) or (npicks, m)
        The tracenumber in the main profile at which we are getting the sample
        from the crossprofile
    np.ndarray (npicks,) or (npicks, m)
        The depths (in sample number) of the layers in the cross profile.
        Note that this essentially assume you are using the same snum/depth
        conversion between the two profiles.

    Raises
    ------
    AttributeError:
        if there are no picks in the cross profile.
    """
    if data_cross.picks is None or data_cross.picks.picknums is None or len(
            data_cross.picks.picknums) == 0 or data_cross.picks.samp1 is None:
        raise AttributeError('We do not have viable cross picks')

    tree = KDTree(np.vstack((
        data_main.x_coord.flatten(), data_main.y_coord.flatten())).transpose())

    if multiple_int:
        # Get the maximum possible size of output
        dist, _ = tree.query(np.vstack(
            (data_cross.x_coord.flatten(),
             data_cross.y_coord.flatten())).transpose())
        maxn = np.sum(dist < cutoff)
        out_tnums = np.zeros((len(data_cross.picks.picknums), maxn), dtype=float)
        out_sns = np.zeros((len(data_cross.picks.picknums), maxn), dtype=float)
        out_tnums[:, :] = np.nan
        out_sns[:, :] = np.nan
        for i, _ in enumerate(out_tnums):
            mask_pick_not_nan = ~np.isnan(data_cross.picks.samp1[i])
            dist, inds = tree.query(np.vstack(
                (data_cross.x_coord[mask_pick_not_nan].flatten(),
                 data_cross.y_coord[mask_pick_not_nan].flatten())).transpose())
            maxn = np.sum(dist < cutoff)
            out_tnums[i, :maxn] = inds[dist < cutoff]
            out_sns[i, :maxn] = data_cross.picks.samp1[i, :][mask_pick_not_nan][dist < cutoff]
    else:
        out_tnums = np.zeros_like(data_cross.picks.picknums, dtype=float)
        out_sns = np.zeros_like(data_cross.picks.picknums, dtype=float)

        for i, _ in enumerate(out_tnums):
            if return_nans:
                mask_pick_not_nan = np.ones_like(
                    data_cross.picks.samp1[i], dtype=bool)
            else:
                mask_pick_not_nan = ~np.isnan(data_cross.picks.samp1[i])
            dist, inds = tree.query(np.vstack(
                (data_cross.x_coord[mask_pick_not_nan].flatten(),
                 data_cross.y_coord[mask_pick_not_nan].flatten())).transpose())

            # need the spot in the cross profile that is closest
            # sequence will be empty if we have a pick that is purely nans
            if len(dist) > 0:
                ind_dat_cross = np.argmin(dist)

                # Where to plot this on the main profile
                out_tnums[i] = inds[ind_dat_cross]

                out_sns[i] = data_cross.picks.samp1[i, :][
                    mask_pick_not_nan][ind_dat_cross].astype(int)
            else:
                out_tnums[i] = np.nan
                out_sns[i] = np.nan

    return out_tnums, out_sns
