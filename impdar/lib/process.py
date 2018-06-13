#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
This is the overarching function that collects information about the processing steps we want to perform
"""
import os.path
from .load import load


def process_and_exit(fn, rev=False, vbp=None, hfilt=None, ahfilt=False, gssi=False, pe=False, **kwargs):
    """Perform one or more processing steps, save, and exit

    Parameters
    ----------
    fn: list of strs
        The filename(s) to process. Assumed to be .mat files unless gssi or pe is True.
    rev: bool, optional
        Reverse the profile orientation. Default is False.
    vbp: 2-tuple, optional
        Vertical bandpass between (vbp1, vbp2) MHz. Default None (no filtering).
    hfilt: 2-tuple, optional
        Horizontal filter subtracting average trace between (hfilt1, hfilt2). Default is None (no hfilt).
    ahfilt: bool, optional
        Adaptively horizontally filter the data.
    gssi: bool, optional
        If True, expect a gssi radar file as input rather than a .mat file.
    pe: bool, optional
        If True, expect a Pulse Ekko radar file as input rather than a .mat file.
    """
        
    if gssi and pe:
        raise ValueError('Input cannot be both pulse-ekko and gssi')
    if gssi:
        radar_data = load('gssi', fn)
    elif pe:
        radar_data = load('pe', fn)
    else:
        radar_data = load('mat', fn)

    processed = process(radar_data, rev=rev, vbp=vbp, hfilt=hfilt, ahfilt=ahfilt, **kwargs)
    if not processed:
        return

    if 'o' in kwargs and kwargs['o'] is not None:
        if len(radar_data) > 1:
            for d, f in zip(radar_data, fn):
                out_fn = kwargs['o'] + '/' + os.path.split(os.path.splitext(f)[0])[1] + '_proc.mat'
                d.save(out_fn)
        else:
            radar_data[0].save(out_fn)
    else:
        for d, f in zip(radar_data, fn):
            out_fn = os.path.splitext(f)[0] + '_proc.mat'
            d.save(out_fn)


def process(RadarDataList, rev=False, vbp=None, hfilt=None, ahfilt=False, nmo=None, crop=None, **kwargs):
    """Perform one or more processing steps on a list of RadarData objects

    Parameters
    ----------
    RadarDataList: list of strs
        The `~impdar.RadarData` objects to process
    rev: bool, optional
        Reverse the profile orientation. Default is False.
    vbp: 2-tuple, optional
        Vertical bandpass between (vbp1, vbp2) MHz. Default None (no filtering).
    hfilt: 2-tuple, optional
        Horizontal filter subtracting average trace between (hfilt1, hfilt2). Default is None (no hfilt).
    ahfilt: bool, optional
        Adaptively horizontally filter the data.

    Returns
    -------
    processed: bool
        If True, we did something, if False we didn't
    """

    # first some argument checking so we don't crash later
    if crop is not None:
        try:
            if crop[0] not in ['top', 'bottom']:
                raise ValueError('First element of crop must be in ["top", "bottom"]')
            if crop[1] not in ['snum', 'twtt', 'depth']:
                raise ValueError('Second element of crop must be in ["snum", "twtt", "depth"]')
            try:
                crop[2] = float(crop[2])
            except ValueError:
                raise ValueError('Third element of crop must be convertible to a float')
        except IndexError:
            raise TypeError('Crop must be subscriptible')

    done_stuff = False
    if rev:
        for dat in RadarDataList:
            dat.reverse()

    if vbp is not None:
        for dat in RadarDataList:
            dat.vertical_band_pass(*vbp)
        done_stuff = True

    if hfilt is not None:
        for dat in RadarDataList:
            dat.hfilt(ftype='hfilt', bounds=hfilt)
        done_stuff = True

    if ahfilt:
        for dat in RadarDataList:
            dat.hfilt(ftype='adaptive')
        done_stuff = True

    if nmo is not None:
        if type(nmo) == float:
            print('One nmo value given. Assuming that this is the separation. Uice=1.6')
            nmo = (nmo, 1.6)
        for dat in RadarDataList:
            dat.nmo(*nmo)
        done_stuff = True

    # Crop after nmo so that we have nmo_depth available for cropping if desired
    if crop is not None:
        for dat in RadarDataList:
            dat.crop(*crop)
        done_stuff = True

    if not done_stuff:
        print('No processing steps performed. Not saving!')
        return False
    return True
