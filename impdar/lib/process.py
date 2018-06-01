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


def process_and_exit(fn, vbp=None, hfilt=None, gssi=False, pe=False, **kwargs):
    if gssi and pe:
        raise ValueError('Input cannot be both pulse-ekko and gssi')
    if gssi:
        radar_data = load('gssi', fn)
    elif pe:
        radar_data = load('pe', fn)
    else:
        radar_data = load('mat', fn)

    done_stuff = False

    if vbp is not None:
        for dat in radar_data:
            dat.vertical_band_pass(*vbp)
        done_stuff = True

    if hfilt is not None:
        for dat in radar_data:
            dat.horizontal_band_pass(*hfilt)
        done_stuff = True

    if not done_stuff:
        print('No processing steps performed. Not saving!')
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
