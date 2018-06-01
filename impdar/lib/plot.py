#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import os.path
import numpy as np
import matplotlib.pyplot as plt
from .load import load


def plot(fn, gssi=False, pe=False, s=False, ftype='png', dpi=300, *args, **kwargs):
    if gssi and pe:
        raise ValueError('Input cannot be both pulse-ekko and gssi')
    if gssi:
        radar_data = load('gssi', fn)
    elif pe:
        radar_data = load('pe', fn)
    else:
        radar_data = load('mat', fn)

    figs = [plot_radargram(dat) for dat in radar_data]
    if s:
        [f.savefig(os.path.splitext(fn0)[0] + '.' + ftype, dpi=dpi) for f, fn0 in zip(figs, fn)]
    else:
        plt.show()


def plot_radargram(dat):
    fig, ax = plt.subplots(figsize=(12, 8))
    X, Y = np.meshgrid(dat.trace_num, dat.travel_time)
    lims = np.percentile(dat.data, (10, 90))

    Y *= -1
    # plt.imshow(data, cmap=plt.cm.gray_r, norm=LogNorm(vmin=data.min(), vmax=data.max()))
    # plt.imshow(data, cmap=plt.cm.bwr, norm=LogNorm(vmin=1.0e-6, vmax=data.max()))
    # plt.imshow(data, cmap=plt.cm.gray_r, vmin=data.min(), vmax=data.max())
    ax.imshow(dat.data, cmap=plt.cm.gray_r, vmin=lims[0], vmax=lims[1])
    return fig
