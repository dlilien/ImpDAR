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
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
from .load import load


def plot(fn, tr=None, gssi=False, pe=False, s=False, ftype='png', dpi=300, xd=False, yd=False, x_range=(0, -1), *args, **kwargs):
    """We have an overarching function here to handle a number of plot types

    Parameters
    ----------
    fn: list of strs
        A list of filenames to plot individually.
    tr: tuple or int, optional
        Plot traces tr[1] to tr[2] (or trace tr) rather than the radargram. Default is None (plot radargram)
    gssi: bool, optional
        If True, fns are .DZT files
    pe: bool, optional
        If True, fns are Pulse Ekko files
    x_range: tuple, optional
        The range of traces to plot in the radargram. Default is (0, -1) (plot all traces)
    """
    if gssi and pe:
        raise ValueError('Input cannot be both pulse-ekko and gssi')
    if gssi:
        radar_data = load('gssi', fn)
    elif pe:
        radar_data = load('pe', fn)
    else:
        radar_data = load('mat', fn)

    if xd:
        xdat = 'dist'
    else:
        xdat = 'tnum'
    if yd:
        ydat = 'depth'
    else:
        ydat = 'twtt'

    if tr is not None:
        figs = [plot_traces(dat, tr, ydat=ydat) for dat in radar_data]
    else:
        figs = [plot_radargram(dat, xdat=xdat, ydat=ydat, x_range=None) for dat in radar_data]

    if s:
        [f[0].savefig(os.path.splitext(fn0)[0] + '.' + ftype, dpi=dpi) for f, fn0 in zip(figs, fn)]
    else:
        plt.show()


def plot_radargram(dat, xdat='tnum', ydat='twtt', x_range=(0, -1), cmap=plt.cm.gray_r, fig=None, ax=None, return_plotinfo=False):
    """This is the function to plot the normal radargrams that we are used to.

    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
    xdat: str, optional


    """
    if xdat not in ['tnum', 'dist']:
        raise ValueError('x axis choices are tnum or dist')
    if ydat not in ['twtt', 'depth']:
        raise ValueError('y axis choices are twtt or depth')

    if x_range is None:
        x_range = (0, -1)
    if x_range[-1] == -1:
        x_range = (x_range[0], dat.tnum)

    lims = np.percentile(dat.data[:, x_range[0]:x_range[-1]][~np.isnan(dat.data[:, x_range[0]:x_range[-1]])], (10, 90))
    clims = [lims[0] * 2 if lims[0] < 0 else lims[0] / 2, lims[1] * 2]

    if fig is not None:
        pass
    else:
        fig, ax = plt.subplots(figsize=(12, 8))

    if hasattr(dat.flags, 'elev') and dat.flags.elev:
        yd = dat.elevation
        ax.set_ylabel('Elevation (m)')
    else:
        ax.invert_yaxis()
        if ydat == 'twtt':
            yd = dat.travel_time
            ax.set_ylabel('Two way travel time (usec)')
        elif ydat == 'depth':
            if dat.nmo_depth is not None:
                yd = dat.nmo_depth
            else:
                yd = dat.travel_time / 2.0 * 1.69e8 * 1.0e-6
            ax.set_ylabel('Depth (m)')

    if xdat == 'tnum':
        xd = np.arange(int(dat.tnum))[x_range[0]:x_range[-1]]
        ax.set_xlabel('Trace number')
    elif xdat == 'dist':
        xd = dat.dist[x_range[0]:x_range[-1]]
        ax.set_xlabel('Distance (km)')

    if hasattr(dat.flags, 'elev') and dat.flags.elev:
        im = ax.imshow(dat.data[:, x_range[0]:x_range[-1]], cmap=cmap, vmin=lims[0], vmax=lims[1], extent=[np.min(xd), np.max(xd), np.min(yd), np.max(yd)], aspect='auto')
    else:
        im = ax.imshow(dat.data[:, x_range[0]:x_range[-1]], cmap=cmap, vmin=lims[0], vmax=lims[1], extent=[np.min(xd), np.max(xd), np.max(yd), np.min(yd)], aspect='auto')
    if not return_plotinfo:
        return fig, ax
    else:
        return im, xd, yd, x_range, lims


def plot_traces(dat, tr, ydat='twtt'):
    #Two options of trace input, a single trace or multiple
    if type(tr) == int:
        tr = (tr, tr + 1)
    if tr[0] == tr[1]:
        tr = (tr[0], tr[0] + 1)

    if ydat not in ['twtt', 'depth']:
        raise ValueError('y axis choices are twtt or depth')
    fig, ax = plt.subplots(figsize=(8, 12))
    lims = np.percentile(dat.data[:, tr[0]:tr[1]], (1, 99))
    ax.invert_yaxis()

    if ydat == 'twtt':
        yd = dat.travel_time
        ax.set_ylabel('Two way travel time (usec)')
    elif ydat == 'depth':
        if dat.nmo_depth is None:
            yd = dat.travel_time / 2.0 * 1.69e8 * 1.0e-6
        else:
            yd = dat.nmo_depth
        ax.set_ylabel('Depth (m)')

    for j in range(*tr):
        ax.plot(dat.data[:, j], yd)

    if lims[0] < 0 and lims[1] > 0:
        ax.set_xlim(lims[0], -lims[0])
    else:
        ax.set_xlim(*lims)
    ax.set_xlabel('Power')
    return fig, ax
