#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
from .load import load


def plot(fns, tr=None, s=False, ftype='png', dpi=300, xd=False, yd=False, x_range=(0, -1), power=None, spectra=False, freq_limit=None, window=None, scale='spectrum', gssi=False, pe=False, gprMax=False, gecko=False, segy=False, *args, **kwargs):
    """We have an overarching function here to handle a number of plot types

    Parameters
    ----------
    fns: list of strs
        A list of filenames to plot individually.
    tr: tuple or int, optional
        Plot traces tr[1] to tr[2] (or trace tr) rather than the radargram. Default is None (plot radargram)
    power: int, optional
        If not None, then plot power returned from this layer
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
        radar_data = load('gssi', fns)
    elif pe:
        radar_data = load('pe', fns)
    elif gecko:
        radar_data = load('gecko', fns)
    elif gprMax:
        radar_data = load('gprMax', fns)
    elif segy:
        radar_data = load('segy', fns)
    else:
        radar_data = load('mat', fns)

    if xd:
        xdat = 'dist'
    else:
        xdat = 'tnum'
    if yd:
        ydat = 'depth'
    else:
        ydat = 'twtt'

    if (tr is not None) and (power is not None):
        raise ValueError('Cannot do both tr and power. Pick one')

    if tr is not None:
        figs = [plot_traces(dat, tr, ydat=ydat) for dat in radar_data]
    elif power is not None:
        # Do it all on one axis if power
        figs = [plot_power(radar_data, power)]
    elif spectra != False:
        #call specdense() here
        figs = [specdense(radar_data, freq_limit, window, scale)]
    else:
        figs = [plot_radargram(dat, xdat=xdat, ydat=ydat, x_range=None) for dat in radar_data]

    for fig, dat in zip(figs, radar_data):
        if dat.fn is not None:
            fig[0].canvas.set_window_title(dat.fn)

    if s:
        [f[0].savefig(os.path.splitext(fn0)[0] + '.' + ftype, dpi=dpi) for f, fn0 in zip(figs, fns)]
    else:
        plt.tight_layout()
        plt.show()


def plot_radargram(dat, xdat='tnum', ydat='twtt', x_range=(0, -1), cmap=plt.cm.gray, fig=None, ax=None, return_plotinfo=False, pick_colors=None):
    """This is the function to plot the normal radargrams that we are used to.

    This function is a little weird since I want to be able to plot on top of existing figures/axes or on new figures an axes. There is therefore an argument `return_plotinfo` that funnels between these options and changes the return types

    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
    xdat: str, optional
        The horizontal axis units. Either tnum or distance.
    ydat: str, optional
        The vertical axis units. Either twtt or or depth. Default twtt.
    x_range: 2-tuple, optional
        The range of values to plot. Default is plot everything (0, -1)
    cmap: matplotlib.pyplot.cm, optional
        The colormap to use
    fig: matplotlib.pyplot.Figure
        Figure canvas that should be plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that should be plotted upon


    Returns
    -------
    If not return_plotinfo

        fig: matplotlib.pyplot.Figure
            Figure canvas that was plotted upon
        ax: matplotlib.pyplot.Axes
            Axes that were plotted upon

    else
        im: pyplot.imshow
            The image object plotted
        xd: np.ndarray
            The x values of the plot
        yd: np.ndarray
            The y values of the plot
        x_range: 2-tuple
            The limits of the x range, after modification to remove negative indices
        clims: 2-tuple
            The limits of the colorbar

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

    if fig is not None:
        if ax is None:
            ax = plt.gca()
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

    if pick_colors is not None:
        plot_picks(dat, xd, yd, fig=fig, ax=ax, colors=pick_colors)
    if not return_plotinfo:
        return fig, ax
    else:
        return im, xd, yd, x_range, lims


def plot_traces(dat, tr, ydat='twtt', fig=None, ax=None):
    """Plot power vs depth or twtt in a trace

    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
    tr: int or 2-tuple
        Either a single trace or a range of traces to plot
    ydat: str, optional
        The vertical axis units. Either twtt or or depth. Default twtt.
    fig: matplotlib.pyplot.Figure
        Figure canvas that should be plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that should be plotted upon

    Returns
    -------
    fig: matplotlib.pyplot.Figure
        Figure canvas that was plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that were plotted upon
    """
    #Two options of trace input, a single trace or multiple
    if hasattr(tr, '__iter__'):
        if not len(tr) == 2:
            raise ValueError('tr must either be a 2-tuple of bounds for the traces or a single trace index')
    if type(tr) == int:
        tr = (tr, tr + 1)
    elif tr[0] == tr[1]:
        tr = (tr[0], tr[0] + 1)

    if ydat not in ['twtt', 'depth']:
        raise ValueError('y axis choices are twtt or depth')
    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(8, 12))
    ax.set_xscale('symlog')
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


def plot_power(dats, idx, fig=None, ax=None):
    """Make a plot of the reflected power along a given pick


    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
    idx: int
        A picknum in the dat.picks.picknum array
    fig: matplotlib.pyplot.Figure
        Figure canvas that should be plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that should be plotted upon

    Returns
    -------
    fig: matplotlib.pyplot.Figure
        Figure canvas that was plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that were plotted upon
    """
    #check to see if user entered an integer pick number
    try:
        idx = int(idx)
    except TypeError:
        raise TypeError('Please enter an integer pick number')

    if type(dats) not in [list, tuple]:
        dats = [dats]

    for dat in dats:
        if (dat.picks is None) or (dat.picks.picknums is None):
            raise ValueError('There are no picks on this radardata, cannot plot return power')

        if idx not in dat.picks.picknums:
            raise ValueError('Pick number {:d} not found in your file'.format(idx))

    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(8, 12))

    # Attempt to plot in projected coordinates
    if dats[0].x_coord is not None:
        if len(dats) > 1:
            lons = np.hstack([dat.x_coord for dat in dats])
            lats = np.hstack([dat.y_coord for dat in dats])
        else:
            lons = dats[0].x_coord
            lats = dats[0].y_coord
    else:
        if len(dats) > 1:
            lons = np.hstack([dat.long for dat in dats])
            lats = np.hstack([dat.lat for dat in dats])
        else:
            lons = dats[0].long
            lats = dats[0].lat

    pick_power = np.hstack([dat.picks.power[dat.picks.picknums.index(idx)].flatten() for dat in dats])

    c = 10 * np.log10(pick_power)
    clims = np.percentile(c[~np.isnan(c)], (1, 99))

    # I think we throw an error if vmin=vmax, but we still want a plot of constant power
    if (clims[0] - clims[1]) / clims[0] < 1.0e-8:
        clims[0] = 0.99 * clims[0]
        clims[1] = 1.01 * clims[1]

    img = ax.scatter(lons.flatten(), lats.flatten(), c=c.flatten(), vmin=clims[0], vmax=clims[1])
    h = fig.colorbar(img)
    h.set_label('dB')
    ax.set_ylabel('Northing')
    ax.set_xlabel('Easting')

    return fig, ax


def plot_picks(rd, xd, yd, colors=None, fig=None, ax=None):
    """Update the plotting of the current pick.

    Parameters
    ----------
    colors: str
        You have choices here. This can be a npicksx3 list, an npicks list of 3-letter strings, a 3 letter string, a single string, or a npicks list. Any of the x3 options are interpretted as top, middle, bottom colors. The others are
    picker:
        argument to pass to plot of cline (if new) for selection tolerance (use if plotting in select mode)
    """

    if ax is None:
        if fig is not None:
            ax = plt.gca()
        else:
            fig, ax = plt.subplots()

    # just do nothing if we have no picks
    if rd.picks is None or rd.picks.samp1 is None:
        return fig, ax

    variable_colors = False
    if colors is None:
        cl = 'mgm'
    else:
        if type(colors) == str:
            if len(colors) == 3:
                cl = colors
            else:
                cl = ('none', colors, 'none')
        elif not len(colors) == rd.picks.samp1.shape[0]:
            raise ValueError('If not a string, must have same length as the picks')
        else:
            variable_colors = True

    for i in range(rd.picks.samp1.shape[0]):
        if variable_colors:
            if len(colors[i]) == 3:
                cl = colors[i]
            else:
                cl = ('none', colors[i], 'none')
        c = np.zeros(xd.shape)
        c[:] = np.nan
        c[~np.isnan(rd.picks.samp2[i, :])] = yd[rd.picks.samp2[i, :][~np.isnan(rd.picks.samp2[i, :])].astype(int)]
        t = np.zeros(xd.shape)
        t[:] = np.nan
        t[~np.isnan(rd.picks.samp1[i, :])] = yd[rd.picks.samp1[i, :][~np.isnan(rd.picks.samp1[i, :])].astype(int)]
        b = np.zeros(xd.shape)
        b[:] = np.nan
        b[~np.isnan(rd.picks.samp3[i, :])] = yd[rd.picks.samp3[i, :][~np.isnan(rd.picks.samp3[i, :])].astype(int)]
        ax.plot(xd, c, color=cl[1])
        ax.plot(xd, t, color=cl[0])
        ax.plot(xd, b, color=cl[2])
    return fig, ax




"""Make a plot of power spectral density across all traces of a radar profile.


    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
    freq_limit: float
        The maximum frequency (in MHz) to limit the y-axis to

    For further information on the 'window' and 'scale' parameters, please see:
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.periodogram.html#scipy.signal.periodogram
    window: string
        Type of window to be used for the signal.periodogram() method.
    scale:
        Whether to plot power spectral density or power spectrum
        'density' or 'spectrum', the default being 'spectrum'
    fig: matplotlib.pyplot.Figure
        Figure canvas that should be plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that should be plotted upon

    Returns
    -------
    fig: matplotlib.pyplot.Figure
        Figure canvas that was plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that were plotted upon
    """

def specdense(dat, freq_limit, window, scale, fig=None, ax=None, **kwargs):
    
    dat = dat[0]

    #get the timestep variable, remove singleton dimension
    timestep = dat.dt

    #calculate frequency information from timestep variable
    fs = 1/timestep

    #extract radar data from matlab file
    data = dat.data

    #shape of data profile should be (samples, traces)
    shape = np.shape(data)
    #extract the number of traces
    traces = shape[1]

    #iterate through traces, record frequencies and powers
    freqs, powers = [], []
    for trace in range(traces):
        #get frequency and power information from trace
        #hanning window will filter out certain frequencies, so it is optional to use it or not
        if window==None:
            f, p = signal.periodogram(data[:, trace], fs=fs, scaling=scale)
        else:
            f, p = signal.periodogram(data[:, trace], fs=fs, window=window, scaling=scale)
        freqs.append(f)
        powers.append(p)

    #extract trace number from matlab file
    x = dat.trace_num

    #frequency range will be the same, so we can select the first element
    #set frequency range to be in MHz
    y = freqs[0]/1e6
    xx, yy = np.meshgrid(x, y)

    #set figure and axis if they are not None
    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(10, 7))

    #plot in MHz
    p = ax.contourf(xx, yy, np.transpose(powers))

    #set colorbar and colorbar label
    cbarlabel = 'Power (Amplitude **2)'
    cbar = plt.colorbar(p, shrink=0.9, orientation='vertical', pad=0.03, ax=ax)
    cbar.set_label(cbarlabel)

    #check to make sure freq_limit is not <= 0 or more than the largest frequency
    if freq_limit is not None:
        if np.logical_or(freq_limit <= 0, freq_limit > np.max(y)):
            raise ValueError('Y-axis limit {} MHz not found in frequencies.'.format(freq_limit))
            return

        #limit y-axis to freq_limit, maximum power output
        #else, no need to do anything
        ax.set_ylim(0, freq_limit)

    #add x and y labels
    ax.set_xlabel('Trace Number')
    ax.set_ylabel('Frequency (MHz)')

    #set title
    title = 'Power Spectral Density as a Function of Trace Number and Frequency'

    #add space between the title and the plot
    ax.set_title(title, pad=20)

    return fig, ax
