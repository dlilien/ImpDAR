#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""Plotting functions for radar data."""
import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from .load import load
from matplotlib.colors import is_color_like

# define a set of non-gray colors (from Paul Tol)
COLORS_NONGRAY = ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE',
                  '#882255', '#44AA99', '#999933', '#AA4499']


def plot(fns, tr=None, s=False, ftype='png', dpi=300, xd=False, yd=False,
         dualy=False, x_range=(0, -1), power=None, spectra=None,
         freq_limit=None, window=None, scaling='spectrum', filetype='mat',
         pick_colors=None, ft=False, hft=False, clims=None, cmap=plt.cm.gray,
         flatten_layer=None, *args, **kwargs):
    """Wrap a number of plot types.

    This should really only be used by the exectuables.
    If you are plotting yourself, just use the individual plotting
    functions that are described below.

    Parameters
    ----------
    fns: list of strs
        A list of filenames to plot individually.
    tr: tuple or int, optional
        Plot traces tr[1] to tr[2] (or trace tr) rather than the radargram.
        Default is None (plot radargram).
    power: int, optional
        If not None, then plot power returned from this layer
    filetype: str, optional
        Type of input file. Default mat.
    x_range: tuple, optional
        The range of traces to plot in the radargram.
        Default is (0, -1) (plot all traces)
    flatten_layer: int, optional
        Distort the radargram so this layer is flat. Default is None (do not distort).
    """
    radar_data = load(filetype, fns)

    if xd:
        xdat = 'dist'
    else:
        xdat = 'tnum'
    if yd:
        if dualy:
            raise ValueError('Only one of yd and dualy can be true')
        ydat = 'depth'
    elif dualy:
        ydat = 'dual'
    else:
        ydat = 'twtt'

    if (tr is not None) and (power is not None):
        raise ValueError('Cannot do both tr and power. Pick one')

    if tr is not None:
        figs = [plot_traces(dat, tr, ydat=ydat) for dat in radar_data]
    elif power is not None:
        # Do it all on one axis if power
        figs = [plot_power(radar_data, power)]
    elif ft:
        figs = [plot_ft(dat) for dat in radar_data]
    elif hft:
        figs = [plot_hft(dat) for dat in radar_data]
    elif spectra:
        figs = [plot_spectrogram(dat, spectra, window=window, scaling=scaling) for dat in radar_data]
    else:
        figs = [plot_radargram(dat,
                               xdat=xdat,
                               ydat=ydat,
                               x_range=None,
                               pick_colors=pick_colors,
                               clims=clims,
                               cmap=cmap,
                               flatten_layer=flatten_layer)
                for dat in radar_data]

    for fig, dat in zip(figs, radar_data):
        if dat.fn is not None:
            fig[0].canvas.set_window_title(dat.fn)

    if s:
        [f[0].savefig(os.path.splitext(fn0)[0] + '.' + ftype, dpi=dpi)
         for f, fn0 in zip(figs, fns)]
    else:
        plt.tight_layout()
        plt.show()


def plot_radargram(dat, xdat='tnum', ydat='twtt', x_range=(0, -1),
                   y_range=(0, -1), cmap=plt.cm.gray, fig=None, ax=None,
                   return_plotinfo=False, pick_colors=None, clims=None,
                   data_name='data', flatten_layer=None, middle_picks_only=False):
    """Plot a radio echogram.

    This function is a little weird since I want to be able to plot on top of
    existing figures/axes or on new figures an axes. There is therefore an
    argument `return_plotinfo` that funnels between these options and changes
    the return types.

    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
    xdat: str, optional
        The horizontal axis units. Either tnum or dist(ance).
    ydat: str, optional
        The vertical axis units. Either twtt or or depth. Default twtt.
    x_range: 2-tuple, optional
        The range of values to plot, in tnum space.
        Default is plot everything (0, -1)
    y_range: 2-tuple, optional
        The range of values to plot, in snum space.
        Default is plot everything (0, -1)
    cmap: matplotlib.pyplot.cm, optional
        The colormap to use
    fig: matplotlib.pyplot.Figure
        Figure canvas that should be plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that should be plotted upon
    data_name: str, optional
        The name of the data attribute. Default 'data'. Must exist.
    flatten_layer: int, optional
        Distort so this layer is flat
    middle_picks_only: bool, optional
        Allows you to specify color triples for plotting picks and not have them misinterptreted.


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
            The limits of the x range,
            after modification to remove negative indices
        clims: 2-tuple
            The limits of the colorbar
    """
    plotting_data = getattr(dat, data_name)
    if xdat not in ['tnum', 'dist']:
        raise ValueError('x axis choices are tnum or dist')
    elif (xdat == 'dist') and dat.dist is None:
        raise ValueError('xdat cannot be dist when the data has no dist')

    if x_range is None:
        x_range = (0, -1)
    if x_range[-1] == -1:
        x_range = (x_range[0], dat.tnum)

    if y_range is None:
        y_range = (0, -1)

    if y_range[-1] == -1:
        y_range = (y_range[0], dat.data.shape[0])

    if dat.data.dtype in [np.complex128]:
        def norm(x):
            return 10.0 * np.log10(np.absolute(x))
    else:
        def norm(x):
            return x

    if clims is None:
        clims = np.percentile(norm(plotting_data[y_range[0]:y_range[-1], x_range[0]:x_range[-1]][~np.isnan(dat.data[y_range[0]:y_range[-1], x_range[0]:x_range[-1]])]), (10, 90))

    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(12, 8))
    if ydat == 'elev':
        if hasattr(dat.flags, 'elev') and dat.flags.elev:
            yd = dat.elevation
            ax.set_ylabel('Elevation (m)')
        else:
            raise ValueError('Elevation plot requested but we have none')
    else:
        ax.invert_yaxis()
        if ydat == 'twtt':
            # we have a chance that there are NaNs after NMO correction...
            y_range = (max(y_range[0], np.min(np.where(~np.isnan(dat.travel_time))[0])), y_range[1])
            yd = dat.travel_time
            ax.set_ylabel('Two way travel time (usec)')
        elif ydat == 'depth':
            if dat.nmo_depth is not None:
                yd = dat.nmo_depth
            else:
                yd = dat.travel_time / 2.0 * (
                    1.69e8 * 1.0e-6)
            ax.set_ylabel('Depth (m)')
        elif ydat == 'dual':
            # we have a chance that there are NaNs after NMO correction...
            y_range = (max(y_range[0], np.min(np.where(~np.isnan(dat.travel_time))[0])), y_range[1])

            yd = dat.travel_time
            ax.set_ylabel('Two way travel time (usec)')
            ax2 = ax.twinx()
            if dat.nmo_depth is not None:
                yd2 = dat.nmo_depth
            else:
                yd2 = dat.travel_time / 2.0 * (
                    1.69e8 * 1.0e-6)
            ax2.set_ylabel('Approximate depth (m)')
            ax2.set_ylim(yd2[y_range[-1] - 1], yd2[y_range[0]])
        else:
            raise ValueError('Unrecognized ydat, choices are elev, twtt, \
                             depth, or dual')

    if xdat == 'tnum':
        xd = np.arange(int(dat.tnum))
        ax.set_xlabel('Trace number')
    elif xdat == 'dist':
        xd = dat.dist
        ax.set_xlabel('Distance (km)')

    if flatten_layer is not None:
        offset, _ = get_offset(dat, flatten_layer)

        # Now construct the data matrix
        tmp_data = np.zeros_like(dat.data)
        tmp_data[:, :] = np.nan
        for j in range(tmp_data.shape[1]):
            if np.isnan(offset[j]):
                continue
            if int(offset[j]) == 0:
                tmp_data[:, j] = dat.data[:, j]
            elif offset[j] < 0 and (abs(offset[j]) < dat.snum):
                tmp_data[:int(offset[j]), j] = dat.data[-int(offset[j]):, j]
            elif (abs(offset[j]) < dat.snum) and offset[j]:
                tmp_data[int(offset[j]):, j] = dat.data[:-int(offset[j]), j]
        im = ax.imshow(norm(tmp_data[:, x_range[0]:x_range[-1]]),
                       cmap=cmap,
                       vmin=clims[0],
                       vmax=clims[1],
                       extent=[np.min(xd[x_range[0]:x_range[-1]]), np.max(xd[x_range[0]:x_range[-1]]), np.max(yd[y_range[0]:y_range[-1]]), np.min(yd[y_range[0]:y_range[-1]])],
                       aspect='auto')
    elif hasattr(dat.flags, 'elev') and dat.flags.elev:
        im = ax.imshow(norm(dat.data[y_range[0]:y_range[-1],
                                     x_range[0]:x_range[-1]]),
                       cmap=cmap,
                       vmin=clims[0],
                       vmax=clims[1],
                       extent=[np.min(xd[x_range[0]:x_range[-1]]), np.max(xd[x_range[0]:x_range[-1]]), np.min(yd[y_range[0]:y_range[-1]]), np.max(yd[y_range[0]:y_range[-1]])],
                       aspect='auto')
    else:
        im = ax.imshow(norm(dat.data[y_range[0]:y_range[-1],
                                     x_range[0]:x_range[-1]]),
                       cmap=cmap,
                       vmin=clims[0],
                       vmax=clims[1],
                       extent=[np.min(xd[x_range[0]:x_range[-1]]), np.max(xd[x_range[0]:x_range[-1]]), np.max(yd[y_range[0]:y_range[-1]]), np.min(yd[y_range[0]:y_range[-1]])],
                       aspect='auto')

    if (pick_colors is not None) and pick_colors:
        plot_picks(dat, xd, yd, fig=fig, ax=ax, colors=pick_colors, flatten_layer=flatten_layer, just_middle=middle_picks_only, x_range=x_range)
    if not return_plotinfo:
        return fig, ax
    else:
        return im, xd, yd, x_range, clims


def plot_ft(dat, fig=None, ax=None, **line_kwargs):
    """Plot the Fourier spectrum of the data in the vertical.

    This will give the power spectral density in terms of the
    frequency (in MHz). We first fft, then average the fft.

    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
    fig: matplotlib.pyplot.Figure
        Figure canvas that should be plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that should be plotted upon
    **line_kwargs
        Arguments passed to the plotting call (e.g. color, linewidth)

    Returns
    -------
    fig: matplotlib.pyplot.Figure
        Figure canvas that was plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that were plotted upon
    """
    fft = np.fft.fft(dat.data, axis=0)
    fft_dat = np.mean(np.abs(fft) ** 2.0, axis=1)
    freq = np.fft.fftfreq(dat.snum) / dat.dt
    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(freq[freq >= 0] / 1.0e6, fft_dat[freq >= 0], **line_kwargs)
    ax.set_xlabel('Freq (MHz)')
    ax.set_ylabel('Power spectral density')
    return fig, ax


def plot_hft(dat, fig=None, ax=None):
    """Plot the Fourier spectrum of the data in the horizontal.

    This will give the power spectral density as a function of the
    horizontal wavelength (in meters). We first fft, then average the fft

    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
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
    fft = np.fft.fft(dat.data, axis=1)
    fft_dat = np.mean(np.abs(fft) ** 2.0, axis=0)

    # approximate as with the hbp
    freq = np.fft.fftfreq(dat.tnum)

    # we expect a divide by zero here
    with np.errstate(divide='ignore', invalid='ignore'):
        wavelength = dat.flags.interp[1] / freq
        wavelength[freq == 0.0] = np.inf

    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(wavelength[freq >= 0], fft_dat[freq >= 0])
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Power spectral density')
    return fig, ax


def plot_traces(dat, tr, ydat='twtt', fig=None, ax=None, linewidth=1.0,
                linestyle='solid'):
    """Plot power vs depth or twtt in a trace.

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
    # Two options of trace input, a single trace or multiple
    if hasattr(tr, '__iter__'):
        if not len(tr) == 2:
            raise ValueError('tr must either be a 2-tuple of bounds for the \
                             traces or a single trace index')
    if type(tr) == int:
        tr = (tr, tr + 1)
    elif tr[0] == tr[1]:
        tr = (tr[0], tr[0] + 1)

    if ydat not in ['twtt', 'depth', 'dual']:
        raise ValueError('y axis choices are twtt or depth')
    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(8, 12))
    # ax.set_xscale('symlog')
    lims = np.percentile(dat.data[:, tr[0]:tr[1]], (1, 99))
    if lims[0] == lims[1]:
        lims[1] = lims[0] + 1.
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
    elif ydat == 'dual':
        # we have a chance that there are NaNs after NMO correction...
        yd = dat.travel_time
        ax.set_ylabel('Two way travel time (usec)')
        ax2 = ax.twinx()
        if dat.nmo_depth is not None:
            yd2 = dat.nmo_depth
        else:
            yd2 = dat.travel_time / 2.0 * (1.69e8 * 1.0e-6)
        ax2.set_ylabel('Approximate depth (m)')
        ax2.set_ylim(yd2[-1], yd2[0])
    else:
        raise ValueError('Unrecognized y scale')

    for j in range(*tr):
        ax.plot(dat.data[:, j], yd, linewidth=linewidth, linestyle=linestyle)
    if lims[0] < 0 and lims[1] > 0:
        ax.set_xlim(lims[0], -lims[0])
    else:
        ax.set_xlim(*lims)
    ax.set_xlabel('Amplitude')
    return fig, ax


def plot_power(dats, idx, fig=None, ax=None, clims=None):
    """Make a plot of the reflected power along a given pick.

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
    # check to see if user entered an integer pick number
    try:
        idx = int(idx)
    except TypeError:
        raise TypeError('Please enter an integer pick number')

    if type(dats) not in [list, tuple]:
        dats = [dats]

    for dat in dats:
        if (dat.picks is None) or (dat.picks.picknums is None):
            raise ValueError('There are no picks on this radardata, \
                             cannot plot return power')

        if idx not in dat.picks.picknums:
            raise ValueError('Pick number {:d} not found in your file'.format(
                idx))

    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(8, 12))

    # Attempt to plot in projected coordinates
    if (dats[0].x_coord is not None) and (dats[0].y_coord is not None):
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

    pick_power = np.hstack([dat.picks.power[dat.picks.picknums.index(idx)
                                            ].flatten() for dat in dats])

    c = 10 * np.log10(pick_power)

    if clims is None:
        clims = np.percentile(c[~np.isnan(c)], (1, 99))

        # I think we throw an error if vmin=vmax
        # but we still want a plot of constant power
        if (clims[0] - clims[1]) / clims[0] < 1.0e-8:
            clims[0] = 0.99 * clims[0]
            clims[1] = 1.01 * clims[1]

    img = ax.scatter(lons.flatten(),
                     lats.flatten(),
                     c=c.flatten(),
                     vmin=clims[0],
                     vmax=clims[1])
    h = fig.colorbar(img)
    h.set_label('dB')
    ax.set_ylabel('Northing')
    ax.set_xlabel('Easting')

    return fig, ax


def plot_picks(rd, xd, yd, colors=None, flatten_layer=None, fig=None, ax=None, just_middle=False, picknums=None, x_range=None, **plotting_kwargs):
    """Update the plotting of the current pick.

    Parameters
    ----------
    colors: str
        You have choices here. This can be a npicksx3 list, an npicks list of
        3-letter strings, a 3 letter string, a single string, or a npicks list.
        Any of the x3 options are interpretted as top, middle, bottom colors.
        If it is a string, the lines are all plotted in this color. If it is
        a list, the different values are used for the different lines.
    flatten_layer: int, optional
        Make this layer flat in the plot. Distorts all layers. Default is no
        distortion.
    """
    if x_range is None:
        x_range = (0, -1)
    if x_range[-1] == -1:
        x_range = (x_range[0], rd.tnum)

    if ax is None:
        if fig is not None:
            ax = plt.gca()
        else:
            fig, ax = plt.subplots()

    # just do nothing if we have no picks
    if rd.picks is None or rd.picks.samp1 is None:
        return fig, ax

    offset, mask = get_offset(rd, flatten_layer)

    if picknums is None:
        if rd.picks.picknums is None:
            return fig, ax
        picknums = rd.picks.picknums

    variable_colors = False
    if not colors:  # may be False or None
        cl = 'mgm'
    else:
        if type(colors) == str:
            if len(colors) == 3:
                cl = colors
            else:
                cl = ('none', colors, 'none')
        elif (type(colors) == bool) and colors:
            colors = (COLORS_NONGRAY * (rd.picks.samp1.shape[0] // len(COLORS_NONGRAY) + 1))[:len(picknums)]
            variable_colors = True
        elif not len(colors) == len(picknums):
            if (len(colors) == 3) and not just_middle:
                cl = colors
            else:
                raise ValueError('If not a string, must have length 3 or length npicks')
        else:
            variable_colors = True

    for j, pn in enumerate(picknums):
        # use i and j so that we can color out of order
        i = rd.picks.picknums.index(pn)
        if variable_colors:
            if hasattr(colors[j], '__len__') and len(colors[j]) == 3 and not just_middle:
                cl = colors[j]
            elif is_color_like(colors[j]):
                cl = ('none', colors[j], 'none')
            else:
                raise ValueError('Color ', colors[j], ' not defined')
        c = np.zeros(xd.shape)
        c[:] = np.nan
        comb_mask = np.logical_or(mask, np.isnan(rd.picks.samp2[i, :]))
        c[~comb_mask] = yd[(rd.picks.samp2[i, :] + offset)[~comb_mask].astype(int)]
        t = np.zeros(xd.shape)
        t[:] = np.nan
        comb_mask = np.logical_or(mask, np.isnan(rd.picks.samp1[i, :]))
        t[~comb_mask] = yd[(rd.picks.samp1[i, :] + offset)[~comb_mask].astype(int)]
        b = np.zeros(xd.shape)
        b[:] = np.nan
        comb_mask = np.logical_or(mask, np.isnan(rd.picks.samp3[i, :]))
        b[~comb_mask] = yd[(rd.picks.samp3[i, :] + offset)[~comb_mask].astype(int)]
        ax.plot(xd[x_range[0]:x_range[1]], c[x_range[0]:x_range[1]], color=cl[1], **plotting_kwargs)
        ax.plot(xd[x_range[0]:x_range[1]], t[x_range[0]:x_range[1]], color=cl[0], **plotting_kwargs)
        ax.plot(xd[x_range[0]:x_range[1]], b[x_range[0]:x_range[1]], color=cl[2], **plotting_kwargs)
    return fig, ax


def plot_spectrogram(dat, freq_limit=None, window=None,
                     scaling='spectrum', fig=None, ax=None, **kwargs):
    """Make a plot of power spectral density across all traces of a radar profile.

    Parameters
    ----------
    dat: impdar.lib.RadarData.Radardata
        The RadarData object to plot.
    freq_limit: tuple
        The minimum and maximum frequency (in MHz) to limit the y-axis to
    window: str, optional
        Type of window to be used for the signal.periodogram() method.

        Default hamming.
        `Further information <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.periodogram.html#scipy.signal.periodogram>`_
    scaling: str, optional
        Whether to plot power spectral density or power spectrum
        'density' or 'spectrum', the default being 'spectrum'.
        `Further information <https://docs.scipy.org/doc/scipy-0.14.0/referenc\
        e/generated/scipy.signal.periodogram.html#scipy.signal.periodogram>`_
    fig: matplotlib.pyplot.Figure, optional
        Figure canvas that should be plotted upon
    ax: matplotlib.pyplot.Axes, optional
        Axes that should be plotted upon

    Returns
    -------
    fig: matplotlib.pyplot.Figure
        Figure canvas that was plotted upon
    ax: matplotlib.pyplot.Axes
        Axes that were plotted upon
    """
    # get the timestep variable, remove singleton dimension
    timestep = dat.dt

    # calculate frequency information from timestep variable
    fs = 1. / timestep

    # shape of data profile should be (samples, traces)
    shape = np.shape(dat.data)
    # extract the number of traces
    traces = shape[1]

    # iterate through traces, record frequencies and powers
    powers = []
    for trace in range(traces):
        # get frequency and power information from trace
        # hanning window will filter out certain frequencies,
        # so it is optional to use it or not
        freq, power = signal.periodogram(dat.data[:, trace],
                                         fs=fs,
                                         window=window,
                                         scaling=scaling)
        powers.append(power)

    # extract trace number from matlab file
    x = dat.trace_num

    # frequency range will be the same, so we can select the first element
    # set frequency range to be in MHz
    y = freq / 1.0e6
    xx, yy = np.meshgrid(x, y)

    # set figure and axis if they are not None
    if fig is not None:
        if ax is None:
            ax = plt.gca()
    else:
        fig, ax = plt.subplots(figsize=(10, 7))

    # plot in MHz
    contours = ax.contourf(xx, yy, np.transpose(powers))

    # set colorbar and colorbar label
    cbarlabel = 'Power (Amplitude **2)'
    cbar = plt.colorbar(contours,
                        shrink=0.9,
                        orientation='vertical',
                        pad=0.03,
                        ax=ax)
    cbar.set_label(cbarlabel)

    # check to make sure freq_limit is not <= smallest freq
    if freq_limit is not None:
        if hasattr(freq_limit, '__len__'):
            if freq_limit[1] < np.nanmin(y):
                raise ValueError('Y-axis limit {} MHz too low.'.format(freq_limit[1]))
            if freq_limit[1] > np.nanmax(y):
                print('Warning: y-axis limit large compared to the frequencies plotted')

            # limit y-axis to freq_limit, maximum power output
            # else, no need to do anything
            ax.set_ylim(freq_limit[0], freq_limit[1])
        else:
            print('Frequency limit should be a tuple of low, high. Ignoring.')

    # add x and y labels
    ax.set_xlabel('Trace Number')
    ax.set_ylabel('Frequency (MHz)')

    # set title
    title = 'PSD(tnum, f)'
    ax.set_title(title)

    return fig, ax


def get_offset(dat, flatten_layer=None):
    if flatten_layer is None:
        offset = np.zeros((dat.data.shape[1]))
        mask = np.zeros((dat.tnum, ), dtype=bool)
    else:
        if flatten_layer not in dat.picks.picknums:
            raise ValueError('That layer is not in existence, cannot flatten')
        layer_ind = dat.picks.picknums.index(flatten_layer)
        layer_depth = dat.picks.samp2[layer_ind, :]
        zero_offset = int(np.nanmean(layer_depth))
        offset = zero_offset - layer_depth
        mask = np.isnan(dat.picks.samp2[layer_ind, :])
    return offset, mask
