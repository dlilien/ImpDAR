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
from matplotlib.widgets import Slider, RadioButtons, TextBox, Button
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
        xdat = 'tracenum'
    if yd:
        ydat = 'depth'
    else:
        ydat = 'twtt'

    if tr is not None:
        figs = [plot_traces(dat, tr) for dat in radar_data]
    else:
        figs = [plot_radargram(dat, interactive=not s, xdat=xdat, ydat=ydat, x_range=None) for dat in radar_data]

    if s:
        [f.savefig(os.path.splitext(fn0)[0] + '.' + ftype, dpi=dpi) for f, fn0 in zip(figs, fn)]
    else:
        plt.show()


def plot_radargram(dat, xdat='tracenum', ydat='twtt', interactive=False, x_range=(0, -1), cmap=plt.cm.gray_r):
    if xdat not in ['tracenum', 'dist']:
        raise ValueError('x axis choices are tracenum or dist')
    if ydat not in ['twtt', 'depth']:
        raise ValueError('y axis choices are twtt or depth')

    if x_range is None:
        x_range = (0, -1)

    lims = np.percentile(dat.data[:, x_range[0]:x_range[-1]][~np.isnan(dat.data[:, x_range[0]:x_range[-1]])], (10, 90))
    clims = [lims[0] * 2 if lims[0] < 0 else lims[0] / 2, lims[1] * 2]

    if not interactive:
        fig, ax = plt.subplots(figsize=(12, 8))
    else:
        fig, (ax, ax_slider1, ax_slider2) = plt.subplots(nrows=3, figsize=(12, 8), gridspec_kw={'height_ratios': (1, 0.02, 0.02)})

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

    if xdat == 'tracenum':
        xd = np.arange(int(dat.snum))[x_range[0]:x_range[-1]]
        ax.set_xlabel('Trace number')
    elif xdat == 'dist':
        xd = dat.dist[x_range[0]:x_range[-1]]
        ax.set_xlabel('Distance (km)')

    if hasattr(dat.flags, 'elev') and dat.flags.elev:
        im = ax.imshow(dat.data[:, x_range[0]:x_range[-1]], cmap=cmap, vmin=lims[0], vmax=lims[1], extent=[np.min(xd), np.max(xd), np.min(yd), np.max(yd)], aspect='auto')
    else:
        im = ax.imshow(dat.data[:, x_range[0]:x_range[-1]], cmap=cmap, vmin=lims[0], vmax=lims[1], extent=[np.min(xd), np.max(xd), np.max(yd), np.min(yd)], aspect='auto')

    if interactive:
        slidermin = Slider(ax_slider1, 'Min', clims[0], clims[1], valinit=lims[0])
        slidermax = Slider(ax_slider2, 'Max', clims[0], clims[1], valinit=lims[1], slidermin=slidermin)

        def update(val):
            im.set_clim(vmin=slidermin.val, vmax=slidermax.val)

        slidermin.on_changed(update)
        slidermax.on_changed(update)
    return fig


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
        yd = dat.nmo_depth
        ax.set_ylabel('Depth (m)')

    for j in range(*tr):
        ax.plot(dat.data[:, j], yd)

    if lims[0] < 0 and lims[1] > 0:
        ax.set_xlim(lims[0], -lims[0])
    else:
        ax.set_xlim(*lims)
    ax.set_xlabel('Power')
    return fig


class interactive_plot():
    
    def __init__(self, dat, xdat='tracenum', ydat='twtt', x_range=(0, -1)):
        plt.ion()
        # line is the matplotlib object of the current pick
        self.line = None
        # pick_pts contains our picked points, which will differ from what we want to save in the file. Need to have one per line.
        self.pick_pts = []
        self.dat = dat
        self.current_layer = 0

        if self.dat.picks is not None:
            self.pick_pts = [p[~np.isnan(p)].tolist() for p in self.dat.picks]

        try:
            if xdat not in ['tracenum', 'dist']:
                raise ValueError('x axis choices are tracenum or dist')
            if ydat not in ['twtt', 'depth']:
                raise ValueError('y axis choices are twtt or depth')

            if x_range is None:
                x_range = (0, -1)
            if x_range[-1] == -1:
                x_range = (x_range[0], self.dat.snum)

            self.lims = np.percentile(dat.data[:, x_range[0]:x_range[-1]], (10, 90))
            self.clims = [self.lims[0] * 2 if self.lims[0] < 0 else self.lims[0] / 2, self.lims[1] * 2]

            self.fig, self.ax = plt.subplots(figsize=(12, 8))

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

            if xdat == 'tracenum':
                self.xd = np.arange(int(self.dat.snum))[x_range[0]:x_range[-1]]
                self.ax.set_xlabel('Trace number')
            elif xdat == 'dist':
                self.xd = self.dat.dist[x_range[0]:x_range[-1]]
                self.ax.set_xlabel('Distance (km)')

            self.im = self.ax.imshow(dat.data[:, x_range[0]:x_range[-1]], cmap=plt.cm.gray_r, vmin=self.lims[0], vmax=self.lims[1], extent=[np.min(self.xd), np.max(self.xd), np.max(self.yd), np.min(self.yd)], aspect='auto')
            self.fig.canvas.mpl_connect('key_press_event', self.press)
            self.fig.canvas.mpl_connect('button_press_event', self.click)
            plt.show()
        except KeyboardInterrupt:
            plt.close('all')

    def press(self, event):
        if event.key == 'f':
            self.press_f()
        elif event.key == ' ':
            self.press_space()
        elif event.key == 'p':
            self.press_p()

    def press_f(self):
        gs = gridspec.GridSpec(5, 5)
        self.fig_fiddle = plt.figure(figsize=(6, 2))

        def press_f_fiddleplot(event):
            if event.key == 'f' or event.key == 'x':
                plt.close(self.fig_fiddle)

        self.fig_fiddle.canvas.mpl_connect('key_press_event', press_f_fiddleplot)

        ax_slider1 = plt.subplot(gs[0,:-2])
        ax_slider2 = plt.subplot(gs[1,:-2])
        ax_bwb = plt.subplot(gs[2:4, 0])
        ax_color = plt.subplot(gs[0:5, 4])
        ax_freq = plt.subplot(gs[4,:-2])

        slidermin = Slider(ax_slider1, 'Min', self.clims[0], self.clims[1], valinit=self.lims[0])
        slidermax = Slider(ax_slider2, 'Max', self.clims[0], self.clims[1], valinit=self.lims[1], slidermin=slidermin)
        slider_freq = Slider(ax_freq, 'Frequency', 0, 50, valinit=4, valstep=1, valfmt='%2.0f MHz')
        radio_bwb = RadioButtons(ax_bwb, ('BWB', 'WBW'))
        radio_color = RadioButtons(ax_color, ('gray', 'bwr', 'viridis', 'plasma', 'inferno', 'magma', 'seismic', 'hsv', 'jet', 'rainbow'))

        def lim_update(val):
            self.im.set_clim(vmin=slidermin.val, vmax=slidermax.val)
            self.lims[0] = slidermin.val
            self.lims[1] = slidermax.val
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        
        def color_update(val):
            self.im.set_cmap(getattr(plt.cm, val))
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()

        def bwb_update(label):
            pass

        def freq_update(val):
            pass

        try:
            slidermin.on_changed(lim_update)
            slidermax.on_changed(lim_update)
        except ValueError:
            pass

        slider_freq.on_changed(freq_update)
        radio_bwb.on_clicked(bwb_update)
        radio_color.on_clicked(color_update)

        plt.show()


    def press_p(self):
        gs = gridspec.GridSpec(2, 3)
        self.fig_select = plt.figure(figsize=(3, 1))

        def press_p_select(event):
            if event.key == 'p' or event.key == 'x':
                plt.close(self.fig_select)

        self.fig_select.canvas.mpl_connect('key_press_event', press_p_select)

        ax_select = plt.subplot(gs[0,0])
        ax_new = plt.subplot(gs[1,0])

        def change_layer(val):
            try:
                val = int(val)
            except:
                raise TypeError('Cannot parse non-int val')

            if val < 0 or val >= len(self.dat.picks):
                raise ValueError('Value not a layer')

        def next_and_exit():
            self.press_space()
            plt.close(self.fig_select)

        tb = TextBox(ax_select, initial=self.current_layer, label='Layer number')
        next = Button(ax_new, 'New')
        next.on_clicked(next_and_exit)
        tb.on_submit(change_layer)


    def press_space(self):
        if self.dat.picks is None:
            self.dat.picks = []
        d = np.zeros((self.dat.snum, ))
        d[d == 0] = np.nan
        self.dat.picks.append(d)
        self.pick_pts.append([])
        
        self.current_layer = len(self.dat.picks) - 1
        self.line = self.ax.plot(self.xd, self.dat.picks[self.current_layer])

    def click(self, event):
        if self.line is None:
            self.press_space()
        if event.button == 1:
            self.pick_pts[self.current_layer].append((event.xdata, event.ydata))
            # closest_coord = np.argmin(np.abs(self.xd - event.xdata))
            plot_pts = np.sort(np.array(self.pick_pts[self.current_layer], dtype=[('x', float), ('y', float)]), order='x', axis=0) 
            # self.dat.picks[self.current_layer][closest_coord] = event.ydata
            self.line[0].set_data(plot_pts['x'], plot_pts['y'])
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
