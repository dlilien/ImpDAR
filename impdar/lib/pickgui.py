#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import sys
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt5.QtWidgets import QFileDialog, QMessageBox
else:
    from matplotlib.backends.backend_qt4agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, RadioButtons, TextBox, Button

import os.path
import numpy as np

from . import plot
from .ui import RawPickGUI


def pick(radardata, guard_save=True):
    app = QtWidgets.QApplication(sys.argv)
    ip = InteractivePicker(radardata)
    ip.show()
    sys.exit(app.exec_())


class InteractivePicker(QtWidgets.QMainWindow, RawPickGUI.Ui_MainWindow):
    
    def __init__(self, dat, xdat='tracenum', ydat='twtt', x_range=(0, -1), guard_save=False):
        # Next line is required for Qt, then give us the layout
        super(self.__class__, self).__init__()
        self.setupUi(self)

        # Connect the menu to actions
        self.actionSave_as.triggered.connect(self._save_as)
        
        # Easy access to normal mpl figure and axes
        self.ax = self.FigCanvasWidget.canvas.ax
        self.fig = self.FigCanvasWidget.canvas.fig
        plt.ion()

        # Two constants to keep track of how to prompt for saves
        self.firstsave = True
        self.saved = True

        # Set defaults
        self.bwb = 'bwb'
        self.freq = 4

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


            if hasattr(dat.flags, 'elev') and dat.flags.elev:
                yd = dat.elevation
                self.ax.set_ylabel('Elevation (m)')
            else:
                self.ax.invert_yaxis()
                if ydat == 'twtt':
                    self.yd = dat.travel_time
                    self.ax.set_ylabel('Two way travel time (usec)')
                elif ydat == 'depth':
                    if dat.nmo_depth is not None:
                        self.yd = dat.nmo_depth
                    else:
                        self.yd = dat.travel_time / 2.0 * 1.69e8 * 1.0e-6
                    self.ax.set_ylabel('Depth (m)')

            if xdat == 'tracenum':
                self.xd = np.arange(int(self.dat.snum))[x_range[0]:x_range[-1]]
                self.ax.set_xlabel('Trace number')
            elif xdat == 'dist':
                self.xd = self.dat.dist[x_range[0]:x_range[-1]]
                self.ax.set_xlabel('Distance (km)')

            self.im = self.ax.imshow(dat.data[:, x_range[0]:x_range[-1]], cmap=plt.cm.gray_r, vmin=self.lims[0], vmax=self.lims[1], extent=[np.min(self.xd), np.max(self.xd), np.max(self.yd), np.min(self.yd)], aspect='auto')
            self.kpid = self.fig.canvas.mpl_connect('key_press_event', self.press)
            self.bpid = self.fig.canvas.mpl_connect('button_press_event', self.click)
            plt.show(self.fig)
        except KeyboardInterrupt:
            plt.close('all')

    def press(self, event):
        if event.key == 'd':
            self.press_f()
        elif event.key == ' ':
            self.press_space()
        elif event.key == 'n':
            self.press_p()

    def press_f(self):
        gs = gridspec.GridSpec(5, 5)
        self.fig_fiddle = plt.figure(figsize=(6, 2))

        def press_f_fiddleplot(event):
            if event.key == 'd' or event.key == 'x':
                plt.close(self.fig_fiddle)

        self.fig.canvas.mpl_disconnect(self.kpid)
        # self.fig.canvas.mpl_disconnect(self.bpid)
        self.kpid = self.fig_fiddle.canvas.mpl_connect('key_press_event', press_f_fiddleplot)

        ax_slider1 = plt.subplot(gs[0,:-2])
        ax_slider2 = plt.subplot(gs[1,:-2])
        ax_bwb = plt.subplot(gs[2:4, 0])
        ax_color = plt.subplot(gs[0:5, 4])
        ax_freq = plt.subplot(gs[4,:-2])

        self.slidermin = Slider(ax_slider1, 'Min', self.clims[0], self.clims[1], valinit=self.lims[0])
        self.slidermax = Slider(ax_slider2, 'Max', self.clims[0], self.clims[1], valinit=self.lims[1], slidermin=self.slidermin)
        self.slider_freq = Slider(ax_freq, 'Frequency', 0, 1000, valinit=4, valstep=1, valfmt='%2.0f MHz')
        self.radio_bwb = RadioButtons(ax_bwb, ('BWB', 'WBW'))
        self.radio_color = RadioButtons(ax_color, ('gray', 'bwr', 'viridis', 'plasma', 'inferno', 'magma', 'seismic', 'hsv', 'jet', 'rainbow'))

        def lim_update(val):
            self.im.set_clim(vmin=self.slidermin.val, vmax=self.slidermax.val)
            self.lims[0] = self.slidermin.val
            self.lims[1] = self.slidermax.val
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        
        def color_update(val):
            self.im.set_cmap(getattr(plt.cm, val))
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()

        def bwb_update(label):
            self.bwb = label

        def freq_update(val):
            self.freq = val

        try:
            self.slidermin.on_changed(lim_update)
            self.slidermax.on_changed(lim_update)
        except ValueError:
            pass

        self.slider_freq.on_changed(freq_update)
        self.radio_bwb.on_clicked(bwb_update)
        self.radio_color.on_clicked(color_update)

        self.fig_fiddle.canvas.mpl_disconnect(self.kpid)
        self.fig.canvas.mpl_connect('key_press_event', self.press)

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
        if self.FigCanvasWidget.mpl_toolbar._active is not None:
            return
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
            self.saved = False
    
    # We have a section here that makes sure you don't close
    def closeEvent(self, event):
        if not self.saved:
            self._save_cancel_close(event)

    def _save_cancel_close(self, event):
        dialog = QMessageBox()
        dialog.setStandardButtons(QMessageBox.Save | QMessageBox.Close | QMessageBox.Cancel)
        result = dialog.exec()
        if result == QMessageBox.Cancel:
            event.ignore()
        elif result == QMessageBox.Close:
            event.accept()
        else:
            if self.firstsave:
                event.accept() if self._save_as(event) else event.ignore()
            else:
                self._save()
                event.accept()

    def _save_inplace(self, evt):
        """Save the file without changing name"""
        self.dat.save(self.dat.fn)
        print('Saved as ' + self.dat.fn)
        self.saved = True
        plt.close(self.fig_save)

    def _save_pick(self, evt):
        """Save with _pick appended"""
        self.dat.save(self.dat.fn[:-4] + '_pick.mat')
        print('Saved as ' + self.dat.fn[:-4] + '_pick.mat')
        self.saved = True
        plt.close(self.fig_save)

    def _save_as(self, event=None):
        """Fancy file handler for gracious exit"""
        fn, test = QFileDialog.getSaveFileName(self, "QFileDialog.getSaveFileName()", self.dat.fn, "All Files (*);;mat Files (*.mat)")
        return fn


class DrawableLine:
    lock = None

    def __init__(self, ax, bg, x=None, y=None):
        self.ax = ax
        self.bg = bg
        if x is not None:
            self.dist = x
            self.y = y
            self.line = self.ax.plot(self.dist, self.y)
        else:
            self.dist = []
            self.y = []
            self.line = None
        self.press = None
        self.bg = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        # if event.inaxes != self.point.axes: return
        if DrawableLine.lock is not self:
            return

        if event.xdata is None or event.ydata is None:
            return

        if event.inaxes != self.ax:
            return

        if self.line is not None:
            self.line.set_animated(True)

        self.oldx = self.dist[:]
        self.oldy = self.y[:]

        self.dist.append(event.xdata)
        self.y.append(event.ydata)

        DrawableLine.lock = self
        self.press = True

    def update_line(self):
        # self.ax.figure.canvas.restore_region(self.bg)
        self.line.set_data(self.dist, self.y)
        self.ax.draw_artist(self.line)
        self.ax.figure.canvas.blit(self.ax.bbox)

    def on_motion(self, event):
        if not self.press or DrawableLine.lock is not self:
            return

        self.dist.append(event.xdata)
        self.y.append(event.ydata)

        if self.line is not None:
            self.update_line()
        else:
            self.line = self.ax.plot(self.dist, self.y)[0]

    def on_release(self, event):
        'on release we reset the press data'
        if DrawableLine.lock is not self:
            return

        self.press = None

        ind = np.argsort(self.dist)
        self.dist = [self.dist[i] for i in ind]
        self.y = [self.y[i] for i in ind]

        if self.line is None:
            self.line = self.ax.plot(self.dist, self.y)[0]
        else:
            self.update_line()
            self.line.set_animated(False)

    def revert(self):
        print('revert called')
        self.dist = self.oldx[:]
        self.y = self.oldy[:]

        self.line.set_animated(True)

        self.line.set_data(self.dist, self.y)

        self.line.figure.canvas.draw()
        self.line.figure.canvas.flush_events()
        self.update_line()
        self.line.set_animated(False)

    def write(self, fn):
        with open(fn, 'w') as fout:
            for x, y in zip(self.dist, self.y):
                fout.write('{:f}, {:f}\n'.format(x, y))


class LineOnRaster(DrawableLine):

    def __init__(self):
        pass


class LineList:

    def __init__(self, ax, lat=None, lon=None, linetype=DrawableLine):
        self.ax = ax
        self.linetype = linetype
        self.lines = []
        self.connect()
        self.waiting_for_renumber = False
        self.renumber = ''
        self.bg = self.ax.figure.canvas.copy_from_bbox(self.ax.bbox)

        self.lat = lat
        self.lon = lon
        self.numtext = None

        self._new()

    def connect(self):
        self.kid = self.ax.figure.canvas.mpl_connect('key_press_event', self.on_key)

    def _new(self, event=None):
        self.lines.append(self.linetype(self.ax, bg=self.bg))
        self.lines[-1].connect()
        self.linetype.lock = self.lines[-1]
        self._update_numtext()

    def _next(self, event):
        i = self.lines.index(self.linetype.lock)
        if i == len(self.lines) - 1:
            self._new(event)
        else:
            self.linetype.lock = self.lines[i + 1]
        self._update_numtext()

    def _previous(self, event):
        i = self.lines.index(self.linetype.lock)
        if i > 0:
            self.linetype.lock = self.lines[i - 1]
        self._update_numtext()

    def _undo(self, event):
        if self.linetype.lock is not None:
            self.linetype.lock.revert()

    def _save(self, event=None):
        for i, line in enumerate(self.lines):
            fn = 'traced_line_{:d}.csv'.format(i)
            line.write(fn)

    def _update_number(self, event):
        try:
            print('Renumbering to {:d}'.format(int(self.renumber)))
            self.linetype.lock = self.lines[int(self.renumber)]
        except:
            print('Failed to pickup {:s}'.format(self.renumber))
        finally:
            self.renumber = ''
            self.waiting_for_renumber = False

    def _update_numtext(self):
        if self.numtext is None or self.linetype.lock is None:
            return
        else:
            if self.linetype.lock is not None and self.linetype.lock.line is not None:
                color = self.linetype.lock.line.get_color()
            else:
                color = 'k'
            self.numtext.set_animated(True)
            self.numtext.set_text(str(self.lines.index(self.linetype.lock)))
            self.numtext.set_color(color)
            self.ax.draw_artist(self.numtext)
            self.numtext.set_animated(False)

    def on_key(self, event):
        if event.key == 'ctrl+n':
            self._new(event)
        elif event.key == 'ctrl+r':
            self.waiting_for_renumber = True
        elif event.key == 'ctrl+x':
            self._save(event)
        elif event.key == 'enter' and self.waiting_for_renumber:
            self._update_number
        elif event.key == 'ctrl+z':
            self._undo(event)
        elif event.key.isnumeric:
            self.renumber += event.key

    def add_buttons(self, button_dict):
        if 'New' in button_dict:
            button_dict['New'].on_clicked(self._new)
        if 'Next' in button_dict:
            button_dict['Next'].on_clicked(self._next)
        if 'Previous' in button_dict:
            button_dict['Previous'].on_clicked(self._previous)
        if 'Undo' in button_dict:
            button_dict['Undo'].on_clicked(self._undo)
        if 'Save' in button_dict:
            button_dict['Save'].on_clicked(self._save)

        if 'Number' in button_dict:
            button_dict['Number'].axis('off')
            self.numtext = button_dict['Number'].text(0.1, 0.9, '0', ha='left', va='top', transform=button_dict['Number'].transAxes, fontsize=24)


class zoom_factory:

    def __init__(self, ax, base_scale=2.):
        self.in_xlim = ax.get_xlim()
        self.in_ylim = ax.get_ylim()

        self.ax = ax
        self.base_scale = base_scale

    def connect(self):
        fig = self.ax.get_figure()  # get the figure of interest
        fig.canvas.mpl_connect('scroll_event', self.zoom_fun)

    def zoom_fun(self, event):
        # get the current x and y limits
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0]) * .5
        cur_yrange = (cur_ylim[1] - cur_ylim[0]) * .5
        xdata = event.xdata  # get event x location
        ydata = event.ydata  # get event y location
        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1. / self.base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = self.base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1.

        if (xdata - cur_xrange * scale_factor < self.in_xlim[0]) or (xdata + cur_xrange * scale_factor > self.in_xlim[1]) or (ydata - cur_yrange * scale_factor < self.in_ylim[0]) or (ydata + cur_yrange * scale_factor > self.in_ylim[1]):
            self.ax.set_xlim(self.in_xlim)
            self.ax.set_ylim(self.in_ylim)
        else:
            self.ax.set_xlim([xdata - cur_xrange * scale_factor, xdata + cur_xrange * scale_factor])
            self.ax.set_ylim([ydata - cur_yrange * scale_factor, ydata + cur_yrange * scale_factor])

        plt.draw()  # force re-draw
