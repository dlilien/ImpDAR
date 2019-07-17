#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
from PyQt5 import QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as Canvas, NavigationToolbar2QT as NavigationToolbar)
import matplotlib

# Ensure using PyQt5 backend since if not QT we will crash
matplotlib.use('QT5Agg')


class MplCanvas(Canvas):
    """The figure canvas itself. I do a sort of dumb 'tight_layout' call but I haven't worked out something better"""

    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        Canvas.updateGeometry(self)
        self.fig.tight_layout(pad=0.1, rect=[0.05, 0.05, 1, 1])


class MplFigCanvasWidget(QtWidgets.QWidget):
    """This is the big box containing the plot and also the MPL toolbar"""

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
        self.canvas = MplCanvas()                  # Create canvas object
        self.canvas.setFocus()
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)

        self.vbl = QtWidgets.QVBoxLayout()         # Set box for plotting
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.mpl_toolbar)
        self.setLayout(self.vbl)
