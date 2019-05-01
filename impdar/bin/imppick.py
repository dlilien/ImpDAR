#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
An executable to start the picker.
"""

import sys
import argparse
from PyQt5 import QtWidgets

from impdar.gui import pickgui

from impdar.lib import load

from matplotlib import rc
rc('text', usetex=False) 


def pick(radardata, guard_save=True, xd=False, yd=False):
    if xd:
        x = 'dist'
    else:
        x = 'tnum'
    if yd:
        y = 'depth'
    else:
        y = 'twtt'

    if not hasattr(radardata, 'picks') or radardata.picks is None:
        radardata.picks = RadarData.Picks(radardata)

    app = QtWidgets.QApplication(sys.argv)
    ip = pickgui.InteractivePicker(radardata, xdat=x, ydat=y)
    ip.show()
    sys.exit(app.exec_())


def main():
    parser = _get_args()
    args = parser.parse_args()
    radardata = load.load('mat', [args.fn])[0]
    ip = pick(radardata, guard_save=True, xd=args.xd, yd=args.yd)


def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fn', type=str, help='The file to pick. At least for now you can only do one at a time.')
    parser.add_argument('-xd', action='store_true', help='Distance on the x')
    parser.add_argument('-yd', action='store_true', help='Depth on the y')
    return parser


if __name__ == '__main__':
    main()
