#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""An executable to start the picker."""

import sys
import argparse
from PyQt5 import QtWidgets

from matplotlib import rc

from impdar.gui import pickgui
from impdar.lib import load, Picks

rc('text', usetex=False)


def pick(radardata, xd=False, yd=False):
    """Fire up the picker."""
    if xd:
        x = 'dist'
    else:
        x = 'tnum'
    if yd:
        y = 'depth'
    else:
        y = 'twtt'

    if not hasattr(radardata, 'picks') or radardata.picks is None:
        radardata.picks = Picks.Picks(radardata)

    app = QtWidgets.QApplication(sys.argv)
    ip = pickgui.InteractivePicker(radardata, xdat=x, ydat=y)
    ip.show()
    sys.exit(app.exec_())


def main():
    """Get arguments, start picking."""
    parser = _get_args()
    args = parser.parse_args(sys.argv[1:])
    radardata = load.load('mat', [args.fn])[0]
    pick(radardata, xd=args.xd, yd=args.yd)


def _get_args():
    """Get the parser for arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('fn',
                        type=str,
                        help='The file to pick. One file at a time.')
    parser.add_argument('-xd', action='store_true', help='Distance on the x')
    parser.add_argument('-yd', action='store_true', help='Depth on the y')
    return parser


if __name__ == '__main__':
    main()
