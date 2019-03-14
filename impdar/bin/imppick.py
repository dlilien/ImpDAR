#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
An executable to start the picker.
"""


import argparse

from impdar.gui import pickgui

from impdar.lib import load

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
from matplotlib import rc
rc('text', usetex=False) 


def main():
    parser = _get_args()
    args = parser.parse_args()
    radardata = load.load('mat', [args.fn])[0]
    ip = pickgui.pick(radardata, guard_save=True, xd=args.xd, yd=args.yd)


def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fn', type=str, help='The file to pick. At least for now you can only do one at a time.')
    parser.add_argument('-xd', action='store_true', help='Distance on the x')
    parser.add_argument('-yd', action='store_true', help='Depth on the y')
    return parser


if __name__ == '__main__':
    main()
