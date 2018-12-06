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

from impdar.lib import pickgui, load

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
from matplotlib import rc
rc('text', usetex=False) 


def main():
    parser = _get_args()
    args = parser.parse_args()
    radardata = load.load('mat', [args.fn])[0]
    ip = pickgui.pick(radardata, guard_save=True)


def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fn', type=str, help='The file to pick. At least for now you can only do one at a time.')
    return parser


if __name__ == '__main__':
    main()
