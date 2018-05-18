#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
A wrapper around the other loading utilities
"""
import os.path
from . import load_gssi, load_pulse_ekko
from .stodeep_fmt import StODeep


def load(filetype, fn, *args, **kwargs):
    if filetype == 'gssi':
        dat = load_gssi.load_gssi(fn)
    elif filetype == 'pe':
        dat = load_pulse_ekko.load_pe(fn)
    elif filetype == 'mat':
        dat = load_mat(fn)
    else:
        raise Exception('Unrecognized filetype')
    return dat


def load_and_exit(filetype, fn, *args, **kwargs):
    dat = load(filetype, fn, *args, **kwargs)

    if 'o' in kwargs and kwargs['o'] is not None:
        out_fn = kwargs['o']
    else:
        out_fn = os.path.splitext(fn)[0] + '_raw.mat'
    dat.save(out_fn)


def load_mat(fn):
    return StODeep(fn)
