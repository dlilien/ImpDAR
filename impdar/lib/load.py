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
from . import load_gssi, load_pulse_ekko, load_gprMax, load_gecko
from .RadarData import RadarData
try:
    from . import load_segy
    segy = True
except ImportError:
    segy = False

def load(filetype, fns):
    """Load a list of files of a certain type

    Parameters
    ----------
    filetype: str
        The type of file to load. Options are:
                        'pe' (pulse ekko)
                        'gssi' (from sir controller)
                        'gprMax' (synthetics)
                        'gecko' (St Olaf Radar)
                        'segy' (SEG Y)
                        'mat' (StODeep matlab format)
   fns: list
        List of files to load

    Returns
    -------
    RadarDataList: list of ~impdar.RadarData (or its subclasses)
        Objects with relevant radar information
    """
    if type(fns) not in {list, tuple}:
        fns = [fns]
    if filetype == 'gssi':
        dat = [load_gssi.load_gssi(fn) for fn in fns]
    elif filetype == 'pe':
        dat = [load_pulse_ekko.load_pe(fn) for fn in fns]
    elif filetype == 'mat':
        dat = [load_mat(fn) for fn in fns]
    elif filetype == 'gprMax':
        dat = [load_gprMax.load_gprMax(fn) for fn in fns]
    elif filetype == 'gecko':
        dat = [load_gecko.load_gecko(fn) for fn in fns]
    elif filetype == 'segy':
        if segy:
            dat = [load_segy.load_segy(fn) for fn in fns]
        else:
            raise ImportError('Failed to import segyio, cannot read segy')
    elif filetype == 'gprMax':
        dat = [load_gprMax.load_gprMax(fn) for fn in fns]
    else:
        raise ValueError('Unrecognized filetype')
    return dat


def load_and_exit(filetype, fn, *args, **kwargs):
    """Load a list of files of a certain type, save them as StODeep mat files, exit

    Parameters
    ----------
    filetype: str
        The type of file to load. Options are:
                        'pe' (pulse ekko)
                        'gssi' (from sir controller)
                        'gprMax' (synthetics)
                        'gecko' (St Olaf Radar)
                        'segy' (SEG Y)
                        'mat' (StODeep matlab format)
    fn: list or str
        List of files to load (or a single file)
    """
    if type(fn) not in {list, tuple}:
        fn = [fn]
    dat = load(filetype, fn)

    if 'o' in kwargs and kwargs['o'] is not None:
        out_fn = kwargs['o']
        if len(dat) > 1:
            raise ValueError('Cannot specify output with multiple inputs. Quitting without saving')
        dat[0].save(out_fn)
    else:
        for d, f in zip(dat, fn):
            out_fn = os.path.splitext(f)[0] + '_raw.mat'
            d.save(out_fn)


def load_mat(fn):
    """Load a .mat with radar info

    Just toss this in here so we have similar naming for

    Parameters
    ----------
    fn: str
        name of matlab file containing relevant variables
    """
    return RadarData(fn)
