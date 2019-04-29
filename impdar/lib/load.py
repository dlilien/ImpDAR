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
from . import load_gssi, load_pulse_ekko, load_olaf
try:
    from . import load_segy
    segy = True
except ImportError:
    segy = False
from .RadarData import RadarData


def load(filetype, fns, nchan=1):
    """Load a list of files of a certain type

    Parameters
    ----------
    filetype: str
        The type of file to load. Options are 'pe' (pulse ekko), 'gssi' (from sir controller) or 'mat' (StODeep matlab format)
    fns: list
        List of files to load
    nchan: int, optional
        Channel number for multichannel data

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
    elif filetype == 'olaf':
        # Slightly different because we assume that we want to concat
        dat = [load_olaf.load_olaf(fns, Channel_Num=nchan)]
    elif filetype == 'segy':
        if segy:
            dat = [load_segy.load_segy(fn) for fn in fns]
        else:
            raise ImportError('Failed to import segyio, cannot read segy')
    else:
        raise ValueError('Unrecognized filetype')
    return dat


def load_and_exit(filetype, fn, *args, nchan=1, **kwargs):
    """Load a list of files of a certain type, save them as StODeep mat files, exit

    Parameters
    ----------
    filetype: str
        The type of file to load. Options are 'pe' (pulse ekko), 'gssi' (from sir controller) or 'mat' (StODeep matlab format
    fn: list or str
        List of files to load (or a single file)
    """
    if type(fn) not in {list, tuple}:
        fn = [fn]
    dat = load(filetype, fn, nchan=nchan)

    if 'o' in kwargs and kwargs['o'] is not None:
        out_fn = kwargs['o']
        if len(dat) > 1:
            raise ValueError('Cannot specify output with multiple inputs. Quitting without saving')
        dat[0].save(out_fn)
    elif filetype == 'olaf' and len(fn) > 1:
        f = fn[0]
        for i in range(1, len(fn)):
            f = common_start(f, fn[i]).rstrip('[')
        out_fn = os.path.splitext(f)[0] + '_raw.mat'
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


def common_start(sa, sb):
    """ returns the longest common substring from the beginning of sa and sb
    
    from https://stackoverflow.com/questions/18715688/find-common-substring-between-two-strings
    """
    def _iter():
        for a, b in zip(sa, sb):
            if a == b:
                yield a
            else:
                return

    return ''.join(_iter())
