#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
A wrapper around the other loading utilities
"""
import os.path
from . import load_gssi, load_pulse_ekko, load_gprMax, load_olaf, load_mcords_nc
from .RadarData import RadarData
try:
    from . import load_segy
    segy = True
except ImportError:
    segy = False


def load(filetype, fns, channel=1):
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
    channel: Receiver channel that the data were recorded on
        This is primarily for the St. Olaf HF data

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
        # Slightly different because we assume that we want to concat
        dat = [load_olaf.load_olaf(fns, channel=channel)]
    elif filetype == 'segy':
        if segy:
            dat = [load_segy.load_segy(fn) for fn in fns]
        else:
            raise ImportError('Failed to import segyio, cannot read segy')
    elif filetype == 'gprMax':
        dat = [load_gprMax.load_gprMax(fn) for fn in fns]
    elif filetype == 'mcords':
        if load_mcords_nc.nc:
            dat = [load_mcords_nc.load_mcords_nc(fn) for fn in fns]
        else:
            raise ImportError('You need netCDF4 in order to read the MCoRDS files')
    else:
        raise ValueError('Unrecognized filetype')
    return dat


def load_and_exit(filetype, fn, channel=1, *args, **kwargs):
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
    channel: Receiver channel that the data were recorded on
        This is primarily for the St. Olaf HF data
    """
    if type(fn) not in {list, tuple}:
        fn = [fn]
    dat = load(filetype, fn, channel=channel)

    if 'o' in kwargs and kwargs['o'] is not None:
        out_fn = kwargs['o']
        if len(dat) > 1:
            raise ValueError('Cannot specify output with multiple inputs. Quitting without saving')
        dat[0].save(out_fn)
    elif filetype == 'gecko' and len(fn) > 1:
        f = fn[0]
        for i in range(1, len(fn)):
            f = _common_start(f, fn[i]).rstrip('[')
        out_fn = os.path.splitext(f)[0] + '_raw.mat'
        dat[0].save(out_fn)
    else:
        for d, f in zip(dat, fn):
            if f[-3:] == 'g00':
                out_fn = os.path.splitext(f)[0] + '_g00_raw.mat'
            else:
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


def _common_start(sa, sb):
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
