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
from . import load_gssi, load_pulse_ekko, load_gprMax, load_olaf, load_mcords_nc, load_mcords_mat, load_segy
from ..RadarData import RadarData

# This should be updated as new functionality arrives
# executables that accept multiple ftypes should use this
# to figure out what the available options are
FILETYPE_OPTIONS = ['mat', 'pe', 'gssi', 'gprMax', 'gecko', 'segy', 'mcords_mat', 'mcords_nc']


def load(filetype, fns_in, channel=1):
    """Load a list of files of a certain type

    Parameters
    ----------
    filetype: str
        The type of file to load. Options are:
                        'mat' (StODeep matlab format)
                        'pe' (pulse ekko)
                        'gssi' (from sir controller)
                        'gprMax' (synthetics)
                        'gecko' (St Olaf Radar)
                        'segy' (SEG Y)
                        'mcords_nc' (MCoRDS netcdf)
                        'mcords_mat' (MCoRDS matlab format)
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
    if not isinstance(fns_in, (list, tuple)):
        fns_in = [fns_in]

    if filetype == 'gssi':
        dat = [load_gssi.load_gssi(fn) for fn in fns_in]
    elif filetype == 'pe':
        dat = [load_pulse_ekko.load_pe(fn) for fn in fns_in]
    elif filetype == 'mat':
        dat = [RadarData(fn) for fn in fns_in]
    elif filetype == 'gprMax':
        if load_gprMax.H5:
            dat = [load_gprMax.load_gprMax(fn) for fn in fns_in]
        else:
            raise ImportError('You need h5py for gprmax')
    elif filetype == 'gecko':
        # Slightly different because we assume that we want to concat
        dat = [load_olaf.load_olaf(fns_in, channel=channel)]
    elif filetype == 'segy':
        if load_segy.SEGY:
            dat = [load_segy.load_segy(fn) for fn in fns_in]
        else:
            raise ImportError('Failed to import segyio, cannot read segy')
    elif filetype == 'gprMax':
        dat = [load_gprMax.load_gprMax(fn) for fn in fns_in]
    elif filetype == 'mcords_nc':
        if load_mcords_nc.NC:
            dat = [load_mcords_nc.load_mcords_nc(fn) for fn in fns_in]
        else:
            raise ImportError('You need netCDF4 in order to read the MCoRDS files')
    elif filetype == 'mcords_mat':
        dat = [load_mcords_mat.load_mcords_mat(fn) for fn in fns_in]
    else:
        raise ValueError('Unrecognized filetype')
    return dat


def load_and_exit(filetype, fns_in, channel=1, *args, **kwargs):
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
                        'mcords_nc' (MCoRDS netcdf)
                        'mcords_mat' (MCoRDS matlab format)
                        'mat' (StODeep matlab format)
    fn: list or str
        List of files to load (or a single file)
    channel: Receiver channel that the data were recorded on
        This is primarily for the St. Olaf HF data
    """
    if not isinstance(fns_in, (list, tuple)):
        fns_in = [fns_in]
    dat = load(filetype, fns_in, channel=channel)

    if filetype == 'gecko' and len(fns_in) > 1:
        f_common = fns_in[0]
        for i in range(1, len(fns_in)):
            f_common = _common_start(f_common, fns_in[i]).rstrip('[')
        fn_out = os.path.splitext(f_common)[0] + '_raw.mat'
        if 'o' in kwargs and kwargs['o'] is not None:
            fn_out = os.path.join(kwargs['o'], os.path.split(fn_out)[-1])
        dat[0].save(fn_out)
    elif 'o' in kwargs and kwargs['o'] is not None:
        if len(fns_in) > 1:
            for d_i, f_i in zip(dat, fns_in):
                fn_out = os.path.join(kwargs['o'], os.path.split(os.path.splitext(f_i)[0] + '_raw.mat')[-1])
                d_i.save(fn_out)
        else:
            fn_out = kwargs['o']
            dat[0].save(fn_out)
    else:
        for d_i, f_i in zip(dat, fns_in):
            fn_out = os.path.splitext(f_i)[0] + '_raw.mat'
            d_i.save(fn_out)


def _common_start(string_a, string_b):
    """ returns the longest common substring from the beginning of sa and sb

    from https://stackoverflow.com/questions/18715688/find-common-substring-between-two-strings
    """
    def _iter():
        for char_a, char_b in zip(string_a, string_b):
            if char_a == char_b:
                yield char_a
            else:
                return

    return ''.join(_iter())
