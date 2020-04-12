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
import numpy as np
from . import load_gssi, load_pulse_ekko, load_gprMax, load_olaf, load_mcords, load_segy, load_UoA_mat, load_ramac, load_bsi
from . import load_delores, load_osu, load_stomat
from ..RadarData import RadarData

# This should be updated as new functionality arrives
# executables that accept multiple ftypes should use this
# to figure out what the available options are
FILETYPE_OPTIONS = ['mat', 'pe', 'gssi','stomat', 'gprMax', 'gecko', 'segy',
                    'mcords_mat', 'mcords_nc', 'UoA_mat', 'ramac', 'bsi', 'delores', 'osu', 'ramac']


def load(filetype, fns_in, channel=1, *args, **kwargs):
    """Load a list of files of a certain type

    Parameters
    ----------
    filetype: str
        The type of file to load.
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
    elif filetype == 'stomat':
        dat = [load_stomat.load_stomat(fn) for fn in fns_in]
    elif filetype == 'gprMax':
        if load_gprMax.H5:
            dat = [load_gprMax.load_gprMax(fn) for fn in fns_in]
        else:
            raise ImportError('You need h5py for gprmax')
    elif filetype == 'bsi':
        # BSI data are slightly different since we may have multiple profiles per file
        if load_bsi.H5:
            data_nestedlist = [load_bsi.load_bsi(fn) for fn in fns_in]
            dat = []
            for data in data_nestedlist:
                dat.extend(data)
        else:
            raise ImportError('You need h5py for bsi')
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
        if load_mcords.NC:
            dat = [load_mcords.load_mcords_nc(fn) for fn in fns_in]
        else:
            raise ImportError('You need netCDF4 in order to read the MCoRDS files')
    elif filetype == 'mcords_mat':
        dat = [load_mcords.load_mcords_mat(fn) for fn in fns_in]
    elif filetype == 'UoA_mat':
        if load_UoA_mat.H5:
            if 'gps_offset' in kwargs:
                gps_offset = kwargs['gps_offset']
            else:
                gps_offset = 0.0
            dat = [load_UoA_mat.load_UoA_mat(fn, gps_offset=gps_offset) for fn in fns_in]
        else:
            raise ImportError('You need h5py for UoA_mat')
    elif filetype == 'delores':
        dat = [load_delores.load_delores(fn, channel=channel) for fn in fns_in]
    elif filetype == 'osu':
        dat = [load_osu.load_osu(fns_in)]
    elif filetype == 'ramac':
        dat = [load_ramac.load_ramac(fn) for fn in fns_in]
    else:
        raise ValueError('Unrecognized filetype')
    return dat


def load_and_exit(filetype, fns_in, channel=1, t_srs=None, *args, **kwargs):
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
    channel: int, optional
        Receiver channel that the data were recorded on
        This is primarily for the St. Olaf HF data
    t_srs: str, optional
        Convert to this coordinate system. Requires GDAL.
    """

    if not isinstance(fns_in, (list, tuple)):
        fns_in = [fns_in]

    # If a pulse ekko project file is input. We need to partition it first
    if filetype =='pe' and np.any(np.array([os.path.splitext(fn)[-1] for fn in fns_in]) == '.GPZ'):
        for fn in fns_in:
            if os.path.splitext(fn)[-1] != '.GPZ':
                print(fn,'is NOT a Pulse Ekko Project File but we are partitioning now.'+\
                        'Load it on its own.')
                continue
            else:
                bn_pe = os.path.splitext(fn)[0]
                print(fn,'is a Pulse Ekko project file.\n'+\
                'We will partition it into the respective data and header files at ./'+\
                bn_pe)
                if not os.path.isdir(bn_pe):
                    os.mkdir(bn_pe)
                os.rename(fn,bn_pe+'/'+fn)
                os.chdir(bn_pe)
                load_pulse_ekko.partition_project_file(fn)
                os.rename(fn,'../'+fn)
        return
    else:
        dat = load(filetype, fns_in, channel=channel, *args, **kwargs)

    if t_srs is not None:
        try:
            for d in dat:
                d.get_projected_coords(t_srs=t_srs)
        except ImportError:
            pass


    if (filetype == 'gecko' or filetype == 'osu') and len(fns_in) > 1:
        f_common = fns_in[0]
        for i in range(1, len(fns_in)):
            f_common = _common_start(f_common, fns_in[i]).rstrip('[')
        fn_out = os.path.splitext(f_common)[0] + '_raw.mat'
        if 'o' in kwargs and kwargs['o'] is not None:
            fn_out = os.path.join(kwargs['o'], os.path.split(fn_out)[-1])
        dat[0].save(fn_out)
    elif 'o' in kwargs and kwargs['o'] is not None:
        if len(fns_in) > 1:
            for d_i in dat:
                fn_out = os.path.join(kwargs['o'], os.path.split(os.path.splitext(d_i.fn)[0] + '_raw.mat')[-1])
                d_i.save(fn_out)
        elif os.path.isdir(kwargs['o']):
            fn_out = kwargs['o'] + os.path.splitext(fns_in[0])[0] + '_raw.mat'
            dat[0].save(fn_out)
        else:
            fn_out = kwargs['o']
            dat[0].save(fn_out)
    else:
        for d_i in dat:
            fn_out = os.path.splitext(d_i.fn)[0] + '_raw.mat'
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
