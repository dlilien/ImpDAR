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
import glob
from . import load_mcords  # needs to be imported first and alone due to opaque h5py/netcdf4 error
from . import load_gssi, load_pulse_ekko, load_gprMax, load_olaf, load_segy, load_UoA
from . import load_delores, load_osu, load_stomat, load_ramac, load_bsi, load_tek
from ..RadarData import RadarData

# This should be updated as new functionality arrives
# executables that accept multiple ftypes should use this
# to figure out what the available options are
FILETYPE_OPTIONS = ['mat', 'pe', 'gssi', 'stomat', 'gprMax', 'gecko', 'segy', 'mcords_mat',
                    'mcords_nc', 'UoA_mat', 'UoA_h5', 'ramac', 'bsi', 'delores', 'osu',
                    'ramac', 'tek']


def load(filetype, fns_in, channel=1, t_srs=None, s_srs=None, *args, **kwargs):
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
        dat = []
        for fn in fns_in:
            # If a pulse ekko project file is input. We need to partition it first.
            if os.path.splitext(fn)[-1] == '.GPZ':
                bn_pe = os.path.splitext(fn)[0]
                print(fn, 'is a Pulse Ekko project file.\n' +
                      'We will partition it into the respective data and header files at ./' +
                      bn_pe)
                if not os.path.isdir(bn_pe):
                    os.mkdir(bn_pe)
                os.rename(fn, os.path.join(bn_pe, fn))
                os.chdir(bn_pe)
                load_pulse_ekko.partition_project_file(fn)
                os.rename(fn, os.path.join('..', fn))
                os.chdir('..')
                # load each file from within the project
                fns_partition = glob.glob(bn_pe+'/*.DT1')
                for fn_i in fns_partition:
                    dat.append(load_pulse_ekko.load_pe(fn_i))
            # If a standard pulse ekko file is input, try to load it
            else:
                try:
                    dat.append(load_pulse_ekko.load_pe(fn))
                except IOError:
                    print('Could not load ', fn, 'as a Pulse Ekko file.')
    elif filetype == 'mat':
        dat = [RadarData(fn) for fn in fns_in]
    elif filetype == 'stomat':
        dat = [load_stomat.load_stomat(fn, **kwargs) for fn in fns_in]
    elif filetype == 'gprMax':
        if load_gprMax.H5:
            dat = [load_gprMax.load_gprMax(fn) for fn in fns_in]
        else:
            raise ImportError('You need h5py for gprmax')
    elif filetype == 'bsi':
        # BSI data are slightly different since we may have multiple profiles per file
        if load_bsi.H5:
            if 'nans' in kwargs:
                nans = kwargs['nans']
            else:
                nans = 'interp'
            data_nestedlist = [load_bsi.load_bsi(fn, nans=nans) for fn in fns_in]
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
    elif filetype in ['UoA_mat', 'UoA_h5']:
        if load_UoA.H5:
            if 'gps_offset' in kwargs:
                gps_offset = kwargs['gps_offset']
            else:
                gps_offset = 0.0
            if filetype == 'UoA_mat':
                dat = [load_UoA.load_UoA_mat(fn, gps_offset=gps_offset) for fn in fns_in]
            elif filetype == 'UoA_h5':
                dat = []
                for fn in fns_in:
                    dat += load_UoA.load_UoA_h5(fn, gps_offset=gps_offset)
        else:
            raise ImportError('You need h5py for UoA')
    elif filetype == 'delores':
        dat = [load_delores.load_delores(fn, channel=channel) for fn in fns_in]
    elif filetype == 'osu':
        dat = [load_osu.load_osu(fns_in)]
    elif filetype == 'ramac':
        dat = [load_ramac.load_ramac(fn) for fn in fns_in]
    elif filetype == 'tek':
        dat = [load_tek.load_tek(fn) for fn in fns_in]
    else:
        raise ValueError('Unrecognized filetype')

    if s_srs is not None:
        try:
            for d in dat:
                d.get_ll(s_srs=s_srs)
        except ImportError:
            pass

    if t_srs is not None:
        try:
            for d in dat:
                d.get_projected_coords(t_srs=t_srs)
        except ImportError:
            pass

    return dat


def load_and_exit(filetype, fns_in, channel=1, t_srs=None, s_srs=None, o=None, *args, **kwargs):
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

    if filetype in ['osu', 'gecko']:
        # In this case, we do all at once since these radars split across files that should be merged
        rd_list = load(filetype, fns_in, channel=channel, t_srs=t_srs, *args, **kwargs)
        _save(rd_list, outpath=o)
    else:
        # Other filetypes can be safely done sequentially
        # If multiple ins, and output is not a dir, we have a problem
        if (len(fns_in) > 1) and (o is not None) and (not os.path.isdir(o)):
            raise FileNotFoundError('The output directory does not exist')

        for fn_i in fns_in:
            rd_list = load(filetype, fn_i, channel=channel, t_srs=t_srs, s_srs=s_srs, *args, **kwargs)
            _save(rd_list, outpath=o)


def _save(rd_list, outpath=None):
    """Save a list of RadarData objects with optional output directory."""
    if outpath is not None:
        if len(rd_list) > 1:
            for rd in rd_list:
                fn_out = os.path.join(outpath, os.path.split(os.path.splitext(rd.fn)[0] + '_raw.mat')[-1])
                rd.save(fn_out)
        elif os.path.isdir(outpath):
            fn_out = outpath + os.path.splitext(rd_list[0].fn)[0] + '_raw.mat'
            rd_list[0].save(fn_out)
        else:
            fn_out = outpath
            rd_list[0].save(fn_out)
    else:
        for rd in rd_list:
            fn_out = os.path.splitext(rd.fn)[0] + '_raw.mat'
            rd.save(fn_out)
