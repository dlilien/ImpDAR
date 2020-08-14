#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 License.
"""Do some filetype conversions. Incomplete."""

import os
from .RadarData import RadarData
from .load import load_gssi, load_pulse_ekko, load_segy, load


def convert(fns_in, out_fmt, t_srs=None, in_fmt=None, *args, **kwargs):
    """Convert between formats. Mainly used to create shps and sgy files."""
    # Basic check on the conversion being implemented.
    # This is really simple because I'm not converting from one proprietary
    # form to another
    if t_srs == 'wgs84':
        t_srs = 'EPSG:4326'

    if out_fmt not in ['shp', 'mat', 'sgy']:
        raise ValueError('Can only convert to shp, mat, or sgy')

    # Treat this like batch input always
    if not isinstance(fns_in, (tuple, list)):
        fns_in = [fns_in]

    # Prepare a list of filetypes so we can take diverse inputs simultaneously
    if in_fmt is None:
        loaders = [0 for i in fns_in]
        for i, f_i in enumerate(fns_in):
            if f_i[-4:] == '.mat':
                loaders[i] = RadarData
            elif f_i[-4:] == '.DZT':
                loaders[i] = load_gssi.load_gssi
            elif f_i[-4:] == '.DT1':
                loaders[i] = load_pulse_ekko.load_pe
            elif f_i[-4:] == '.sgy':
                if not load_segy.SEGY:
                    raise ImportError('You cannot use segy without segyio installed!')
                loaders[i] = load_segy.load_segy
            else:
                raise ValueError('Unrecognized file extension {:s}'.format(f_i[-4:]))
    else:
        loaders = [lambda x: load(in_fmt, x)[0] for i in fns_in]

    # Now actually load the data
    for loader, fn_i in zip(loaders, fns_in):
        data = loader(fn_i)

        # Convert
        if out_fmt == 'mat':
            # Guard against silly re-write
            if loader == RadarData and out_fmt == 'mat':
                raise ValueError('You are trying a blank conversion that will cause an overwrite...')
            fn_out = os.path.splitext(data.fn)[0] + '.mat'
            data.save(fn_out)
        elif out_fmt == 'shp':
            fn_out = os.path.splitext(data.fn)[0] + '.shp'
            data.output_shp(fn_out, t_srs=t_srs)
        elif out_fmt == 'sgy':
            if not load_segy.SEGY:
                raise ImportError('You cannot use segy without segyio installed!')
            fn_out = os.path.splitext(data.fn)[0] + '.sgy'
            data.save_as_segy(fn_out)


if __name__ == '__main__':
    pass
