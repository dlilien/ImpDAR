#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 License.

"""
Do some filetype conversions. Created mainly to have a .DZG to .shp convertsion
"""

import os
from .load import load_mat
from . import load_gssi, load_pulse_ekko

def convert(fn, out_fmt, t_srs='wgs84', in_fmt=None, *args, **kwargs):
    # Basic check on the conversion being implemented. This is really simple because I'm not yet allowing conversion from one proprietary form to another
    if out_fmt not in ['shp', 'mat']:
        raise ValueError('Can only convert to shp or mat')

    # Treat this like batch input always
    if type(fn) == 'str':
        fn = [fn]

    # Prepare a list of filetypes so we can take diverse inputs simultaneously
    if in_fmt is None:
        loaders = [0 for i in fn]
        for i, f in enumerate(fn):
            if f[-4:] == '.mat':
                loaders[i] = load_mat
            if f[-4:] == '.DZT':
                loaders[i] = load_gssi.load_gssi
            else:
                raise ValueError('Unrecognized file extension {:s}'.format(f[-4:]))
    else:
        if in_fmt == 'mat':
            loaders = [load_mat for i in fn]
        if in_fmt == 'gssi':
            loaders = [load_gssi.load_gssi for i in fn]
        if in_fmt == 'pe':
            loaders = [load_pulse_ekko.load_pe for i in fn]


    # Now actually load the data
    data = [loader(f) for loader, f in zip(loaders, fn)]

    # Convert
    if out_fmt == 'mat':
        for loader, fn, dat in zip(loaders, fn, data):
            # Guard against silly re-write
            if loader == load_mat and out_fmt == 'mat':
                continue
            out_fn = os.path.splitext(f)[0] + '.mat'
            dat.save(out_fn)
    elif out_fmt == 'shp':
        for loader, fn, dat in zip(loaders, fn, data):
            out_fn = os.path.splitext(f)[0] + '.shp'
            dat.output_shp(out_fn, t_srs=t_srs)




if __name__ == '__main__':
    pass
