#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make an executable for single actions of impulse radar processing.

All functionality probably overlaps with impdar, but the call is much cleaner. However, you are limited to one processing step.
"""

import os.path
import argparse
from impdar.lib.load import load


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Choose a processing step')

    # Horizontal window filter
    parser_hfilt = add_procparser(subparsers, 'hfilt', 'Horizontally filter the data by subtracting the average trace from a window', hfilt, defname='hfilted')
    parser_hfilt.add_argument('start_trace', type=int, help='First trace of representative subset')
    parser_hfilt.add_argument('end_trace', type=int, help='Last trace of representative subset')
    add_def_args(parser_hfilt)

    # Adaptive horizontal filter
    add_simple_procparser(subparsers, 'ahfilt', 'Horizontally filter the data adaptively', ahfilt, defname='ahfilt')

    # Simply reverse the files
    add_simple_procparser(subparsers, 'rev', 'Reverse the data', rev, defname='rev')

    # Vertical bandpass
    parser_vbp = add_procparser(subparsers, 'vbp', 'Vertically bandpass the data', vbp, defname='bandpassed')
    parser_vbp.add_argument('low_MHz', type=float, help='Lowest frequency passed (in MHz)')
    parser_vbp.add_argument('high_MHz', type=float, help='Highest frequency passed (in MHz)')
    add_def_args(parser_vbp)

    # Crop in the vertical 
    parser_crop = add_procparser(subparsers, 'crop', 'Crop the data in the vertical', crop, defname='cropped')
    parser_crop.add_argument('top_or_bottom', choices=['top', 'bottom'], help='Remove from the top or bottom')
    parser_crop.add_argument('dimension', choices=['snum', 'twtt', 'depth'], help='Set the bound in terms of snum (sample number), twtt (two way travel time in microseconds), or depth (m, calculated using the nmo_depth or a light speed of 1.69e8m/s if it doesn\'t')
    parser_crop.add_argument('lim', type=float, help='The cutoff value')
    add_def_args(parser_crop)

    # Normal move-out
    parser_nmo = add_procparser(subparsers, 'nmo', 'Normal move-out correction', nmo, defname='nmo')
    parser_nmo.add_argument('ant_sep', type=float, help='Antenna separation')
    parser_nmo.add_argument('--uice', type=float, default=1.69e8, help='Speed of light in ice in m/s (default 1.69e8)')
    parser_nmo.add_argument('--uair', type=float, default=3.0e8, help='Speed of light in air in m/s (default 3.0e8)')
    add_def_args(parser_nmo)
    return parser


def add_simple_procparser(subparsers, name, helpstr, func, defname='proc'):
    parser = add_procparser(subparsers, name, helpstr, func, defname=defname)
    add_def_args(parser)
    return parser


def add_procparser(subparsers, name, helpstr, func, defname='proc'):
    parser = subparsers.add_parser(name, help=helpstr)
    parser.set_defaults(func=func, name=defname)
    return parser


def add_def_args(parser):
    parser.add_argument('fns', type=str, nargs='+', help='The files to process')
    parser.add_argument('-o', type=str, help='Output to this file (or folder if multiple inputs)')
    parser.add_argument('-pe', action='store_true', help='Inputs are pulse ekko files')
    parser.add_argument('-gssi', action='store_true', help='Inputs are gssi files')


def main():
    parser = _get_args()
    args = parser.parse_args()
    print(args)
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])
        return

    if args.gssi and args.pe:
        raise ValueError('Input cannot be both pulse-ekko and gssi')
    if args.gssi:
        radar_data = load('gssi', args.fns)
    elif args.pe:
        radar_data = load('pe', args.fns)
    else:
        radar_data = load('mat', args.fns)
    
    for dat in radar_data:
        args.func(dat, **vars(args))

    if args.o is not None:
        if len(radar_data) > 1:
            for d, f in zip(radar_data, args.fns):
                out_fn = args.o + '/' + os.path.split(os.path.splitext(f)[0])[1] + '_{:s}.mat'.format(args.name)
                d.save(out_fn)
        else:
            radar_data[0].save(out_fn)
    else:
        for d, f in zip(radar_data, args.fns):
            out_fn = os.path.splitext(f)[0] + '_{:s}.mat'.format(args.name)
            d.save(out_fn)


def hfilt(dat, start_trace=0, end_trace=-1, **kwargs):
    dat.hfilt(ftype='hfilt', bounds=(start_trace, end_trace))


def ahfilt(dat, **kwargs):
    dat.hfilt(ftype='adaptive')


def rev(dat, **kwargs):
    dat.reverse()


def vbp(dat, low=1, high=10000, **kwargs):
    dat.vertical_band_pass(low, high)


def crop(dat, top_or_bottom, dimension, lim, **kwargs):
    dat.crop(top_or_bottom, dimension, lim)


def nmo(dat, ant_sep=0.0, uice=1.69e8, uair=3.0e8, **kwargs):
    dat.nmo(ant_sep, uice=uice, uair=uair)


if __name__ == '__main__':
    main()
