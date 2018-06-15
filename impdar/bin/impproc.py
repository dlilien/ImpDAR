#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make an executable for single actions of impulse radar processing.

All functionality probably overlaps with impdar, but the call is much cleaner. You get a lot more flexibility on things like keyword arguments. However, you are limited to one processing step.
"""

import os.path
import argparse
from impdar.lib.load import load
from impdar.lib.process import concat
from impdar.lib.gpslib import interp as interpdeep


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

    # Concatenate
    add_simple_procparser(subparsers, 'cat', 'Concatenate the data', concat, defname='cat')

    # Restack
    parser_restack = add_procparser(subparsers, 'restack', 'Restack to interval', restack, defname='restacked')
    parser_restack.add_argument('traces', type=int, help='Number of traces to stack. Must be an odd number')
    add_def_args(parser_restack)

    # Range gain
    parser_rgain = add_procparser(subparsers, 'rgain', 'Add a range gain', rgain, defname='rgain')
    parser_rgain.add_argument('-slope', type=float, default=0.1, help='Slope of the linear range gain. Default 0.1')
    add_def_args(parser_rgain)

    # Range gain
    parser_agc = add_procparser(subparsers, 'agc', 'Add an automatic gain', agc, defname='agc')
    parser_agc.add_argument('-window', type=int, default=50, help='Number of samples to average')
    add_def_args(parser_agc)

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

    # Reinterpolate GPS
    parser_interp = add_procparser(subparsers, 'interp', 'Reinterpolate GPS', interp, defname='interp')
    parser_interp.add_argument('spacing', type=float, help='New spacing of radar traces, in meters')
    parser_interp.add_argument('gps_fn', type=str, help='CSV or mat file containing the GPS information. .csv and .txt files are assumed to be csv, .mat are mat')
    parser_interp.add_argument('--offset', type=float, default=0.0, help='Offset from GPS time to radar time')
    parser_interp.add_argument('--minmove', type=float, default=1.0e-2, help='Minimum movement to not be stationary')
    add_def_args(parser_interp)

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

    if args.name == 'cat':
        radar_data = concat(radar_data)
        bn = os.path.splitext(args.fns[0])[0]
        args.fns = [bn + '.mat']
    elif args.name == 'interp':
        interp(radar_data, **vars(args))
        
    else:
        for dat in radar_data:
            args.func(dat, **vars(args))

    if args.o is not None:
        if len(radar_data) > 1:
            for d, f in zip(radar_data, args.fns):
                bn = os.path.split(os.path.splitext(f)[0])[1]
                if bn[-4:] == '_raw':
                    bn = bn[:-4]
                out_fn = args.o + '/' + bn + '_{:s}.mat'.format(args.name)
                d.save(out_fn)
        else:
            out_fn = args.o
            radar_data[0].save(out_fn)
    else:
        for d, f in zip(radar_data, args.fns):
            bn = os.path.splitext(f)[0]
            if bn[-4:] == '_raw':
                bn = bn[:-4]
            out_fn = bn + '_{:s}.mat'.format(args.name)
            d.save(out_fn)


def hfilt(dat, start_trace=0, end_trace=-1, **kwargs):
    dat.hfilt(ftype='hfilt', bounds=(start_trace, end_trace))


def ahfilt(dat, **kwargs):
    dat.hfilt(ftype='adaptive')


def rev(dat, **kwargs):
    dat.reverse()


def vbp(dat, low=1, high=10000, **kwargs):
    dat.vertical_band_pass(low, high)


def crop(dat, lim=0, top_or_bottom='top', dimension='snum', **kwargs):
    dat.crop(lim, top_or_bottom=top_or_bottom, dimension=dimension)


def nmo(dat, ant_sep=0.0, uice=1.69e8, uair=3.0e8, **kwargs):
    dat.nmo(ant_sep, uice=uice, uair=uair)


def restack(dat, traces=1, **kwargs):
    dat.restack(traces)


def rgain(dat, slope=0.1, **kwargs):
    dat.rangegain(slope)


def agc(dat, window=50, scale_factor=50, **kwargs):
    dat.agc(window=window, scaling_factor=scale_factor)


def interp(dats, spacing, gps_fn, offset=0.0, minmove=1.0e-2, **kwargs):
    interpdeep(dats, spacing, fn=gps_fn, offset=offset, min_movement=minmove)


if __name__ == '__main__':
    main()
