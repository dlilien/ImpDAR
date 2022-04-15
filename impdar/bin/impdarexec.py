#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
#
"""The primary impdar executable, called as `impdar`."""
import sys
import argparse
from impdar.lib import load, process, plot, convert


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_load = subparsers.add_parser('load', help='Load data')
    parser_load.set_defaults(func=load.load_and_exit)
    parser_load.add_argument('filetype', type=str,
                             help='Type of file',
                             choices=load.FILETYPE_OPTIONS)
    parser_load.add_argument('fns_in',
                             type=str,
                             nargs='+',
                             help='File(s) to load')
    parser_load.add_argument('-channel', type=int, default=1,
                             help='Receiver channel to load this is primarily for the St. Olaf HF data.')
    parser_load.add_argument('-gps_offset',
                             type=float,
                             help='Offset of GPS and data times for UoA_mat',
                             default=0.0)
    parser_load.add_argument('-t_srs', type=str, default=None,
                             help='Convert to this coordinate reference system. (GDAL required), default UTM')
    parser_load.add_argument('-s_srs', type=str, default=None,
                             help='Convert from this system. (GDAL required), default UTM')
    parser_load.add_argument('-o', type=str, help='Write to this filename')
    parser_load.add_argument('--nans', type=str, choices=['interp', 'delete'], default=None,
                             help='Interpolate or delete bad GPS. Only used by BSI.')
    parser_load.add_argument('-dname',
                             type=str,
                             help='Name of data field',
                             default='data')

    # Options for processing data
    parser_proc = subparsers.add_parser('proc', help='Process data')
    parser_proc.set_defaults(func=process.process_and_exit)
    parser_load.add_argument('--filetype',
                             type=str,
                             help='Type of file',
                             default='mat',
                             choices=load.FILETYPE_OPTIONS)
    parser_proc.add_argument('-cat',
                             action='store_true',
                             help='Concatenate the files')
    parser_proc.add_argument('-vbp',
                             nargs=2,
                             type=float,
                             help='Bandpass the data vertically at \
                                 low (MHz) and high (MHz)')
    parser_proc.add_argument('-hfilt',
                             nargs=2,
                             type=int,
                             help='Remove the average trace \
                                 (average between hfilt0 and hfilt1)')
    parser_proc.add_argument('-ahfilt',
                             nargs=1,
                             type=int,
                             help='Adaptive horizontal filtering')
    parser_proc.add_argument('-rev',
                             action='store_true',
                             help='Reverse profile')
    parser_proc.add_argument('-nmo',
                             nargs=2,
                             type=float,
                             help='Normal moveout correction. \
                                     First argument is the \
                                     transmitter-receiver separation. \
                                     Second argument is the velocity \
                                     of the radar wave (in m/s).')
    parser_proc.add_argument('-crop',
                             nargs=3,
                             type=str,
                             help='Crop the radar data in the travel-time \
                                    direction. Args are the limit, whether \
                                    to crop off ["top", "bottom"], with limit \
                                    defined in terms of \
                                    ["snum", "twtt", "depth"]')
    parser_proc.add_argument('-hcrop',
                             nargs=3,
                             type=str,
                             help='Crop the radar data in the horizontal. \
                                     Arguments are the limit, whether to crop \
                                     off ["left", "right], with limit defined \
                                     in terms of ["tnum", "dist"]')
    parser_proc.add_argument('-restack',
                             nargs=1,
                             type=int,
                             help='Restack to this (odd) number of traces')
    parser_proc.add_argument('-interp',
                             nargs=2,
                             type=str,
                             help='Reinterpolate GPS. \
                                     First argument is the new spacing, in \
                                     meters. Second argument is the filename \
                                     (csv or mat) with the new GPS data')
    parser_proc.add_argument('-denoise',
                             nargs=2,
                             type=int,
                             help='Denoising filter vertical and horizontal (scipy wiener for now)')
    parser_proc.add_argument('-migrate',
                             type=str,
                             help='Migrate with the indicated routine.')
    parser_proc.add_argument('fn',
                             type=str,
                             nargs='+',
                             help='File(s) to process')
    parser_proc.add_argument('-o', type=str, help='Write to this filename')

    # plotting
    parser_plot = subparsers.add_parser('plot', help='Plot data')
    parser_plot.set_defaults(func=plot.plot)
    parser_plot.add_argument('fns',
                             type=str,
                             nargs='+',
                             help='File(s) to plot')
    parser_plot.add_argument('-s',
                             action='store_true',
                             help='Save file (do not plt.show())')
    parser_plot.add_argument('-yd', action='store_true',
                             help='Plot the depth rather than travel time')
    parser_plot.add_argument('-xd', action='store_true',
                             help='Plot the dist rather than the trace num')
    parser_plot.add_argument('-tr', nargs=2, type=int, default=None,
                             help='Plot the traces in this range (line plot)')
    parser_plot.add_argument('-power', type=int, default=None, help='Input a picked layer number to plot the RMS power for each trace in map view.')
    parser_plot.add_argument('-spectra', nargs=2, type=float, default=None,
                             help='Plot power spectral density across traces of radar profile. Input frequency bounds (MHz).')
    parser_plot.add_argument('-o', type=str, help='Write to this filename')
    parser_plot.add_argument('-freq_limit',
                             type=float,
                             default=None,
                             help='Maximum frequeny to plot power spectral \
                                     density to')
    parser_plot.add_argument('-window',
                             type=str,
                             default='hanning',
                             help='Type of window function to be used for the singal.periodogram() method')
    parser_plot.add_argument('-scaling',
                             type=str,
                             default='spectrum',
                             help='Whether to plot power spectral density or power spectrum: default is spectrum')

    parser_convert = subparsers.add_parser('convert',
                                           help='Convert filetype (lossy)')
    parser_convert.set_defaults(func=convert.convert)
    parser_convert.add_argument('fns_in',
                                type=str,
                                nargs='+',
                                help='File(s) to convert')
    parser_convert.add_argument('out_fmt',
                                type=str,
                                choices=['shp', 'mat', 'segy'])
    parser_convert.add_argument('-in_fmt',
                                type=str,
                                default=None,
                                choices=load.FILETYPE_OPTIONS,
                                help='Input format type. If none, guess from extension, but  be warned, we are bad at guessing!')
    parser_convert.add_argument('-t_srs', type=str, default=None,
                                help='Target srs, in a format recognized by gdal. Default None (write raw input)')
    return parser


def main():
    """Call impdar exec."""
    parser = _get_args()
    args = parser.parse_args(sys.argv[1:])
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])
        return None
    return args.func(**vars(args))


if __name__ == '__main__':
    main()
