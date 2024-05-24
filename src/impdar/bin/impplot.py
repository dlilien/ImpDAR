#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
#
"""The impdar plotter."""

import sys
import argparse
from impdar.lib import plot
from impdar.lib.load import FILETYPE_OPTIONS


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    rg_parser = _add_simple_procparser(subparsers,
                                       'rg',
                                       'Plot radargram',
                                       plot_radargram,
                                       defname='radargram',
                                       xd=True,
                                       yd=True,
                                       dualy=True)
    rg_parser.add_argument('-picks',
                           action='store_true',
                           help='Plot picks')
    rg_parser.add_argument('-clims',
                           nargs=2,
                           type=float,
                           help='Color limits')
    rg_parser.add_argument('-flatten_layer',
                           type=int,
                           default=None,
                           help='Distort plot so this layer is flat')
    rg_parser.add_argument('-cmap',
                           type=str,
                           default='gray',
                           help='Color map name')

    _add_simple_procparser(subparsers,
                           'ft',
                           'Plot ft',
                           plot_ft,
                           defname='spec')

    _add_simple_procparser(subparsers,
                           'hft',
                           'Plot ft',
                           plot_hft,
                           defname='spec')

    trace_parser = _add_simple_procparser(subparsers,
                                          'traces',
                                          'Plot traces vs depth',
                                          plot_traces,
                                          defname='traces',
                                          xd=False,
                                          yd=True,
                                          dualy=True)
    trace_parser.add_argument('t_start',
                              type=int,
                              help='Starting trace number')
    trace_parser.add_argument('t_end',
                              type=int,
                              help='Ending trace number')

    power_parser = _add_simple_procparser(subparsers,
                                          'power',
                                          'Plot power on a layer',
                                          plot_power,
                                          defname='power',
                                          other_ftypes=False)
    power_parser.add_argument('layer',
                              type=int,
                              help='Layer upon which to plot the power')

    spec_parser = _add_simple_procparser(subparsers,
                                         'spectrogram',
                                         'Plot spectrogram for all traces',
                                         plot_spectrogram,
                                         defname='spectrogram',
                                         other_ftypes=False)
    spec_parser.add_argument('freq_lower',
                             type=float,
                             help='Lower frequency bound')
    spec_parser.add_argument('freq_upper',
                             type=float,
                             help='Uppwer frequency bound')
    return parser


def _add_simple_procparser(subparsers, name, helpstr, func, defname='proc',
                           xd=False, yd=False, dualy=False, other_ftypes=True):
    """Add a simple subparser."""
    parser = _add_procparser(subparsers, name, helpstr, func, defname=defname)
    _add_def_args(parser, xd=xd, yd=yd, dualy=dualy)
    return parser


def _add_procparser(subparsers, name, helpstr, func, defname='proc'):
    """Add a subparser."""
    parser = subparsers.add_parser(name, help=helpstr)
    parser.set_defaults(func=func, name=defname)
    return parser


def _add_def_args(parser, xd=False, yd=False, dualy=False, other_ftypes=True):
    """Set up common arguments for different types of commands."""
    parser.add_argument('fns',
                        type=str,
                        nargs='+',
                        help='The files to process')
    parser.add_argument('-o',
                        type=str,
                        help='Output to this file (folder if multiple inputs)')

    parser.add_argument('-s',
                        action='store_true',
                        help='Save file (do not plt.show())')
    parser.add_argument('--o_fmt',
                        type=str,
                        default='png',
                        help='Save file with this extension (default png)')
    parser.add_argument('-dpi',
                        type=int,
                        default=300,
                        help='Save file with this resolution (default 300)')

    if xd:
        parser.add_argument('-xd',
                            action='store_true',
                            help='Plot the dist rather than the trace number')
    if yd:
        parser.add_argument('-yd',
                            action='store_true',
                            help='Plot the depth rather than travel time')

    if dualy:
        parser.add_argument('-dualy',
                            action='store_true',
                            help='Primary y axis is TWTT, secondary is depth')

    if other_ftypes:
        parser.add_argument('--in_fmt', type=str,
                            help='Type of file',
                            default='mat',
                            choices=FILETYPE_OPTIONS)


def plot_radargram(fns=None, s=False, o=None, xd=False, yd=False, o_fmt='png',
                   dpi=300, in_fmt='mat', picks=False, clims=None, cmap='gray',
                   flatten_layer=None, dualy=False, **kwargs):
    """Plot data as a radio echogram."""
    plot.plot(fns, xd=xd, yd=yd, s=s, o=o, ftype=o_fmt, dpi=dpi,
              filetype=in_fmt, pick_colors=picks, cmap=cmap, clims=clims,
              flatten_layer=flatten_layer, dualy=dualy)


def plot_ft(fns=None, s=False, o=None, xd=False, yd=False, o_fmt='png',
            dpi=300, in_fmt='mat', **kwargs):
    """
    Plot the fourier spectrum of the data.

    Can be useful if you have mystery data of unknown frequency.
    """
    plot.plot(fns, xd=xd, yd=yd, s=s, o=o, ftype=o_fmt, dpi=dpi,
              filetype=in_fmt, ft=True)


def plot_hft(fns=None, s=False, o=None, xd=False, yd=False, o_fmt='png',
             dpi=300, in_fmt='mat', **kwargs):
    """
    Plot the fourier spectrum of the data in the horizontal.

    Might be useful for guessing how to horizontally filter.
    """
    plot.plot(fns, xd=xd, yd=yd, s=s, o=o, ftype=o_fmt, dpi=dpi,
              filetype=in_fmt, hft=True)


def plot_power(fns=None, layer=None, s=False, o=None, o_fmt='png',
               dpi=300, in_fmt='mat', **kwargs):
    """Plot the return power of a particular layer."""
    plot.plot(fns, power=layer, s=s, o=o, ftype=o_fmt, dpi=dpi,
              filetype=in_fmt)


def plot_traces(fns=None, t_start=None, t_end=None, yd=False, dualy=False,
                s=False, o=None, o_fmt='png', dpi=300, in_fmt='mat', **kwargs):
    """Plot traces in terms of amplitude vs some vertical variable."""
    plot.plot(fns, tr=(t_start, t_end), yd=yd, s=s, o=o, ftype=o_fmt, dpi=dpi,
              dualy=dualy, filetype=in_fmt)


def plot_spectrogram(fns=None, freq_lower=None, freq_upper=None, window=None,
                     scaling='spectrum', yd=False, s=False, o=None,
                     o_fmt='png', dpi=300, in_fmt='mat', **kwargs):
    """Plot a spectrogram."""
    plot.plot(fns,
              spectra=(freq_lower, freq_upper),
              window=window,
              scaling=scaling,
              yd=yd,
              s=s,
              o=o,
              ftype=o_fmt,
              dpi=dpi,
              filetype=in_fmt)


def main():
    """Get arguments, plot data."""
    parser = _get_args()
    args = parser.parse_args(sys.argv[1:])
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])
        return
    args.func(**vars(args))


if __name__ == '__main__':
    main()
