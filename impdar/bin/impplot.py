#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
#
# Legacy header: Created: B. Welch, modified by S. Harris, J. Olson, and B. Youngblood.

import sys
import argparse
from impdar.lib import plot
from impdar.lib.load import FILETYPE_OPTIONS


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    add_simple_procparser(subparsers, 'rg', 'Plot radargram', plot_radargram, defname='radargram', xd=True, yd=True)

    trace_parser = add_simple_procparser(subparsers, 'traces', 'Plot traces vs depth', plot_traces, defname='traces', xd=False, yd=True)
    trace_parser.add_argument('t_start', type=int, help='Starting trace number')
    trace_parser.add_argument('t_end', type=int, help='Ending trace number')

    power_parser = add_simple_procparser(subparsers, 'power', 'Plot power on a layer', plot_power, defname='power', xd=False, yd=False, other_ftypes=False)
    power_parser.add_argument('layer', type=int, help='Layer upon which to plot the power')
    return parser


def add_simple_procparser(subparsers, name, helpstr, func, defname='proc', xd=False, yd=False, other_ftypes=True):
    parser = add_procparser(subparsers, name, helpstr, func, defname=defname)
    add_def_args(parser, xd=xd, yd=yd)
    return parser


def add_procparser(subparsers, name, helpstr, func, defname='proc'):
    parser = subparsers.add_parser(name, help=helpstr)
    parser.set_defaults(func=func, name=defname)
    return parser


def add_def_args(parser, xd=False, yd=False, other_ftypes=True):
    parser.add_argument('fns', type=str, nargs='+', help='The files to process')
    parser.add_argument('-o', type=str, help='Output to this file (or folder if multiple inputs)')

    parser.add_argument('-s', action='store_true', help='Save file (do not plt.show())')
    parser.add_argument('--o_fmt', type=str, default='png', help='Save file with this extension (default png)')
    parser.add_argument('-dpi', type=int, default=300, help='Save file with this resolution (default 300)')

    if xd:
        parser.add_argument('-xd', action='store_true', help='Plot the distance rather than the trace number')
    if yd:
        parser.add_argument('-yd', action='store_true', help='Plot the depth rather than travel time')

    if other_ftypes:
        parser.add_argument('--in_fmt', type=str,
                            help='Type of file',
                            default='mat',
                            choices=FILETYPE_OPTIONS)


def plot_radargram(fns=None, s=False, o=None, xd=False, yd=False, o_fmt='png', dpi=300, in_fmt='mat', **kwargs):
    plot.plot(fns, xd=xd, yd=yd, s=s, o=o, ftype=o_fmt, dpi=dpi, filetype=in_fmt)


def plot_power(fns=None, layer=None, s=False, o=None, o_fmt='png', dpi=300, in_fmt='mat', **kwargs):
    plot.plot(fns, power=layer, s=s, o=o, ftype=o_fmt, dpi=dpi, filetype=in_fmt)


def plot_traces(fns=None, t_start=None, t_end=None, yd=False, s=False, o=None, o_fmt='png', dpi=300, in_fmt='mat', **kwargs):
    plot.plot(fns, tr=(t_start, t_end), yd=yd, s=s, o=o, ftype=o_fmt, dpi=dpi, filetype=in_fmt)


def main():
    parser = _get_args()
    args = parser.parse_args(sys.argv[1:])
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])
        return
    args.func(**vars(args))


if __name__ == '__main__':
    main()
