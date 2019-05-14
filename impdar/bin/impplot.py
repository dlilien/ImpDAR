#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
#
# Legacy header:
#	Created: B. Welch - 10/15/01
#	Modification History:
#       1)  Added ability to load processed in Stodeep - S. Harris 6/5/02
# 		2)	Converted to new structure-based flagging format - J. Olson 7/10/08
#       3)  Added call for new batchdeep.m shell - B. Youngblood 7/12/08

import argparse
from matplotlib.pyplot import show
from impdar import plot


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    radargram_parser = add_simple_procparser(subparsers, 'rg', 'Plot radargram', plot_radargram, defname='radargram', xd=True, yd=True)

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
    parser.add_argument('-ftype', type=str, default='png', help='Save file with this extension (default png)')
    parser.add_argument('-dpi', type=int, default=300, help='Save file with this resolution (default 300)')

    if xd:
        parser.add_argument('-xd', action='store_true', help='Plot the distance rather than the trace number')
    if yd:
        parser.add_argument('-yd', action='store_true', help='Plot the depth rather than travel time')

    if other_ftypes:
        parser.add_argument('-pe', action='store_true', help='Inputs are pulse ekko files')
        parser.add_argument('-gssi', action='store_true', help='Inputs are gssi files')
        parser.add_argument('-gprMax', action='store_true', help='Inputs are gprMax files')
        parser.add_argument('-gecko', action='store_true', help='Inputs are gecko files')
        parser.add_argument('-segy', action='store_true', help='Inputs are segy files')


def plot_radargram(fns=None, s=False, o=None, xd=False, yd=False, ftype='png', dpi=300, gssi=False, pe=False, gprMax=False, gecko=False, segy=False, **kwargs):
    plot.plot(fns, xd=xd, yd=yd, s=s, o=o, ftype=ftype, dpi=dpi, gssi=gssi, pe=pe, gprMax=gprMax, gecko=gecko, segy=segy)


def plot_power(fns=None, layer=None, s=False, o=None, ftype='png', dpi=300, **kwargs):
    plot.plot(fns, power=layer, s=s, o=o, ftype=ftype, dpi=dpi)


def plot_traces(fns=None, t_start=None, t_end=None, yd=False, s=False, o=None, ftype='png', dpi=300, gssi=False, pe=False, gprMax=False, gecko=False, segy=False, **kwargs):
    plot.plot(fns, tr=(t_start, t_end), yd=yd, s=s, o=o, ftype=ftype, dpi=dpi, gssi=gssi, pe=pe, gprMax=gprMax, gecko=gecko, segy=segy)


def main():
    parser = _get_args()
    args = parser.parse_args()
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])
        return
    args.func(**vars(args))


if __name__ == '__main__':
    main()
