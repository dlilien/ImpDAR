#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien@berens>
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
from impdar import load, process, plot


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_load = subparsers.add_parser('load', help='Load data')
    parser_load.set_defaults(func=load.load_and_exit)
    parser_load.add_argument('filetype', type=str, help='Type of file', choices=['gssi', 'pe', 'mat'])
    parser_load.add_argument('fn', type=str, nargs='+', help='File(s) to load')
    parser_load.add_argument('-o', type=str, help='Write to this filename')

    # Options for processing data
    parser_proc = subparsers.add_parser('proc', help='Process data')
    parser_proc.set_defaults(func=process.process_and_exit)
    parser_proc.add_argument('-gssi', action='store_true', help='Indicates that the file(s) are gssi output')
    parser_proc.add_argument('-pe', action='store_true', help='Indicates that the file(s) are pulse ekko output')
    parser_proc.add_argument('-cat', action='store_true', help='Concatenate the files')
    parser_proc.add_argument('-vbp', nargs=2, type=float, help='Bandpass the data vertically at low (MHz) and high (MHz)')
    parser_proc.add_argument('-hfilt', nargs=2, type=int, help='Remove the average trace (average between hfilt0 and hfilt1)')
    parser_proc.add_argument('-ahfilt', action='store_true', help='Adaptive horizontal filtering')
    parser_proc.add_argument('-rev', action='store_true', help='Reverse profile')
    parser_proc.add_argument('-nmo', nargs=2, type=float, help='Normal moveout correction. First argument is the transmitter-receiver separation. Second argument is the velocity of the radar wave (in m/s).')
    parser_proc.add_argument('-crop', nargs=3, type=str, help='Crop the radar data in the travel-time direction. Arguments are to crop off ["top", "bottom"], with limit defined in terms of ["snum", "twtt", "depth"] and then the limit itself')
    parser_proc.add_argument('-restack', nargs=1, type=int, help='Restack to this (odd) number of traces')
    parser_proc.add_argument('-interp', nargs=2, type=str, help='Reinterpolate GPS. First argument is the new spacing, in meters. Second argument is the filename (csv or mat) with the new GPS data')
    parser_proc.add_argument('fn', type=str, nargs='+', help='File(s) to process')
    parser_proc.add_argument('-o', type=str, help='Write to this filename')

    parser_plot = subparsers.add_parser('plot', help='Plot data')
    parser_plot.set_defaults(func=plot.plot)
    parser_plot.add_argument('fn', type=str, nargs='+', help='File(s) to plot')
    parser_plot.add_argument('-s', action='store_true', help='Save file (do not plt.show())')
    parser_plot.add_argument('-yd', action='store_true', help='Plot the depth rather than travel time')
    parser_plot.add_argument('-xd', action='store_true', help='Plot the distance rather than the trace number')
    parser_plot.add_argument('-o', type=str, help='Write to this filename')
    return parser


def main():
    parser = _get_args()
    args = parser.parse_args()
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])
        return
    return args.func(**vars(args))


if __name__ == '__main__':
    main()
