#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make an executable for single actions of ApRES handling.
"""

import sys
import os.path
import argparse

from impdar.lib.ApresData import load_apres, load_diff
from impdar.lib.ApresData import FILETYPE_OPTIONS, ApresDiffData

from impdar.lib import plot

def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Choose a processing step')

    # Loading functionality
    parser_load = _add_procparser(subparsers,
                                  'load',
                                  'load apres data',
                                  lambda x: x,
                                  defname='load')
    _add_def_args(parser_load)

    # Initial range conversion (pulse compression)
    parser_fullproc = _add_procparser(subparsers,
                                      'proc',
                                      'full processing flow on the apres data object',
                                      full_processing,
                                      'proc')
    parser_fullproc.add_argument('-max_range',
                             type=int,
                             help='maximum range for the pulse compression')
    parser_fullproc.add_argument('-num_chirps',
                             type=int,
                              help='number of chirps to stack (default: stack all)')
    parser_fullproc.add_argument('noise_bed_range',
                             type=int,
                              help='bed range under which \
                                    the noise phasor will be calculated')
    _add_def_args(parser_fullproc)

    # Initial range conversion (pulse compression)
    parser_range = _add_procparser(subparsers,
                                  'range',
                                  'convert the recieved waveform to a \
                                        range-amplitude array',
                                  pulse_compression,
                                  'range')
    parser_range.add_argument('-max_range',
                             type=int,
                             help='maximum range for the pulse compression')
    _add_def_args(parser_range)

    # Stacking
    parser_stack = _add_procparser(subparsers,
                                  'stack',
                                  'stack apres chirps into a single array',
                                  stack,
                                  'stacked')
    parser_stack.add_argument('-num_chirps',
                             type=int,
                              help='number of chirps to stack (default: stack all)')
    _add_def_args(parser_stack)

    # Uncertainty
    parser_unc = _add_procparser(subparsers,
                                  'uncertainty',
                                  'calculate the phase uncertainty for all samples',
                                  uncertainty,
                                  'uncertainty')
    parser_unc.add_argument('-noise_bed_range',
                             type=int,
                              help='bed range under which \
                                    the noise phasor will be calculated')
    _add_def_args(parser_unc)

    # Load Differencing Object from two impdar acquisitions
    parser_diffload = _add_procparser(subparsers,
                                    'diffload',
                                    'create an ApresDiff object',
                                    lambda x: x,
                                    defname='diffload')
    _add_def_args(parser_diffload)

    # Full Difference Processing
    parser_diffproc = _add_procparser(subparsers,
                                    'diffproc',
                                    'create an ApresDiff object and then execuate \
                                        the full differencing processing flow',
                                    full_differencing,
                                    'diffproc')
    parser_diffproc.add_argument('-window',
                             type=int,
                              help='window size over which the cross correlation is done')
    parser_diffproc.add_argument('-step',
                             type=int,
                              help='step size in samples for the moving window')
    parser_diffproc.add_argument('-thresh',
                             type=int,
                              help='step size in samples for the moving window')
    parser_diffproc.add_argument('-strain_window',
                             type=tuple,
                              help='step size in samples for the moving window')
    parser_diffproc.add_argument('-w_surf',
                             type=float,
                              help='surface vertical velocity (ice equivalent accumulation rate)')
    parser_diffproc.set_defaults(window=20, step=20,thresh=0.95,strain_window=(200,1000),w_surf=-0.15)
    _add_def_args(parser_diffproc)

    # Phase Differencing
    parser_pdiff = _add_procparser(subparsers,
                                  'pdiff',
                                  'unwrap the differenced phase profile \
                                       from top to bottom',
                                  phase_differencing)
    parser_pdiff.add_argument('-window',
                             type=int,
                              help='window size over which the cross correlation is done')
    parser_pdiff.add_argument('-step',
                             type=int,
                              help='step size in samples for the moving window')
    parser_pdiff.set_defaults(window=20, step=20)
    _add_def_args(parser_pdiff)

    # Phase Unwrap
    parser_unwrap = _add_procparser(subparsers,
                                  'unwrap',
                                  'unwrap the differenced phase profile \
                                       from top to bottom',
                                  unwrap)
    _add_def_args(parser_unwrap)

    # Range Differencing
    parser_rdiff = _add_procparser(subparsers,
                                  'rdiff',
                                  'convert the differenced phase profile to range',
                                  range_differencing)
    _add_def_args(parser_rdiff)

    # Plot Single Apres Acquisition
    parser_plot = _add_procparser(subparsers,
                                  'plot',
                                  'plot apres data from a single acquisition',
                                  plot_single)
    parser_plot.add_argument('-s',
                             action='store_true',
                             help='save file (do not plt.show())')
    parser_plot.add_argument('-yd', action='store_true',
                             help='plot the depth rather than travel time')
    _add_def_args(parser_plot)

    # Plot Differenced Apres Acquisitions
    parser_plotdiff = _add_procparser(subparsers,
                                  'plot_diff',
                                  'plot profiles for differenced apres acquisitions',
                                  plot_differenced)
    parser_plotdiff.add_argument('-s',
                             action='store_true',
                             help='Save file (do not plt.show())')
    parser_plotdiff.add_argument('-yd', action='store_true',
                             help='Plot the depth rather than travel time')
    _add_def_args(parser_plotdiff)

    return parser


def _add_simple_procparser(subparsers, name, helpstr, func, defname='proc'):
    """Add a sub parser that can do a simple thing with no arguments."""
    parser = _add_procparser(subparsers, name, helpstr, func, defname=defname)
    _add_def_args(parser)
    return parser


def _add_procparser(subparsers, name, helpstr, func, defname='proc'):
    """Wrap adding subparser because we mostly want the same args."""
    parser = subparsers.add_parser(name, help=helpstr)
    parser.set_defaults(func=func, name=defname)
    return parser


def _add_def_args(parser):
    """Set some default arguments common to the different processing types."""
    parser.add_argument('fns',
                        type=str,
                        nargs='+',
                        help='The files to process')
    parser.add_argument('-o',
                        type=str,
                        help='Output to this file (folder if multiple inputs)')


def main():
    """Get arguments, process data, handle saving."""
    parser = _get_args()
    args = parser.parse_args(sys.argv[1:])
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])

    if args.name == 'diffload':
        apres_data = load_diff.load_diff(args.fns[0],args.fns[1])
    else:
        try:
            apres_data = load_apres.load_apres(args.fns)
        except:
            apres_data = ApresDiffData(args.fns[0])


    if args.name == 'load':
        name = 'raw'
        pass
    else:
        name = args.name
        args.func(apres_data, **vars(args))

    if args.o is not None:
        out_fn = args.o
        apres_data.save(out_fn)
    else:
        bn = os.path.splitext(args.fns[0])[0]
        if bn[-4:] == '_raw':
            bn = bn[:-4]
        out_fn = bn + '_{:s}.mat'.format(name)
        apres_data.save(out_fn)


def full_processing(dat, p=2, max_range=4000, num_chirps=0, noise_bed_range=3000, **kwargs):
    """Full processing flow for ApresData object.
    Range conversion, stacking, uncertainty."""
    dat.apres_range(p,max_range)
    if num_chirps == 0:
        dat.stacking()
    else:
        dat.stacking(num_chirps)
    dat.phase_uncertainty(noise_bed_range)


def pulse_compression(dat, p=2, max_range=4000, **kwargs):
    """Range conversion."""
    dat.apres_range(p,max_range)


def stack(dat, num_chirps=0, **kwargs):
    """Stack chirps."""
    if num_chirps == 0:
        dat.stacking()
    else:
        dat.stacking(num_chirps)


def uncertainty(dat,noise_bed_range=3000, **kwargs):
    """Calculate uncertainty."""
    dat.phase_uncertainty(noise_bed_range)


def full_differencing(diffdat, win=20, step=20, thresh=0.95,
                      strain_window=(200,1000), w_surf=-0.15, **kwargs):
    diffdat.phase_diff(win,step)
    diffdat.phase_unwrap(win,thresh)
    diffdat.range_diff()
    diffdat.strain_rate(strain_window=strain_window,w_surf=w_surf)
    diffdat.bed_pick()


def phase_differencing(diffdat, win=20, step=20, **kwargs):
    diffdat.phase_diff(win,step)


def unwrap(diffdat,win=20,thresh=.95, **kwargs):
    diffdat.phase_unwrap(win,thresh)


def range_differencing(diffdat, **kwargs):
    diffdat.range_diff()


def plot_single(dat, s=False, o=None, o_fmt='png',
               dpi=300, **kwargs):
    """Plot the return power of a particular layer."""
    plot.plot_apres(dat, s=s, o=o, ftype=o_fmt, dpi=dpi)

def plot_differenced(dat, s=False, o=None, o_fmt='png',
                     dpi=300, **kwargs):
    plot.plot_apres_diff(dat, s=s, o=o, ftype=o_fmt, dpi=dpi)


if __name__ == '__main__':
    main()
