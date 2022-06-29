#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Benjamin Hills <bhills@uw.edu>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make an executable for single actions of QuadPol apres handling.
"""

import sys
import os.path
import argparse
import numpy as np

from impdar.lib.ApresData import FILETYPE_OPTIONS
from impdar.lib.ApresData import load_apres, QuadPol

from impdar.lib import plot

def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Choose a processing step')

    # Loading functionality
    parser_load = _add_procparser(subparsers,
                                  'load',
                                  'load quadpol data',
                                  lambda x: x,
                                  defname='load')
    _add_def_args(parser_load)

    # Initial range conversion (pulse compression)
    parser_fullproc = _add_procparser(subparsers,
                                      'proc',
                                      'full processing flow on the quadpol data object',
                                      full_processing,
                                      'proc')
    parser_fullproc.add_argument('-nthetas',
                                 type=int,
                                 help='number of theta values to rotate into')
    parser_fullproc.add_argument('-dtheta',
                                 type=int,
                                 help='size of coherence window in theta direction')
    parser_fullproc.add_argument('-drange',
                                 type=int,
                                 help='size of coherence window in range direction')
    _add_def_args(parser_fullproc)

    # Initial range conversion (pulse compression)
    parser_rotate = _add_procparser(subparsers,
                                    'rotated',
                                    'use a rotational transform to find data at all azimuths',
                                    rotate,
                                    'rotated')
    parser_rotate.add_argument('-nthetas',
                               type=int,
                               help='number of theta values to rotate into')
    _add_def_args(parser_rotate)

    # Coherence 2D
    parser_coherence = _add_procparser(subparsers,
                                       'chhvv',
                                       '2-dimensional coherence between HH and VV polarizations',
                                       coherence,
                                       'chhvv')
    parser_coherence.add_argument('-dtheta',
                                  type=int,
                                  help='size of coherence window in theta direction')
    parser_coherence.add_argument('-drange',
                                  type=int,
                                  help='size of coherence window in range direction')
    _add_def_args(parser_coherence)

    # Cross-Polarized Extinction
    parser_cpe = _add_procparser(subparsers,
                                 'cpe',
                                 'find the depth-azimuth profile for cross-polarized extinction',
                                 cross_polarized_extinction,
                                 'cpe')
    parser_cpe.add_argument('-Wn',
                            type=float,
                            help='filter frequency')
    parser_cpe.add_argument('-fs',
                            type=float,
                            help='sampling frequency')
    _add_def_args(parser_cpe)

    # Plot Quad Pol Apres Acquisition
    parser_plot = _add_procparser(subparsers,
                                  'plot',
                                  'plot apres data from a single acquisition',
                                  plot_quadpol)
    parser_plot.add_argument('-s',
                             action='store_true',
                             help='save file (do not plt.show())')
    parser_plot.add_argument('-yd', action='store_true',
                             help='plot the depth rather than travel time')
    _add_def_args(parser_plot)

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

    try:
        qp_data = load_apres.load_quadpol(args.fns)
    except:
        qp_data = QuadPol(args.fns[0])


    if args.name == 'load':
        name = 'raw'
        pass
    else:
        name = args.name
        args.func(apres_data, **vars(args))

    if args.o is not None:
        out_fn = args.o
        qp_data.save(out_fn)
    else:
        bn = os.path.splitext(args.fns[0])[0]
        if bn[-4:] == '_raw':
            bn = bn[:-4]
        out_fn = bn + '_{:s}.mat'.format(name)
        qp_data.save(out_fn)


def full_processing(dat, nthetas=100, dtheta=20.0*np.pi/180.,
                    drange=100., Wn=0., fs=0., **kwargs):
    """Full processing flow for QuadPol object.
    Range conversion, stacking, uncertainty."""
    dat.rotational_transform(n_thetas=nthetas, **kwargs)
    dat.coherence2d(delta_theta=dtheta, delta_range=drange)
    dat.power_anomaly()
    dat.lowpass(Wn=Wn, fs=fs)


def rotate(dat, nthetas=100, **kwargs):
    """Range conversion."""
    dat.rotational_transform(n_thetas=nthetas, **kwargs)


def coherence(dat, dtheta=20.0*np.pi/180., drange=100., **kwargs):
    """Stack chirps."""
    dat.coherence2d(delta_theta=dtheta, delta_range=drange)


def cross_polarized_extinction(dat,
                               Wn=0., fs=0., **kwargs):
    """Calculate uncertainty."""
    dat.power_anomaly()
    dat.lowpass(Wn=Wn, fs=fs)


def plot_quadpol(dat, s=False, o=None, o_fmt='png',
                 dpi=300, **kwargs):
    plot.plot_quadpol(dat, s=s, o=o, ftype=o_fmt, dpi=dpi)


if __name__ == '__main__':
    main()
