#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make an executable for single actions of QuadPole pRES handling.
"""
import sys
import os.path
import argparse

import numpy as np

from impdar.lib.ApresData import load_quadpol


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Choose a processing step')

    # Load and exit
    parser_load = _add_procparser(subparsers,
                                  'load',
                                  'load apres data',
                                  lambda x: x,
                                  defname='load')
    _add_def_args(parser_load)

    # Rotational Transform
    parser_rotate = _add_procparser(subparsers,
                                   'rotate',
                                   'Rotate to full azimuthal dependency',
                                   rotate,
                                   defname='rotated')
    parser_rotate.add_argument('--theta0',
                              default=0.0,
                              help='Minimum azimuth (in radians)',
                              type=float)
    parser_rotate.add_argument('--theta1',
                              default=np.pi,
                              help='Maximum azimuth (in radians)',
                              type=float)
    parser_rotate.add_argument('--n',
                              help='Number of columns (azimuths) in matrix',
                              default=100,
                              type=int)
    _add_def_args(parser_rotate)

    # Calculation of the HHVV coherence
    parser_coherence = _add_procparser(subparsers,
                                       'coherence',
                                       'Calculate HHVV coherence',
                                       coherence,
                                       defname='coherence')
    parser_coherence.add_argument('--dtheta',
                              default=20.0 * np.pi / 180.0,
                              help='Window size in azimuth, in radians. Default equivalent to 20 degrees.',
                              type=float)
    parser_coherence.add_argument('--drange',
                              default=100.,
                              help='Window in the vertical, in meters. Default 100.',
                              type=float)
    _add_def_args(parser_coherence)
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


def rotate(dat, theta0=0, theta1=np.pi, n=100, **kwargs):
    return dat.rotational_transform(theta_start=theta0, theta_end=theta1, n_thetas=n)


def coherence(dat, dtheta=0.1, drange=100.0, **kwargs):
    return dat.coherence2d(delta_theta=dtheta, delta_range=drange)


def main():
    """Get arguments, process data, handle saving."""
    parser = _get_args()
    args = parser.parse_args(sys.argv[1:])
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])

    if args.name == 'load':
        apres_data = [load_quadpol.load_quadpol(fn) for fn in args.fns]
        pass
    else:
        apres_data = [load_quadpol.load_quadpol(fn, load_single_pol=False) for fn in args.fns]
        for dat in apres_data:
            args.func(dat, **vars(args))

    if args.name == 'load':
        name = 'qp'
    else:
        name = args.name

    if args.o is not None:
        if ((len(apres_data) > 1) or (args.o[-1] == '/')):
            for d, f in zip(apres_data, args.fns):
                bn = os.path.split(os.path.splitext(f)[0])[1]
                if bn[-3:] == '_qp':
                    bn = bn[:-3]
                out_fn = os.path.join(args.o, bn + '_{:s}.h5'.format(name))
                d.save(out_fn)
        else:
            out_fn = args.o
            apres_data[0].save(out_fn)
    else:
        for d, f in zip(apres_data, args.fns):
            bn = os.path.splitext(f)[0]
            if bn[-4:] == '_qp':
                bn = bn[:-4]
            out_fn = bn + '_{:s}.h5'.format(name)
            d.save(out_fn)



if __name__ == '__main__':
    main()
