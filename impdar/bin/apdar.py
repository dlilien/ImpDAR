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
import numpy as np

from impdar.lib.ApresData import load_apres, load_time_diff, load_quadpol
from impdar.lib.ApresData import FILETYPE_OPTIONS, ApresData, ApresTimeDiff, ApresQuadPol
from impdar.lib.ApresData.ApresFlags import ApresFlags, TimeDiffFlags, QuadPolFlags

from impdar.lib import plot


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Choose a processing step')

    # Loading functionality
    parser_load = _add_procparser(subparsers,
                                  'load',
                                  'load apres data',
                                  load,
                                  defname='load')
    parser_load.add_argument('-acq_type',
                                 type=str,
                                 help='Acquisition type',
                                 default='single',
                                 choices=['single', 'timediff', 'quadpol'])
    _add_def_args(parser_load)

    # Full processing flow for a single ApRES acquisition
    parser_singleproc = _add_procparser(subparsers,
                                      'proc',
                                      'full processing flow on the apres data object',
                                      single_processing,
                                      'proc')
    parser_singleproc.add_argument('-max_range',
                                 type=float,
                                 help='maximum range for the range conversion')
    parser_singleproc.add_argument('-num_chirps',
                                 type=int,
                                 help='number of chirps to stack (default: stack all)')
    parser_singleproc.add_argument('-noise_bed_range',
                                 type=float,
                                 help='bed range under which the noise phasor will be calculated')
    parser_singleproc.set_defaults(max_range=4000., num_chirps=0, noise_bed_range=3000.)
    _add_def_args(parser_singleproc)

    # Full Difference Processing
    parser_diffproc = _add_procparser(subparsers,
                                    'diffproc',
                                    'create an ApresDiff object and then execute the full differencing processing flow',
                                    time_diff_processing,
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

    # Full quadpol processing
    parser_qpproc = _add_procparser(subparsers,
                                      'qpproc',
                                      'full processing flow on the quadpol data object',
                                      quadpol_processing,
                                      'qpproc')
    parser_qpproc.add_argument('-nthetas',
                                 type=int,
                                 help='number of theta values to rotate into')
    parser_qpproc.add_argument('-dtheta',
                                 type=float,
                                 help='size of coherence window in theta direction')
    parser_qpproc.add_argument('-drange',
                                 type=float,
                                 help='size of coherence window in range direction')
    parser_qpproc.add_argument('-cross_pol_flip',
                                 type=str,
                                 help='flip the sign on one of the cross polarized terms.')
    parser_qpproc.set_defaults(nthetas=100, dtheta=20.*np.pi/180., drange=100, cross_pol_flip=False)
    _add_def_args(parser_qpproc)

    # Initial range conversion (deramping)
    parser_range = _add_procparser(subparsers,
                                  'range',
                                  'convert the recieved waveform to a range-amplitude array',
                                  range_conversion,
                                  'range')
    parser_range.add_argument('-max_range',
                             type=float,
                             help='maximum range for the range conversion',
                             default=4000.)
    _add_def_args(parser_range)

    # Stacking
    parser_stack = _add_procparser(subparsers,
                                   'stack',
                                   'stack apres chirps into a single array',
                                   stack,
                                   'stacked')
    parser_stack.add_argument('-num_chirps',
                              type=int,
                              help='number of chirps to stack (default: stack all)',
                              default=0)
    _add_def_args(parser_stack)

    # Uncertainty
    parser_unc = _add_procparser(subparsers,
                                  'uncertainty',
                                  'calculate the phase uncertainty for all samples',
                                  uncertainty,
                                  'uncertainty')
    parser_unc.add_argument('-noise_bed_range',
                            type=float,
                            help='bed range under which the noise phasor will be calculated',
                            default=3000.)
    _add_def_args(parser_unc)

    # Phase Differencing
    parser_pdiff = _add_procparser(subparsers,
                                  'pdiff',
                                  'calculate correlation between the two acquisitions',
                                  phase_differencing,
                                  'pdiff')
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
                                  'unwrap the differenced phase profile from top to bottom',
                                  unwrap)
    _add_def_args(parser_unwrap)

    # Range Differencing
    parser_rdiff = _add_procparser(subparsers,
                                  'rdiff',
                                  'convert the differenced phase profile to range',
                                  range_differencing)
    _add_def_args(parser_rdiff)

    # Rotation, fill in data at all azimuths
    parser_rotate = _add_procparser(subparsers,
                                    'rotate',
                                    'use a rotational transform to find data at all azimuths',
                                    rotate,
                                    'rotated')
    parser_rotate.add_argument('-nthetas',
                               type=int,
                               help='number of theta values to rotate into',
                               default=100)
    parser_rotate.add_argument('-cross_pol_flip',
                             type=str,
                              help='flip the sign on one of the cross polarized terms.',
                              default=False)
    _add_def_args(parser_rotate)

    # Coherence 2D
    parser_coherence = _add_procparser(subparsers,
                                       'coherence',
                                       '2-dimensional coherence between HH and VV polarizations',
                                       coherence,
                                       'chhvv')
    parser_coherence.add_argument('-dtheta',
                                  type=float,
                                  help='size of coherence window in theta direction')
    parser_coherence.add_argument('-drange',
                                  type=float,
                                  help='size of coherence window in range direction')
    parser_coherence.set_defaults(dtheta=20.*np.pi/180., drange=100.)
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

    # Plot Apres Acquisition
    parser_plot = _add_procparser(subparsers,
                                  'plot',
                                  'plot apres data from a single acquisition',
                                  plot_apres,
                                  'plot')
    parser_plot.add_argument('-acq_type',
                                 type=str,
                                 help='Acquisition type',
                                 default=None,
                                 choices=['single', 'timediff', 'quadpol'])
    parser_plot.add_argument('-s',
                             action='store_true',
                             help='Save file (do not plt.show())')
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

    if args.name == 'load':
        apres_data, name = args.func(**vars(args))
    else:
        apres_data, empty_name = load(**vars(args))
        name = args.name
        args.func(apres_data, **vars(args))

    if args.name == 'plot':
        return
    elif args.o is not None:
        out_fn = args.o
        apres_data.save(out_fn)
    else:
        bn = os.path.splitext(args.fns[0])[0]
        if bn[-3:] == 'raw':
            bn = bn[:-6]
        out_fn = bn + '_{:s}.mat'.format(name)
        apres_data.save(out_fn)

# --------------------------------------------------------
### Loading ###
# --------------------------------------------------------

def load(fns='', acq_type=None,**kwargs):
    if acq_type == 'single':
        apres_data = load_apres.load_apres(fns)
        name = 'apraw'
    elif acq_type == 'timediff':
        if len(fns) == 1:
            apres_data = load_time_diff.load_time_diff(fns[0],load_single_acquisitions=False)
        else:
            apres_data = load_time_diff.load_time_diff(fns)
        name = 'tdraw'
    elif acq_type == 'quadpol':
        if len(fns) == 1:
            apres_data = load_quadpol.load_quadpol(fns[0],load_single_pol=False)
        else:
            apres_data = load_quadpol.load_quadpol(fns)
        name = 'qpraw'

    if acq_type is None:
        try:
            apres_data = load_apres.load_apres(fns)
            name = 'apraw'
        except:
            Warning('Cannot load as single Apres acquisition, trying alternative acq_type...')
        try:
            if len(fns) == 1:
                apres_data = load_time_diff.load_time_diff(fns[0],load_single_acquisitions=False)
            else:
                apres_data = load_time_diff.load_time_diff(fns)
            name = 'tdraw'
        except:
            Warning('Cannot load as timediff Apres acquisition, trying alternative acq_type...')
        try:
            if len(fns) == 1:
                apres_data = load_quadpol.load_quadpol(fns[0],load_single_pol=False)
            else:
                apres_data = load_quadpol.load_quadpol(fns)
            name = 'qpraw'
        except:
            Warning('Cannot load as quadpol Apres acquisition, trying alternative acq_type...')

    return apres_data, name

# --------------------------------------------------------
### Full Processing Flow ###
# --------------------------------------------------------

def single_processing(dat, p=2, max_range=4000., num_chirps=0., noise_bed_range=3000., **kwargs):
    """Full processing flow for ApresData object.
    Range conversion, stacking, uncertainty."""
    dat.apres_range(p,max_range)
    if num_chirps == 0.:
        dat.stacking()
    else:
        dat.stacking(num_chirps)
    dat.phase_uncertainty(noise_bed_range)


def time_diff_processing(diffdat, win=20, step=20, thresh=0.95,
                      strain_window=(200,1000), w_surf=-0.15, **kwargs):
    diffdat.phase_diff(win,step)
    diffdat.phase_unwrap(win,thresh)
    diffdat.range_diff()
    diffdat.strain_rate(strain_window=strain_window,w_surf=w_surf)
    diffdat.bed_pick()


def quadpol_processing(dat, nthetas=100, dtheta=20.0*np.pi/180.,
                    drange=100., Wn=0., fs=0., cross_pol_flip=False, **kwargs):
    """Full processing flow for ApresQuadPol object.
    Range conversion, stacking, uncertainty."""
    dat.rotational_transform(n_thetas=nthetas, cross_pol_flip=cross_pol_flip)
    dat.find_cpe()
    dat.coherence2d(delta_theta=dtheta, delta_range=drange)

# --------------------------------------------------------
### Individual Processing Functions ###
# --------------------------------------------------------

def range_conversion(dat, p=2, max_range=4000, **kwargs):
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


def phase_differencing(diffdat, win=20, step=20, **kwargs):
    diffdat.phase_diff(win,step)


def unwrap(diffdat,win=20,thresh=.95, **kwargs):
    diffdat.phase_unwrap(win,thresh)


def range_differencing(diffdat, **kwargs):
    diffdat.range_diff()


def rotate(dat, nthetas=100, cross_pol_flip=False, **kwargs):
    """Range conversion."""
    dat.rotational_transform(n_thetas=nthetas, cross_pol_flip=cross_pol_flip)


def coherence(dat, dtheta=20.0*np.pi/180., drange=100., **kwargs):
    """Stack chirps."""
    dat.coherence2d(delta_theta=dtheta, delta_range=drange)


def cross_polarized_extinction(dat,
                               Wn=0., fs=0., **kwargs):
    """Calculate uncertainty."""
    dat.find_cpe(Wn=Wn)

# --------------------------------------------------------
### Plotting ###
# --------------------------------------------------------

def plot_apres(dat, acq_type=None, s=False, o=None, o_fmt='png',
                 dpi=300, **kwargs):
    if type(dat.flags) is ApresFlags:
        plot.plot_apres(dat, s=s, o=o, ftype=o_fmt, dpi=dpi)
    elif type(dat.flags) is TimeDiffFlags:
        plot.plot_apres_diff(dat, s=s, o=o, ftype=o_fmt, dpi=dpi)
    elif type(dat.flags) is QuadPolFlags:
        plot.plot_apres_quadpol(dat, s=s, o=o, ftype=o_fmt, dpi=dpi)


if __name__ == '__main__':
    main()
