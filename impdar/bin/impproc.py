#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make an executable for single actions of impulse radar processing.

All functionality probably overlaps with impdar, but the call is much cleaner.
You get a lot more flexibility on things like keyword arguments.
However, you are limited to one processing step.
"""
import sys
import os.path
import argparse

from impdar.lib.load import load, FILETYPE_OPTIONS
from impdar.lib.process import concat
from impdar.lib.gpslib import interp as interpdeep


def _get_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Choose a processing step')

    # Horizontal window filter
    parser_hfilt = _add_procparser(subparsers,
                                   'hfilt',
                                   'Horizontally filter the data by \
                                       subtracting the average trace \
                                       from a window',
                                   hfilt,
                                   defname='hfilted')
    parser_hfilt.add_argument('start_trace',
                              type=int,
                              help='First trace of representative subset')
    parser_hfilt.add_argument('end_trace',
                              type=int,
                              help='Last trace of representative subset')
    _add_def_args(parser_hfilt)

    # Adaptive horizontal filter
    parser_ahfilt = _add_procparser(subparsers,
                           'ahfilt',
                           'Horizontally filter the data adaptively',
                           ahfilt,
                           defname='ahfilt')
    parser_ahfilt.add_argument('win',
                                type=int,
                                help='Number of traces to include in the moving average')
    _add_def_args(parser_ahfilt)

    # Simply reverse the files
    _add_simple_procparser(subparsers,
                           'rev',
                           'Reverse the data',
                           rev,
                           defname='rev')

    # Concatenate
    _add_simple_procparser(subparsers,
                           'cat',
                           'Concatenate the data',
                           concat,
                           defname='cat')

    # Elevation correction
    _add_simple_procparser(subparsers,
                           'elev',
                           'Elevation correct',
                           elev,
                           defname='elev')

    # Restack
    parser_restack = _add_procparser(subparsers,
                                     'restack',
                                     'Restack to interval',
                                     restack,
                                     defname='restacked')
    parser_restack.add_argument('traces',
                                type=int,
                                help='Number of traces to stack. \
                                        Must be an odd number')
    _add_def_args(parser_restack)

    # Range gain
    parser_rgain = _add_procparser(subparsers,
                                   'rgain',
                                   'Add a range gain',
                                   rgain,
                                   defname='rgain')
    parser_rgain.add_argument('-slope',
                              type=float,
                              default=0.1,
                              help='Slope of linear range gain. Default 0.1')
    _add_def_args(parser_rgain)

    # Automatic range gain
    parser_agc = _add_procparser(subparsers,
                                 'agc',
                                 'Add an automatic gain',
                                 agc,
                                 defname='agc')
    parser_agc.add_argument('-window',
                            type=int,
                            default=50,
                            help='Number of samples to average')
    _add_def_args(parser_agc)

    # Vertical bandpass
    parser_vbp = _add_procparser(subparsers,
                                 'vbp',
                                 'Vertically bandpass the data',
                                 vbp,
                                 defname='bandpassed')
    parser_vbp.add_argument('low_MHz',
                            type=float,
                            help='Lowest frequency passed (in MHz)')
    parser_vbp.add_argument('high_MHz',
                            type=float,
                            help='Highest frequency passed (in MHz)')
    _add_def_args(parser_vbp)

    # Horizontal bandpass
    parser_hbp = _add_procparser(subparsers,
                                 'hbp',
                                 'Horizontally bandpass the data',
                                 hbp,
                                 defname='hbp')
    parser_hbp.add_argument('low',
                            type=float,
                            help='Lowest frequency passed (in wavelength)')
    parser_hbp.add_argument('high',
                            type=float,
                            help='Highest frequency passed (in wavelength)')
    _add_def_args(parser_hbp)

    parser_lp = _add_procparser(subparsers,
                                'lp',
                                'Horizontally lowpass the data',
                                lp,
                                defname='lp')
    parser_lp.add_argument('low',
                           type=float,
                           help='Lowest frequency passed (in wavelength)')
    _add_def_args(parser_lp)

    # Crop in the vertical
    parser_crop = _add_procparser(subparsers,
                                  'crop',
                                  'Crop the data in the vertical',
                                  crop,
                                  defname='cropped')
    parser_crop.add_argument('top_or_bottom',
                             choices=['top', 'bottom'],
                             help='Remove from the top or bottom')
    parser_crop.add_argument('dimension',
                             choices=['snum', 'twtt', 'depth', 'pretrig'],
                             help='Set the bound in terms of snum \
                                    (sample number), twtt (two way travel \
                                    time in microseconds), depth (m, \
                                    calculated using the nmo_depth or a light \
                                    speed of 1.69e8m/s if it doesn\'t, or \
                                    pretrig (the recorded trigger sample)')
    parser_crop.add_argument('lim',
                             type=float,
                             help='The cutoff value')
    _add_def_args(parser_crop)

    # Crop in the horizontal
    parser_hcrop = _add_procparser(subparsers,
                                   'hcrop',
                                   'Crop the data in the horizontal',
                                   hcrop,
                                   defname='hcropped')
    parser_hcrop.add_argument('left_or_right',
                              choices=['left', 'right'],
                              help='Remove from the left or right')
    parser_hcrop.add_argument('dimension',
                              choices=['tnum', 'dist'],
                              help='Set the bound in terms of tnum \
                                      (trace number, 1 indexed) or dist \
                                      (distance in km)')
    parser_hcrop.add_argument('lim',
                              type=float,
                              help='The cutoff value')
    _add_def_args(parser_hcrop)

    # Normal move-out
    parser_nmo = _add_procparser(subparsers,
                                 'nmo',
                                 'Normal move-out correction',
                                 nmo,
                                 defname='nmo')
    parser_nmo.add_argument('ant_sep',
                            type=float,
                            help='Antenna separation')
    parser_nmo.add_argument('--uice',
                            type=float,
                            default=1.69e8,
                            help='Speed of light in ice in m/s \
                                    (default 1.69e8)')
    parser_nmo.add_argument('--uair',
                            type=float,
                            default=3.0e8,
                            help='Speed of light in air in m/s \
                                   (default 3.0e8)')
    parser_nmo.add_argument('--const_firn_offset',
                            type=float,
                            default=None,
                            help='A constant value added to depth to account for firn. Default None (0.0).')
    parser_nmo.add_argument('--rho_profile',
                            type=str,
                            default=None,
                            help='Filename for a depth density profile to \
                                    correct wave velocity.')
    _add_def_args(parser_nmo)

    # Reinterpolate GPS
    parser_interp = _add_procparser(subparsers,
                                    'interp',
                                    'Reinterpolate GPS',
                                    interp,
                                    defname='interp')
    parser_interp.add_argument('spacing',
                               type=float,
                               help='New spacing of radar traces, \
                                       in meters')
    parser_interp.add_argument('--gps_fn',
                               type=str,
                               help='CSV or mat file containing the GPS \
                                       information. .csv and .txt files are \
                                       assumed to be csv, .mat are mat. \
                                       Default is None--use associated \
                                       (presumably non-precision) GPS',
                               default=None)
    parser_interp.add_argument('--offset',
                               type=float,
                               default=0.0,
                               help='Offset from GPS time to radar time')
    parser_interp.add_argument('--minmove',
                               type=float,
                               default=1.0e-2,
                               help='Minimum movement to not be stationary')
    parser_interp.add_argument('--extrapolate',
                               action='store_true',
                               help='Extrapolate GPS data beyond bounds')
    _add_def_args(parser_interp)

    # GPS
    parser_geolocate = _add_procparser(subparsers,
                                       'geolocate',
                                       'GPS control',
                                       geolocate,
                                       defname='geolocate')
    parser_geolocate.add_argument('gps_fn',
                                  type=str,
                                  help='CSV or mat file containing the GPS \
                                         information. .csv and .txt files are \
                                         assumed to be csv, .mat are mat. \
                                         Default is None--use associated \
                                         (presumably non-precision) GPS')
    parser_geolocate.add_argument('--extrapolate',
                                  action='store_true',
                                  help='Extrapolate GPS data beyond bounds')
    parser_geolocate.add_argument('--guess',
                                  action='store_true',
                                  help='Guess at offset')
    _add_def_args(parser_geolocate)

    # Denoise
    parser_denoise = _add_procparser(subparsers,
                                     'denoise',
                                     'Denoising filter for the data image',
                                     denoise,
                                     defname='denoise')
    parser_denoise.add_argument('vert_win',
                                type=int,
                                help='Size of filtering window in vertical \
                                    (number of samples)')
    parser_denoise.add_argument('hor_win',
                                type=int,
                                help='Size of filtering window in horizontal \
                                    (number of traces)')
    _add_def_args(parser_denoise)

    # Migration
    parser_mig = _add_procparser(subparsers,
                                 'migrate',
                                 'Migration',
                                 mig,
                                 defname='migrated')
    parser_mig.add_argument('--mtype',
                            type=str,
                            default='phsh',
                            choices=['stolt', 'kirch', 'phsh', 'tk',
                                     'sumigtk', 'sustolt', 'sumigffd'],
                            help='Migration routines.')
    parser_mig.add_argument('--vel',
                            type=float,
                            default=1.69e8,
                            help='Speed of light in dielectric medium m/s  \
                                    default is for ice, 1.69e8)')
    parser_mig.add_argument('--vel_fn',
                            type=str,
                            default=None,
                            help='Filename for inupt velocity array. \
                                Column 1: velocities, \
                                Column 2: z locations, \
                                Column 3: x locations (optional)')
    parser_mig.add_argument('--nearfield',
                            action='store_true',
                            help='Boolean for nearfield operator in \
                                    Kirchhoff migration.')
    parser_mig.add_argument('--htaper',
                            type=int,
                            default=100,
                            help='Number of samples for horizontal taper')
    parser_mig.add_argument('--vtaper',
                            type=int,
                            default=1000,
                            help='Number of samples for vertical taper')
    parser_mig.add_argument('--nxpad',
                            type=int,
                            default=100,
                            help='Number of traces to pad with zeros for FFT')
    parser_mig.add_argument('--tmig',
                            type=int,
                            default=0,
                            help='Times for velocity profile')
    parser_mig.add_argument('--verbose',
                            type=int,
                            default=1,
                            help='Print output from SeisUnix migration')
    _add_def_args(parser_mig)

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
    parser.add_argument('--ftype',
                        type=str,
                        default='mat',
                        help='Type of file to load (default ImpDAR mat)',
                        choices=FILETYPE_OPTIONS)


def main():
    """Get arguments, process data, handle saving."""
    parser = _get_args()
    args = parser.parse_args(sys.argv[1:])
    if not hasattr(args, 'func'):
        parser.parse_args(['-h'])

    radar_data = load(args.ftype, args.fns)

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
        if ((len(radar_data) > 1) or (args.o[-1] == '/')):
            for d, f in zip(radar_data, args.fns):
                bn = os.path.split(os.path.splitext(f)[0])[1]
                if bn[-4:] == '_raw':
                    bn = bn[:-4]
                out_fn = os.path.join(args.o,
                                      bn + '_{:s}.mat'.format(args.name))
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
    """Perform some horizontal filtering."""
    dat.hfilt(ftype='hfilt', bounds=(start_trace, end_trace))


def ahfilt(dat, window_size=1000, **kwargs):
    """Adaptive horizontal filter."""
    dat.hfilt(ftype='adaptive', window_size=window_size)


def rev(dat, **kwargs):
    """Flip the data horizontally."""
    dat.reverse()


def elev(dat, **kwargs):
    """Move the data to start at the surface elevation (DO LAST)."""
    dat.elev_correct()


def vbp(dat, low_MHz=1, high_MHz=10000, **kwargs):
    """Vertically bandpass the data."""
    dat.vertical_band_pass(low_MHz, high_MHz)


def hbp(dat, low=1, high=10, **kwargs):
    """Horizontally band pass the data."""
    dat.horizontal_band_pass(low, high)


def lp(dat, low=1, **kwargs):
    """Low pass filter the data."""
    dat.lowpass(low)


def crop(dat, lim=0, top_or_bottom='top', dimension='snum', **kwargs):
    """Crop in the vertical."""
    dat.crop(lim, top_or_bottom=top_or_bottom, dimension=dimension)


def hcrop(dat, lim=0, left_or_right='left', dimension='tnum', **kwargs):
    """Crop in the horizontal."""
    dat.hcrop(lim, left_or_right=left_or_right, dimension=dimension)


def nmo(dat, ant_sep=0.0, uice=1.69e8, uair=3.0e8, rho_profile=None, **kwargs):
    """Move-out correction to account for antenna spacing."""
    dat.nmo(ant_sep, uice=uice, uair=uair, rho_profile=rho_profile)


def restack(dat, traces=1, **kwargs):
    """Restack to reduce size/noise."""
    dat.restack(traces)


def rgain(dat, slope=0.1, **kwargs):
    """Set range gain."""
    dat.rangegain(slope)


def agc(dat, window=50, scale_factor=50, **kwargs):
    """Automatically control gain."""
    dat.agc(window=window, scaling_factor=scale_factor)


def interp(dats, spacing, gps_fn, offset=0.0, minmove=1.0e-2,
           extrapolate=False, **kwargs):
    """Move data to constant spacing."""
    interpdeep(dats,
               spacing,
               fn=gps_fn,
               offset=offset,
               min_movement=minmove,
               extrapolate=extrapolate)


def geolocate(dats, gps_fn, extrapolate=False, guess=False, **kwargs):
    """Attach precision gps, with option to constant space."""
    interpdeep(dats,
               spacing=None,
               fn=gps_fn,
               extrapolate=extrapolate,
               guess_offset=guess)


def denoise(dat, vert_win=1, hor_win=10, noise=None, filter_type='wiener', **kwargs):
    """Despeckle."""
    dat.denoise(vert_win=vert_win, hor_win=hor_win, noise=noise, ftype=filter_type)


def mig(dat, mtype='stolt', vel=1.69e8, vtaper=100, htaper=100, tmig=0,
        verbose=0, vel_fn=None, nxpad=1, nearfield=False, **kwargs):
    """Migrate data."""
    dat.migrate(mtype,
                vel=vel,
                vtaper=vtaper,
                htaper=htaper,
                tmig=tmig,
                verbose=verbose,
                vel_fn=vel_fn,
                nxpad=nxpad,
                nearfield=nearfield)


if __name__ == '__main__':
    main()
