#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

Migration routines for ImpDAR

Much of this code is either directly referencing or written from older scripts in SeisUnix:
https://github.com/JohnWStockwellJr/SeisUnix/wiki
Options are:
    Kirchhoff (diffraction summation)
    Stolt (frequency wavenumber, constant velocity)
    Gazdag (phase shift, either constant or depth-varying velocity)
    SeisUnix (reference su routines directly)

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Mar 12 2019

"""

from __future__ import print_function
import os
import subprocess as sp
import numpy as np


def migrationSeisUnix(dat,
                      mtype='sumigtk',
                      vel=1.69e8,
                      vel_fn=None,
                      tmig=0,
                      verbose=1,
                      nxpad=100,
                      htaper=100,
                      vtaper=1000,
                      nz=None,
                      dz=None,
                      quiet=False
                      ):
    """

    Migration through Seis Unix. For now only three options:
    ---------
    1) sumigtk - Migration via T-K domain method for common-midpoint stacked data
    2) sugffd - Fourier finite difference migration for zero-offset data. This method is a hybrid
                migration which combines the advantages of phase shift and finite difference migrations.
    3) sustolt - Stolt migration for stacked data or common-offset gathers


    Parameters
    ---------
    dat: data as a class in the ImpDAR format
    vel: wave velocity, default is for ice
    vfile=         name of file containing velocities
    dx:  distance between successive
    fmax=Nyquist            maximum frequency
    tmig=0.0                times corresponding to interval velocities in vmig
    nxpad=0                 number of cdps to pad with zeros before FFT
    ltaper=0                length of linear taper for left and right edges
    verbose=0               =1 for diagnostic print
    dt=from header(dt) or  .004    time sampling interval
    ft=0.0                 first time sample
    ntau=nt(from data)     number of migrated time samples
    dtau=dt(from header)   migrated time sampling interval
    ftau=ft                first migrated time sample
    Q=1e6                  quality factor
    ceil=1e6               gain ceiling beyond which migration ceases

    Output
    ---------
    dat: data as a class in the ImpDAR format (with dat.data now being migrated data)

    """

    if sp.Popen(['which', mtype]).wait() != 0:
        raise FileNotFoundError('Cannot find chosen SeisUnix migration routine,' + mtype + '. Either install or choose a different migration routine.')

    # save to seisunix format for migration with SU routines
    out_fn = os.path.splitext(dat.fn)[0] + '.sgy'
    dat.save_as_segy(out_fn)

    # Get the trace spacing
    if np.mean(dat.trace_int) <= 0:
        Warning("The trace spacing, variable 'dat.trace_int', should be greater than 0. Using gradient(dat.dist) instead.")
        trace_int = np.gradient(dat.dist)
    else:
        trace_int = dat.trace_int
    dx = np.mean(trace_int)
    if nz is None:
        nz = dat.snum
    if dz is None:
        dz = 169 * dat.travel_time[-1] / 2 / dat.snum

    # Do the migration through the command line
    segy_name = os.path.splitext(dat.fn)[0]
    bin_fn = os.path.splitext(dat.fn)[0] + '_mig.bin'

    if quiet:
        stderr = sp.PIPE
    else:
        stderr = None

    ps1 = sp.Popen(['segyread', 'tape=' + segy_name + '.sgy'], stdout=sp.PIPE, stderr=stderr)
    ps2 = sp.Popen(['segyclean'], stdin=ps1.stdout, stdout=sp.PIPE, stderr=stderr)
    # Time Wavenumber
    if mtype == 'sumigtk':
        ps3 = sp.Popen(['sumigtk',
                        'tmig={:f}'.format(tmig),
                        'vmig={:f}'.format(vel * 1.e-6),
                        'verbose=' + str(verbose),
                        'nxpad={:d}'.format(int(nxpad)),
                        'ltaper={:d}'.format(htaper),
                        'dxcdp={:f}'.format(dx)],
                       stdout=sp.PIPE,
                       stderr=stderr,
                       stdin=ps2.stdout)

    # Fourier Finite Difference
    elif mtype == 'sumigffd':
        if vel_fn is None:
            raise ValueError('vel_fn needed for gffd')
        ps3 = sp.Popen(['sumigffd',
                        'vfile=' + vel_fn,
                        'nz={:d}'.format(nz),
                        'dz={:f}'.format(dz),
                        'dt={:f}'.format(dat.dt * 1.0e-6),
                        'dx={:f}'.format(dx)],
                       stdout=sp.PIPE,
                       stderr=stderr,
                       stdin=ps2.stdout)
    # Stolt
    elif mtype == 'sustolt':
        ps3 = sp.Popen(['sustolt',
                        'tmig={:f}'.format(tmig),
                        'vmig={:f}'.format(vel * 1.0e-6),
                        'verbose=' + str(verbose),
                        'lstaper={:d}'.format(htaper),
                        'lbtaper={:d}'.format(vtaper),
                        'dxcdp={:f}'.format(dx),
                        'cdpmin=0',
                        'cdpmax={:d}'.format(dat.tnum)],
                       stdout=sp.PIPE,
                       stderr=stderr,
                       stdin=ps2.stdout)
    # Stolt
    else:
        ps1.stdout.close()
        ps2.communicate()

        raise ValueError('The SeisUnix migration routine', mtype, 'has not been implemented in ImpDAR. Optionally, use ImpDAR to convert to SegY and run the migration in the command line.')


    ps4 = sp.Popen(['sustrip', segy_name + '_' + mtype + '.sgy'], stdin=ps3.stdout, stderr=stderr, stdout=sp.PIPE)
    with open(bin_fn, 'wb') as fout:
        fout.write(ps4.communicate()[0])
    with open(bin_fn, 'rb') as fid:
        data_flat = np.fromfile(fid, np.float32)
    for ps in [ps1, ps2, ps3, ps4]:
        ps.wait()
        ps.stdout.close()
        try:
            ps.stderr.close()
        except AttributeError:
            pass

    dat.data = np.transpose(np.reshape(data_flat, (dat.tnum, dat.snum)))

    # Cleanup but don't worry if we fail
    try:
        os.remove(bin_fn)
        os.remove('header')
        os.remove('binary')
        os.remove(segy_name + '.sgy')
    except FileNotFoundError:
        pass
    return dat


def _check_data_shape(dat):
    if np.size(dat.data, 1) != dat.tnum or np.size(dat.data, 0) != dat.snum:
        raise ValueError('The input array must be of size (tnum,snum)')
