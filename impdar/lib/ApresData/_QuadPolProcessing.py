#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright Â© 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Process quad-polarized ApRES data

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Oct 13 2020

"""
import sys
import time

import numpy as np
from ._TimeDiffProcessing import coherence
from scipy.signal import butter, filtfilt

from ..ImpdarError import ImpdarError

try:
    USE_C = True
    from .coherence import coherence2d_loop
except ImportError:
    USE_C = False


def rotational_transform(self, theta_start=0, theta_end=np.pi, n_thetas=100,
                         cross_pol_exception=False,
                         cross_pol_flip=False,
                         flip_force=False):
    """
    Azimuthal (rotational) shift of principal axes
    at the transmitting and receiving antennas
    Mott, 2006

    Parameters
    ---------
    self: class
        ApresQuadPol object
    theta_start: float
        Starting point for array of azimuths
    theta_end: float
        Ending point for array of azimuths
    n_thetas: int
        number of thetas to rotate through
    cross_pol_exception: bool
        continue even if the cross-polarized terms do not look correct
    cross_pol_flip: bool or string
        Which cross-polarized term to flip if they look inverted
    """

    # Antennas are not symmetric on 180 degree
    # Issues have come up with flipped antennas, so check
    if abs(np.sum(np.imag(self.shv)+np.imag(self.svh))) < \
            abs(np.sum(np.imag(self.shv)-np.imag(self.svh))) or \
            abs(np.sum(np.real(self.shv)+np.real(self.svh))) < \
            abs(np.sum(np.real(self.shv)-np.real(self.svh))) or \
            flip_force:
        if cross_pol_exception:
            pass
        elif cross_pol_flip == 'HV':
            Warning('Flipping sign of cross-polarized term HV')
            self.shv *= -1.
        elif cross_pol_flip == 'VH':
            Warning('Flipping sign of cross-polarized term VH')
            self.svh *= -1.
        else:
            raise ValueError('Cross-polarized terms are of the opposite sign, check and update.')

    self.thetas = np.linspace(theta_start, theta_end, n_thetas)

    self.HH = np.empty((len(self.range), len(self.thetas))).astype(np.cdouble)
    self.HV = np.empty((len(self.range), len(self.thetas))).astype(np.cdouble)
    self.VH = np.empty((len(self.range), len(self.thetas))).astype(np.cdouble)
    self.VV = np.empty((len(self.range), len(self.thetas))).astype(np.cdouble)

    for i, theta in enumerate(self.thetas):
        self.HH[:, i] = self.shh*np.cos(theta)**2. + \
            (self.svh + self.shv)*np.sin(theta)*np.cos(theta) + \
            self.svv*np.sin(theta)**2
        self.HV[:, i] = self.shv*np.cos(theta)**2. + \
            (self.svv - self.shh)*np.sin(theta)*np.cos(theta) - \
            self.svh*np.sin(theta)**2
        self.VH[:, i] = self.svh*np.cos(theta)**2. + \
            (self.svv - self.shh)*np.sin(theta)*np.cos(theta) - \
            self.shv*np.sin(theta)**2
        self.VV[:, i] = self.svv*np.cos(theta)**2. - \
            (self.svh + self.shv)*np.sin(theta)*np.cos(theta) + \
            self.shh*np.sin(theta)**2

    self.flags.rotation = np.array([1, n_thetas])


def coherence2d(self, delta_theta=20.0*np.pi/180., delta_range=100.,
                force_python=False):
    """
    Coherence between two 2-d images (e.g. c_hhvv).
    Jordan et al. (2019) eq. 19

    Parameters
    ---------
    self: class
        ApresQuadPol object
    delta_theta : float
            window size in the azimuthal dimension; default: 20 degrees
    delta_range: float
            window size in the vertical; default: 100.
    """

    if self.flags.rotation[0] != 1:
        raise ImpdarError('Rotate the quad-pol acquisition before \
                          calling this function.')

    nrange = int(delta_range//abs(self.range[0]-self.range[1]))
    ntheta = int(delta_theta//abs(self.thetas[0]-self.thetas[1]))

    HH_start = self.HH[:, :ntheta]
    HH_end = self.HH[:, -ntheta:]
    HH_ = np.hstack((HH_end, self.HH, HH_start))

    VV_start = self.VV[:, :ntheta]
    VV_end = self.VV[:, -ntheta:]
    VV_ = np.hstack((VV_end, self.VV, VV_start))

    chhvv = np.nan*np.ones_like(HH_).astype(np.cdouble)
    range_bins, azimuth_bins = HH_.shape[0], HH_.shape[1]

    t0 = time.time()
    if USE_C and not force_python:
        chhvv = coherence2d_loop(np.ascontiguousarray(chhvv, dtype=np.cdouble),
                                 np.ascontiguousarray(HH_, dtype=np.cdouble),
                                 np.ascontiguousarray(VV_, dtype=np.cdouble),
                                 nrange,
                                 ntheta,
                                 range_bins,
                                 azimuth_bins)
        t1 = time.time()
        print('Execution with c code={:s} took {:6.2f}'.format(str(USE_C), t1 - t0))
    else:
        print('Beginning iteration through {:d} azimuths'.format(azimuth_bins - 2 * ntheta))
        print('Azimuth bin: ', end='')
        sys.stdout.flush()
        for i in range(azimuth_bins):
            if (i < ntheta) or (i > azimuth_bins - ntheta - 1):
                continue
            if (i - ntheta) % 10 == 0:
                print('{:d}'.format(i - ntheta), end='')
            else:
                print('.', end='')
            sys.stdout.flush()
            for j in range(range_bins):
                imin, imax = i - ntheta, i + ntheta
                jmin, jmax = max(0, j - nrange), min(range_bins-1, j + nrange)
                chhvv[j, i] = coherence(HH_[jmin:jmax,imin:imax].flatten(),
                                        VV_[jmin:jmax,imin:imax].flatten())
        print('coherence calculation done')
        t1 = time.time()
        print('Execution with c code={:s} took {:6.2f}'.format(str(False), t1 - t0))

    self.chhvv = chhvv[:, ntheta:-ntheta]

    # If the cpe axis has already been identified, get coherence along it
    if self.flags.cpe is True:
        self.chhvv_cpe = self.chhvv[np.arange(self.snum), self.cpe_idxs]

    self.flags.coherence = np.array([1, delta_theta, delta_range])

# --------------------------------------------------------------------------------------------

def phase_gradient2d(self, filt=None, Wn=0):
    """
    Depth-gradient of hhvv coherence image.
    Jordan et al. (2019) eq 23

    Parameters
    ---------
    self: class
        ApresQuadPol object
    filt : string
            filter type; default lowpass
    Wn : float
            filter frequency
    """

    if self.flags.coherence[0] != 1:
        raise ImpdarError('Calculate coherence before calling this function.')

    # Real and imaginary parts of the hhvv coherence
    R_ = np.real(self.chhvv).copy()
    I_ = np.imag(self.chhvv).copy()

    # Filter image before gradient calculation
    if filt is not None:
        if filt == 'lowpass':
            R_ = lowpass(R_, Wn, 1./self.dt)
            I_ = lowpass(I_, Wn, 1./self.dt)
        else:
            raise TypeError('Filter: %s has \
                            not been implemented yet.' % filt)

    # Depth gradient for each component
    dRdz = np.gradient(R_, self.range, axis=0)
    dIdz = np.gradient(I_, self.range, axis=0)

    # Phase-depth gradient from Jordan et al. (2019) eq. 23
    self.dphi_dz = (R_*dIdz-I_*dRdz)/(R_**2.+I_**2.)

    # If the cpe axis has already been identified, get phase_gradient along it
    if self.flags.cpe is True:
        self.dphi_dz_cpe = self.dphi_dz[np.arange(self.snum), self.cpe_idxs]

    self.flags.phasegradient = True


def find_cpe(self, Wn=50, rad_start=np.pi/4., rad_end=3.*np.pi/4.,
             *args, **kwargs):
    """
    Find the cross-polarized extinction axis.
    Ershadi et al. (2022)

    Parameters
    ---------
    self: class
        ApresQuadPol object
    Wn: float
        Filter frequency
    rad_start: float
        Starting point for the azimuthal window to look within for cpe axis
    rad_end: float
        Ending point for the azimuthal window to look within for cpe axis
    """

    if self.flags.rotation[0] != 1:
        raise ImpdarError('Rotate the quad-pol acquisition before \
                          calling this function.')

    # Power anomaly
    HV_pa = power_anomaly(self.HV.copy())
    # Lowpass filter
    HV_pa = lowpass(HV_pa, Wn, 1./self.dt)

    # Find the index of the cross-polarized extinction axis
    idx_start = np.argmin(abs(self.thetas-rad_start))
    idx_stop = np.argmin(abs(self.thetas-rad_end))
    CPE_idxs = np.empty_like(self.range).astype(int)
    for i in range(len(CPE_idxs)):
        CPE_idxs[i] = np.argmin(HV_pa[i, idx_start:idx_stop])
    CPE_idxs += idx_start

    self.cpe_idxs = CPE_idxs
    self.cpe = np.array([self.thetas[i] for i in CPE_idxs]).astype(float)

    # If the coherence calculation has already been done
    # get it along the cpe axis
    if self.flags.coherence[0] == 1.:
        self.chhvv_cpe = self.chhvv[np.arange(self.snum), self.cpe_idxs]
    # If the phase gradient calculation has already been done
    # get it along the cpe axis
    if self.flags.phasegradient:
        self.dphi_dz_cpe = self.dphi_dz[np.arange(self.snum), self.cpe_idxs]

    self.flags.cpe = True

# --------------------------------------------------------------------------------------------

def phase_gradient_to_fabric(self, c=300e6, fc=300e6, delta_eps=0.035, eps=3.12):
    """
    Calculate the fabric strength from phase gradient
    Jordan et al. (2019)

    Parameters
    ---------
    self: class
        ApresQuadPol object
    c: float
        speed of light
    fc: float
        center frequency
    delta_eps: float
        ice crystal permittivity anisotropy
    eps: float
        permittivity of ice
    """

    if not hasattr(self, 'dphi_dz_cpe'):
        raise AttributeError("Get the phase gradient along CPE axis \
                             before calling this function.")

    self.e2e1 = (c/(4.*np.pi*fc))*(2.*np.sqrt(eps)/delta_eps)*self.dphi_dz_cpe

# --------------------------------------------------------------------------------------------

def power_anomaly(data):
    """
    Calculate the power anomaly from mean for a given image.
    Ershadi et al. (2021) eq 21

    Parameters
    --------
    data : array
            2-d array of azimuth-depth return
    """

    # Calculate Power
    P = 10.*np.log10((data)**2.)
    # Remove the mean for each row
    Pa = np.transpose(np.transpose(P) - np.nanmean(P, axis=1))

    return Pa

# --------------------------------------------------------------------------------------------

def lowpass(data, Wn, fs, order=3):
    """
    lowpass filter

    Parameters
    --------
    data: array
            2-d array of azimuth-depth return
    Wn: float
            filter frequency
    fs: float
            sample frequency
    order: int
            order of filter
    """

    # Subset the array around nan values
    nan_idx = next(k for k, value in enumerate(data[:, 1]) if ~np.isnan(value))
    if nan_idx != 0:
        data_sub = data[nan_idx:-nan_idx+1]
    else:
        data_sub = data.copy()

    # Get the filter coefficients
    b, a = butter(order, Wn, btype='low', fs=fs)
    # Filter (need to transpose to filter along depth axis)
    data_filtered = filtfilt(b, a, data_sub, axis=0)
    # Insert into original array
    if nan_idx != 0:
        data[nan_idx:-nan_idx+1] = data_filtered
        return data
    else:
        return data_filtered

# --------------------------------------------------------------------------------------------

def azimuthal_rotation(data,thetas,azi):
    """
    Rotate a quad-pol image based on some known antenna orientation.

    Parameters
    --------
    data : array
            2-d array of azimuth-depth return
    thetas : array
            series of azimuths for the image
    azi : array
            azimuth of the antenna orientation on measurement

    Output
    --------
    data : array
            rotated image
    """

    thetas += azi
    if azi < 0:
        idx_clip = np.argwhere(thetas > 0)[0][0]
        hold = data[:, idx_clip:]
        data = np.append(hold, data[:, :idx_clip], axis=1)
    elif azi > 0:
        idx_clip = np.argwhere(thetas > np.pi)[0][0]
        hold = data[:, idx_clip:]
        data = np.append(hold, data[:, :idx_clip], axis=1)
    thetas -= azi

    return data
