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
from ._ApresDataProcessing import coherence
from scipy.signal import butter, filtfilt

try:
    USE_C = True
    from .coherence import coherence2d_loop
except ImportError:
    USE_C = False


def rotational_transform(self,theta_start=0,theta_end=np.pi,n_thetas=100):
    """
    Azimuthal (rotational) shift of principal axes
    at the transmitting and receiving antennas
    Mott, 2006

    Parameters
    --------
    S : array
        2-d array with [[shh,svh][shv,svv]] of complex numbers
    theta : complex
            rotational offset
    """

    self.thetas = np.linspace(theta_start,theta_end,n_thetas)

    self.HH = np.empty((len(self.range),len(self.thetas))).astype(np.complex)
    self.HV = np.empty((len(self.range),len(self.thetas))).astype(np.complex)
    self.VH = np.empty((len(self.range),len(self.thetas))).astype(np.complex)
    self.VV = np.empty((len(self.range),len(self.thetas))).astype(np.complex)

    for i,theta in enumerate(self.thetas):
        self.HH[:,i] = self.shh*np.cos(theta)**2.+(self.svh+self.shv)*np.sin(theta)*np.cos(theta)+self.svv*np.sin(theta)**2
        self.HV[:,i] = self.shv*np.cos(theta)**2.+(self.svv-self.shh)*np.sin(theta)*np.cos(theta)-self.svh*np.sin(theta)**2
        self.VH[:,i] = self.svh*np.cos(theta)**2.+(self.svv-self.shh)*np.sin(theta)*np.cos(theta)-self.shv*np.sin(theta)**2
        self.VV[:,i] = self.svv*np.cos(theta)**2.-(self.svh+self.shv)*np.sin(theta)*np.cos(theta)+self.shh*np.sin(theta)**2

    self.flags.rotation = np.array([1,n_thetas])


def coherence2d(self, delta_theta=20.0*np.pi/180., delta_range=100.):
    """
    Coherence between two 2-d images (e.g. c_hhvv).
    Jordan et al. (2019) eq. 19

    Parameters
    --------
    delta_theta : float
            window size in the azimuthal dimension; default: 20 degrees
    delta_range: float
            window size in the vertical; default: 100.
    """
    nrange = int(delta_range//abs(self.range[0]-self.range[1]))
    ntheta = int(delta_theta//abs(self.thetas[0]-self.thetas[1]))

    # THs, Rs = self.thetas, self.range
    # theta_start = THs[:ntheta]
    # theta_end = THs[-ntheta:]
    # THs = np.hstack((theta_end - np.pi, THs, theta_start + np.pi))

    HH_start = self.HH[:, :ntheta]
    HH_end = self.HH[:, -ntheta:]
    HH_ = np.hstack((HH_end, self.HH, HH_start))

    VV_start = self.VV[:, :ntheta]
    VV_end = self.VV[:, -ntheta:]
    VV_ = np.hstack((VV_end, self.VV, VV_start))

    chhvv = np.nan*np.ones_like(HH_).astype(np.complex)
    range_bins, azimuth_bins = HH_.shape[0], HH_.shape[1]

    t0 = time.time()
    if USE_C:
        chhvv = coherence2d_loop(np.ascontiguousarray(chhvv, dtype=np.cdouble),
                                 np.ascontiguousarray(HH_, dtype=np.cdouble),
                                 np.ascontiguousarray(VV_, dtype=np.cdouble),
                                 nrange,
                                 ntheta,
                                 range_bins,
                                 azimuth_bins)
    else:
        print('Beginning iteration through {:d} azimuths'.format(azimuth_bins - 2 * ntheta))
        print('Azimuth bin: ', end='')
        sys.stdout.flush()
        for i in range(HH_.shape[1]):
            if (i < ntheta) or (i > HH_.shape[1] - ntheta - 1):
                continue
            if (i - ntheta) % 10 == 0:
                print('{:d}'.format(i - ntheta), end='')
            else:
                print('.', end='')
            sys.stdout.flush()
            for j in range(HH_.shape[0]):
                imin, imax = i - ntheta, i + ntheta
                jmin, jmax = max(0, j - nrange), min(HH_.shape[0], j + nrange)
                chhvv[j,i] = coherence(VV_[j-nrange:j+nrange,i-ntheta:i+ntheta].flatten(),
                                        HH_[j-nrange:j+nrange,i-ntheta:i+ntheta].flatten())
        print('coherence calculation done')
    t1 = time.time()
    print('Execution with c code={:s} took {:6.2f}'.format(str(USE_C), t1 - t0))

    self.chhvv = chhvv[:, ntheta:-ntheta]

    self.flags.coherence = np.array([1, delta_theta, delta_range])

# --------------------------------------------------------------------------------------------

def phase_gradient2d(self,filt=None,Wn=0):
    """
    Depth-gradient of hhvv coherence image.
    Jordan et al. (2019) eq 23

    Parameters
    --------
    filt : string
            filter type; default lowpass
    Wn : float
            filter frequency
    """

    # Real and imaginary parts of the hhvv coherence
    R = np.real(self.chhvv).copy()
    I = np.imag(self.chhvv).copy()

    if filt is not None:
        if filt == 'lowpass':
            R = lowpass(R,Wn,1./self.dt)
            I = lowpass(I,Wn,1./self.dt)
        else:
            raise TypeError('Filter: %s has not been implemented yet.'%filt)

    dRdz = np.gradient(R,self.range,axis=0)
    dIdz = np.gradient(I,self.range,axis=0)

    self.dphi_dz = (R*dIdz-I*dRdz)/(R**2.+I**2.)

    self.flags.phasegradient = True

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
    Pa = np.transpose(np.transpose(P) - np.mean(P,axis=1))

    return Pa

# --------------------------------------------------------------------------------------------

def lowpass(data, Wn, fs, order=3):
    """
    lowpass filter

    Parameters
    --------
    data : array
            2-d array of azimuth-depth return
    Wn : float
            filter frequency
    fs : float
            sample frequency
    order : int
            order of filter
    """

    # Subset the array around nan values
    nan_idx = next(k for k, value in enumerate(data[:,0]) if ~np.isnan(value))
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

def birefringent_phase_shift(z,freq=300e6,eps_bi=0.00354,eps=3.15,c=3e8):
    """
    Two-way birefringent phase shift
    Jordan et al. (2019)

    Parameters
    ---------
    z: float
        depth
    freq: float
        center frequency
    eps_bi: float
        birefringent permittivity difference (i.e. eps_parallel - eps_perpendicular)
    eps: float
        mean permittivity (relative)
    c: float
        light speed in vacuum
    """

    #TODO: This function is not updated for ImpDAR yet
    delta = 4.*np.pi*freq/c*(z*eps_bi/(2.*np.sqrt(eps)))

    return delta

# --------------------------------------------------------------------------------------------

def phase_gradient_to_fabric(self,c=300e6,fc=300e6,delta_eps=0.035,eps=3.12):
    """
    """

    #TODO: This function is not updated for ImpDAR yet
    max_idx = np.argmax(self.dphi_dz,axis=1)
    E2E1 = (c/(4.*np.pi*fc))*(2.*np.sqrt(eps)/delta_eps)*self.dphi_dz[max_idx,np.arange(len(self.range))]

    return E2E1
