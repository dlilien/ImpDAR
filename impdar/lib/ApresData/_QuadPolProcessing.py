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

import numpy as np
from ._ApresDataProcessing import coherence
from scipy.signal import butter, filtfilt

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

# --------------------------------------------------------------------------------------------

def coherence2d(self,delta_theta=20*np.pi/180.,delta_range=100.):
    """

    Parameters
    --------
    delta_theta : float
            window size in the azimuthal dimension; default: 20 degrees
    delta_range: float
            window size in the vertical; default: 100.
    """

    # Create a me
    THs,Rs = np.meshgrid(self.thetas,self.range)

    nrange = int(delta_range//abs(self.range[0]-self.range[1]))
    ntheta = int(delta_theta//abs(self.thetas[0]-self.thetas[1]))

    theta_start = THs[:,1:ntheta+1]
    theta_end = THs[:,-ntheta-1:-1]
    THs = np.append(THs,theta_start+np.pi,axis=1)
    THs = np.append(theta_end-np.pi,THs,axis=1)

    R_start = Rs[:,1:ntheta+1]
    R_end = Rs[:,-ntheta-1:-1]
    Rs = np.append(Rs,R_start,axis=1)
    Rs = np.append(R_end,Rs,axis=1)

    HH_start = self.HH[:,1:ntheta+1]
    HH_end = self.HH[:,-ntheta-1:-1]
    HH_ = np.append(self.HH,HH_start,axis=1)
    HH_ = np.append(HH_end,HH_,axis=1)

    VV_start = self.VV[:,1:ntheta+1]
    VV_end = self.VV[:,-ntheta-1:-1]
    VV_ = np.append(self.VV,VV_start,axis=1)
    VV_ = np.append(VV_end,VV_,axis=1)

    chhvv = np.nan*np.ones_like(HH_).astype(np.complex)
    for i,θ in enumerate(THs[0]):
        for j in range(len(Rs[:,0])):
            if j < nrange or j > len(Rs[:,0])-nrange or i < ntheta or i > len(THs[0])-ntheta-1:
                continue
            else:
                chhvv[j,i] = coherence(HH_[j-nrange:j+nrange,i-ntheta:i+ntheta].flatten(),
                                         VV_[j-nrange:j+nrange,i-ntheta:i+ntheta].flatten())

    self.chhvv = chhvv[:,ntheta:-ntheta]

    self.flags.coherence = np.array([1,delta_theta,delta_range])

# --------------------------------------------------------------------------------------------

def phase_gradient2d(self,filt=None,Wn=0):

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

def power_anomaly(self):
    """
    """

    # Calculate Power
    P = 10.*np.log10((self.HH)**2.)
    # Remove the mean for each row
    Pa = np.transpose(np.transpose(P) - np.mean(P,axis=1))

    return Pa

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

# --------------------------------------------------------------------------------------------

def lowpass(data, Wn, fs, order=3):
    """
    """

    # Subset the array around nan values
    nan_idx = next(k for k, value in enumerate(data[:,0]) if ~np.isnan(value))
    data_sub = data[nan_idx:-nan_idx+1]

    # Get the filter coefficients
    b, a = butter(order, Wn, btype='low', fs=fs)
    # Filter (need to transpose to filter along depth axis)
    data_filtered = filtfilt(b, a, data_sub, axis=0)
    # Insert into original array
    data[nan_idx:-nan_idx+1] = data_filtered

    return data
