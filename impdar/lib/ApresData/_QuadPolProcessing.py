#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright Â© 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Process ApRES data

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 23 2019

"""

import numpy as np
from ._ApresDataProcessing import coherence

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

# --------------------------------------------------------------------------------------------

def copolarized_coherence(self,delta_theta=20*np.pi/180.,delta_range=100):

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

# --------------------------------------------------------------------------------------------

def copolarized_phase_gradient(self,delta_theta=20*np.pi/180.,delta_range=100):

    THs,Rs = np.meshgrid(self.thetas,self.range)
    nrange = int(delta_range//abs(self.range[0]-self.range[1]))

    R = np.real(self.chhvv)
    I = np.imag(self.chhvv)

    dRdz = np.nan*np.ones_like(self.chhvv).astype(np.float)
    dIdz = np.nan*np.ones_like(self.chhvv).astype(np.float)

    for i,theta in enumerate(self.thetas):
        print('i:',i,'theta:',np.round(theta,2))
        for j in range(len(self.range)):
            if j < nrange or j > len(self.range)-nrange:
                continue
            else:
                Rs_ij = R[j-nrange:j+nrange,i]
                Is_ij = I[j-nrange:j+nrange,i]
                ds_ij = Rs[j-nrange:j+nrange,i]

                Ridxs = ~np.isnan(Rs_ij) & ~np.isnan(ds_ij)
                Iidxs = ~np.isnan(Is_ij) & ~np.isnan(ds_ij)

                if len(Rs_ij[Ridxs])<3 or len(Is_ij[Iidxs])<3:
                    dRdz[j,i] = np.nan
                    dIdz[j,i] = np.nan
                else:
                    dRdz[j,i] = np.polyfit(ds_ij[Ridxs],Rs_ij[Ridxs],1)[0]
                    dIdz[j,i] = np.polyfit(ds_ij[Iidxs],Is_ij[Iidxs],1)[0]

    self.dphi_dz = (R*dIdz-I*dRdz)/(R**2.+I**2.)

# --------------------------------------------------------------------------------------------

def birefringent_phase_shift(z,freq=200e6,eps_bi=0.00354,eps=3.15,c=3e8):
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
    delta = 4.*np.pi*freq/c*(z*eps_bi/(2.*np.sqrt(eps)))
    return delta

# --------------------------------------------------------------------------------------------

def phase_gradient_to_fabric(self,c = 300e6,fc = 200e6,Δϵ = 0.035,ϵ = 3.12):
    max_idx = np.argmax(self.dϕdz,axis=1)
    min_idx = np.argmin(self.dϕdz,axis=1)

    E2E1 = (c/(4.*np.pi*fc))*(2.*np.sqrt(ϵ)/Δϵ)*self.dϕdz[max_idx,np.arange(len(d))]
    return E2E1
