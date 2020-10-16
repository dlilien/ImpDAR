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
from _ApresDataProcessing import coherence

def rotational_transform(S,theta):
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

    shh = S[0,0]
    svh = S[0,1]
    shv = S[1,0]
    svv = S[1,1]

    S_ = np.empty_like(S)
    S_[0,0] = shh*np.cos(theta)**2.+(svh+shv)*np.sin(theta)*np.cos(theta)+svv*np.sin(theta)**2
    S_[0,1] = shv*np.cos(theta)**2.+(svv-shh)*np.sin(theta)*np.cos(theta)-svh*np.sin(theta)**2
    S_[1,0] = svh*np.cos(theta)**2.+(svv-shh)*np.sin(theta)*np.cos(theta)-shv*np.sin(theta)**2
    S_[1,1] = svv*np.cos(theta)**2.-(svh+shv)*np.sin(theta)*np.cos(theta)+shh*np.sin(theta)**2

    return S_

# --------------------------------------------------------------------------------------------

def copolarized_coherence(self,theta_start=0,theta_end=np.pi,n_thetas=100,
                            delta_theta = 20*np.pi/180.,delta_d = 100):

    thetas = np.linspace(theta_start,theta_end,n_thetas)
    Thetas,Ds = np.meshgrid(thetas,d)

    HH = np.empty((len(d),len(thetas))).astype(np.complex)
    VH = np.empty((len(d),len(thetas))).astype(np.complex)
    HV = np.empty((len(d),len(thetas))).astype(np.complex)
    VV = np.empty((len(d),len(thetas))).astype(np.complex)

    for i,θ in enumerate(thetas):
        S = np.array([[self.hh,self.vh],[self.hv,self.vv]])
        S_ = rotational_transform(S,θ)
        HH[:,i] = S_[0,0]
        HV[:,i] = S_[0,1]
        VH[:,i] = S_[1,0]
        VV[:,i] = S_[1,1]

    nd = int(delta_d//abs(d[0]-d[1]))
    ntheta = int(delta_theta//abs(thetas[0]-thetas[1]))

    theta_start = thetas[:,1:ntheta+1]
    theta_end = thetas[:,-ntheta-1:-1]
    theta_ = np.append(thetas,theta_start+np.pi,axis=1)
    theta_ = np.append(theta_end-np.pi,theta_,axis=1)

    D_start = Ds[:,1:ntheta+1]
    D_end = Ds[:,-ntheta-1:-1]
    D_ = np.append(Ds,D_start,axis=1)
    D_ = np.append(D_end,D_,axis=1)

    HH_start = HH[:,1:ntheta+1]
    HH_end = HH[:,-ntheta-1:-1]
    HH_ = np.append(HH,HH_start,axis=1)
    HH_ = np.append(HH_end,HH_,axis=1)

    VV_start = VV[:,1:ntheta+1]
    VV_end = VV[:,-ntheta-1:-1]
    VV_ = np.append(VV,VV_start,axis=1)
    VV_ = np.append(VV_end,VV_,axis=1)

    chhvv = np.nan*np.ones_like(HH_).astype(np.complex)

    for i,θ in enumerate(theta_[0]):
        for j in range(len(D_[:,0])):
            if j < nd or j > len(D_[:,0])-nd or i < ntheta or i > len(theta_[0])-ntheta-1:
                continue
            else:
                chhvv[j,i] = np.corrcoef(HH_[j-nd:j+nd,i-ntheta:i+ntheta].flatten(),
                                         VV_[j-nd:j+nd,i-ntheta:i+ntheta].flatten())[1,0]

    chhvv = chhvv[:,ntheta:-ntheta]

# --------------------------------------------------------------------------------------------

def copolarized_phase_gradient(self):
    R = np.real(chhvv)
    I = np.imag(chhvv)

    Δd = 100
    Δθ = 20*np.pi/180.
    nd = int(Δd//abs(d[0]-d[1]))
    nθ = int(Δθ//abs(θs[0]-θs[1]))

    dRdz = np.nan*np.ones_like(chhvv).astype(np.float)
    dIdz = np.nan*np.ones_like(chhvv).astype(np.float)

    for i,θ in enumerate(θs):
        print('i:',i,'θ:',np.round(θ,2))
        for j in range(len(d)):
            if j < nd or j > len(d)-nd:
                continue
            else:
                Rs_ij = R[j-nd:j+nd,i]
                Is_ij = I[j-nd:j+nd,i]
                ds_ij = Ds[j-nd:j+nd,i]

                Ridxs = ~np.isnan(Rs_ij) & ~np.isnan(ds_ij)
                Iidxs = ~np.isnan(Is_ij) & ~np.isnan(ds_ij)

                if len(Rs_ij[Ridxs])<3 or len(Is_ij[Iidxs])<3:
                    dRdz[j,i] = np.nan
                    dIdz[j,i] = np.nan
                else:
                    dRdz[j,i] = np.polyfit(ds_ij[Ridxs],Rs_ij[Ridxs],1)[0]
                    dIdz[j,i] = np.polyfit(ds_ij[Iidxs],Is_ij[Iidxs],1)[0]

    dϕdz = (R*dIdz-I*dRdz)/(R**2.+I**2.)

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

def phase_gradient_to_cof(self,c = 300e6,fc = 200e6,Δϵ = 0.035,ϵ = 3.12):
    max_idx = np.argmax(self.dϕdz,axis=1)
    min_idx = np.argmin(self.dϕdz,axis=1)

    E2E1 = (c/(4.*np.pi*fc))*(2.*np.sqrt(ϵ)/Δϵ)*self.dϕdz[max_idx,np.arange(len(d))]
    return E2E1
