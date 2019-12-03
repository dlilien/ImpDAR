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

def apres_range(self,p,max_range=4000,winfun='blackman'):
    """

    Parameters
    ---------
    self: class
        data object
    p: int
        pad factor, level of interpolation for fft
    winfun: str
        window function for fft

    Output
    --------
    Rcoarse: array
        range to bin centres (m)
    Rfine: array
        range to reflector from bin centre (m)
    spec_cor: array
        spectrum corrected. positive frequency half of spectrum with
        ref phase subtracted. This is the complex signal which can be used for
        cross-correlating two shot segements.



    ### Original Matlab File Notes ###
    Phase sensitive processing of FMCW radar data based on Brennan et al. 2013

    Based on Paul's scripts but following the nomenclature of:
    "Phase-sensitive FMCW radar imaging system for high precision Antarctic
    ice shelf profile monitoring"
    Brennan, Lok, Nicholls and Corr, 2013

    Summary: converts raw FMCW radar voltages into a range for

    Craig Stewart
    2013 April 24
    Modified frequencies 10 April 2014
    """

    if self.flags.range != 0:
        raise TypeError('The range filter has already been done on these data.')

    # Processing settings
    nf = int(np.floor(p*self.snum/2))    # number of frequencies to recover
    # TODO: implement all the other window functions
    if winfun not in ['blackman','bartlett','hamming','hanning','kaiser']:
        raise TypeError('Window must be in: blackman, bartlett, hamming, hanning, kaiser')
    elif winfun == 'blackman':
        win = np.blackman(self.snum)

    # Get the coarse range
    self.Rcoarse = np.arange(nf)*self.header.ci/(2*self.header.bandwidth*p)

    # Calculate phase of each range bin centre for correction
    # eq 17: phase for each range bin centre (measured at t=T/2), given that tau = n/(B*p)
    self.phiref = 2.*np.pi*self.header.fc*np.arange(nf)/(self.header.bandwidth*p) - \
            (self.header.chirp_grad*np.arange(nf)**2)/(2*self.header.bandwidth**2*p**2)

    # Measure the sampled IF signal: FFT to measure frequency and phase of IF
    #deltaf = 1/(T*p); # frequency step of FFT
    #f = [0:deltaf:fs/2-deltaf]; # frequencies measured by the fft - changed 16 April 2014, was #f = [0:deltaf:fs/2];
    #Rcoarse = f*ci*T/(2*B); # Range at the centre of each range bin: eq 14 (rearranged) (p is accounted for inf)
    #Rcoarse = [0:1/p:T*fs/2-1/p]*ci/(2*B); # Range at the centre of each range bin: eq 14 (rearranged) (p is accounted for inf)

    # --- Loop through for each chirp in burst --- #

    # preallocate
    spec = np.zeros((self.bnum,self.cnum,nf)).astype(np.cdouble)
    spec_cor = np.zeros((self.bnum,self.cnum,nf)).astype(np.cdouble)

    for ib in range(self.bnum):
        for ic in range(self.cnum):
            # isolate the chirp and preprocess before transform
            chirp = self.data[ib,ic,:].copy()
            chirp = chirp-np.mean(chirp) # de-mean
            chirp *= win # windowed

            # fourier transform
            fft_chirp = (np.sqrt(2.*p)/len(chirp))*np.fft.fft(chirp,p*self.snum) # fft and scale for padding
            fft_chirp /= np.sqrt(np.mean(win**2.)) # scale with rms of window

            # output
            spec[ib,ic,:] = fft_chirp[:nf] # positive frequency half of spectrum up to (nyquist minus deltaf)
            comp = np.exp(-1j*(self.phiref)) # unit phasor with conjugate of phiref phase
            spec_cor[ib,ic,:] = comp*fft_chirp[:nf] # positive frequency half of spectrum with ref phase subtracted

    self.data = spec_cor.copy()
    self.spec = spec.copy()

    self.Rfine = phase2range(np.angle(self.data),self.header.lambdac,
            np.tile(self.Rcoarse,(self.bnum,self.cnum,1)),
            self.header.chirp_grad,self.header.ci)

    # Crop output variables to useful depth range only
    n = np.argmin(self.Rcoarse<=max_range)
    self.Rcoarse = self.Rcoarse[:n]
    self.Rfine = self.Rfine[:,:,:n]
    self.data = self.data[:,:,:n]
    self.spec = self.spec[:,:,:n]
    self.snum = n

    self.flags.range = max_range

# --------------------------------------------------------------------------------------------

def phase2range(phi,lambdac,rc=None,K=None,ci=None):
    """
    Convert phase difference to range for FMCW radar

    Parameters
    ---------
    phi: float or array
        phase (radians), must be of spectrum after bin center correction
    lambdac: float
        wavelength (m) at center frequency
    rc: float; optional
        coarse range of bin center (m)
    K:  float; optional
        chirp gradient (rad/s/s)
    ci: float; optional
        propagation velocity (m/s)

    Output
    --------
    r: float or array
        range (m)

    ### Original Matlab File Notes ###
    Craig Stewart
    2014/6/10
    """

    if not all([K,ci]) or rc is None:
        # First order method
        r = lambdac*phi/(4.*np.pi)
    else:
        # Precise
        r = phi/((4.*np.pi/lambdac) - (4.*rc*K/ci**2.))
    return r

# --------------------------------------------------------------------------------------------

def range_diff(self,acq1,acq2,win,step):
    """
    Calculate the vertical motion using a correlation coefficient.

    Parameters
    ---------
    self: class
        data object
    acq1: array
        first acquisition for comparison
    acq2: array
        second acquisition for comparison
    win: int
        window size over which to do the correlation coefficient calculation
    step: int
        step size for the window to move between calculations

    Output
    --------
    ds: array
        depths at which the correlation coefficient is calculated
    phase_diff: array
        correlation coefficient between acquisitions
        amplitude indicates how well reflection packets match between acquisitions
        phase is a measure of the vertical motion
    range_diff: array
        vertical motion in meters
    """

    if np.shape(acq1) != np.shape(acq2):
        raise TypeError('Acquisition inputs must be of the same shape.')

    idxs = np.arange(0,(len(acq1)-win),step)
    ds = self.Rcoarse[idxs]
    phase_diff = np.empty_like(ds)
    for i,idx in enumerate(idxs):
        # index two sub_arrays to compare
        arr1 = acq1[idx:idx+win]
        arr2 = acq2[idx:idx+win]
        # correlation coefficient to get the motion
        # the amplitude indicates how well the reflections match between acquisitions
        # the phase is a measure of the offset
        phase_diff[i] = np.corrcoef(arr1,arr2)[1,0]

    # convert the phase offset to a distance vector
    range_diff = phase2range(np.angle(phase_diff),
            self.header.lambdac,
            self.Rcoarse,
            self.header.chirp_grad,
            self.header.ci)

    return ds, phase_diff, range_diff

# --------------------------------------------------------------------------------------------

def stacking(self,num_chirps=None):
    """
    Stack traces/chirps together to beat down the noise.

    Parameters
    ---------
    num_chirps: int
        number of chirps to average over
    """

    if num_chirps == None:
        num_chirps = self.cnum*self.bnum

    num_chirps = int(num_chirps)

    if num_chirps == self.cnum:
        self.data = np.reshape(np.mean(self.data,axis=1),(self.bnum,1,self.snum))
        self.cnum = 1
    else:
        # reshape to jump across bursts
        data_hold = np.reshape(self.data,(1,self.cnum*self.bnum,self.snum))
        # take only the first set of chirps
        data_hold = data_hold[:,:num_chirps,:]

        self.data = np.array([np.mean(data_hold,axis=1)])
        self.bnum = 1
        self.cnum = 1

    self.flags.stack = num_chirps
