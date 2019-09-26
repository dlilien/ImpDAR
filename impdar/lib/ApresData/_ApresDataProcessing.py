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

def apres_range(self,p,max_range=2000,winfun='blackman'):
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
    spec = np.zeros((nf,self.cnum)).astype(np.cdouble)
    spec_cor = np.zeros((nf,self.cnum)).astype(np.cdouble)

    for ii in range(self.cnum):
        # isolate the chirp and preprocess before transform
        chirp = self.data[:,ii].copy()
        chirp = chirp-np.mean(chirp) # de-mean
        chirp *= win # windowed

        # fourier transform
        fft_chirp = (np.sqrt(2.*p)/len(chirp))*np.fft.fft(chirp,p*self.snum) # fft and scale for padding
        fft_chirp /= np.sqrt(np.mean(win**2.)) # scale with rms of window

        # output
        spec[:,ii] = fft_chirp[:nf] # positive frequency half of spectrum up to (nyquist minus deltaf)
        comp = np.exp(-1j*(self.phiref)) # unit phasor with conjugate of phiref phase
        spec_cor[:,ii] = comp*fft_chirp[:nf] # positive frequency half of spectrum with ref phase subtracted

    self.data = spec_cor.copy()

    self.Rfine = phase2range(np.angle(self.data),self.header.lambdac,
            np.transpose(np.tile(self.Rcoarse,(self.cnum,1))),
            self.header.chirp_grad,self.header.ci)

    # Crop output variables to useful depth range only
    n = np.argmin(self.Rcoarse<=max_range)
    self.Rcoarse = self.Rcoarse[:n]
    self.Rfine = self.Rfine[:n,:]
    self.data = self.data[:n,:]
    self.snum = n

    self.flags.range = max_range

# --------------------------------------------------------------------------------------------

def phase2range(phi,lambdac,rc=None,K=None,ci=None):
    """
    # r = fmcw_phase2range(phi,lambdac,rc,K,ci)
    #
    # Convert phase difference to range for FMCW radar
    #
    # args:
    # phi: phase (radians), must be of spectrum after bin centre correction
    # lambdac: wavelength (m) at centre frequency
    #
    # optional args: (used for precise method)
    # rc: coarse range of bin centre (m)
    # K = chirp gradient (rad/s/s)
    # ci = propagation velocity (m/s)
    #
    # Craig Stewart
    # 2014/6/10
    """

    if not all([K,ci]) or rc is None:
        # First order method
        r = lambdac*phi/(4.*np.pi)
    else:
        # Precise
        r = phi/((4.*np.pi/lambdac) - (4.*rc*K/ci**2.))

    return r

# --------------------------------------------------------------------------------------------

def stacking(self,num_chirps=None):
    """
    Parameters
    ---------
    num_chirps: int
        number of chirps to average
    """

    #TODO: update for stacking across multiple bursts

    if self.flags.range == 0:
        raise TypeError('Do the range conversion before stacking chirps.')

    if num_chirps == None:
        num_chirps = self.cnum
    elif num_chirps > self.cnum:
        Warning('Number of chirps given for average is greater than the number of chirps in the burst.\
                reducing to self.cnum.')
        num_chirps = self.cnum

    if num_chirps == self.cnum:
        self.stacked_data = np.mean(np.real(self.data),axis=1)
    else:
        data_hold = np.empty((self.cnum//num_chirps,self.snum))
        for i in range(self.cnum//num_chirps):
            data_hold[i,:] = np.mean(np.real(self.data[:,i*num_chirps:(i+1)*num_chirps]),axis=1)
        self.stacked_data = data_hold

    self.flags.stack = num_chirps
