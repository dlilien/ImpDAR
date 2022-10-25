#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Benjamin Hills <bhills@uw.edu>
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


def apres_range(self, p, max_range=4000, winfun='blackman'):
    """
    Range conversion.

    Parameters
    ---------
    self: class
        ApresData object
    p: int
        pad factor, level of interpolation for fft
    max_range: float
        cut off after some maximum range
    winfun: str
        window function for fft


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
    # window for fft
    if winfun not in ['blackman', 'bartlett', 'hamming', 'hanning', 'kaiser']:
        raise TypeError('Window must be in: blackman, bartlett, hamming, hanning, kaiser')
    elif winfun == 'blackman':
        win = np.blackman(self.snum)
    elif winfun == 'bartlett':
        win = np.bartlett(self.snum)
    elif winfun == 'hamming':
        win = np.hamming(self.snum)
    elif winfun == 'hanning':
        win = np.hanning(self.snum)
    elif winfun == 'kaiser':
        win = np.kaiser(self.snum)

    # round-trip delay Brennan et al. (2014) eq. 18
    tau = np.arange(nf)/(self.header.bandwidth*p)

    # Get the coarse range
    self.Rcoarse = tau*self.header.ci/2.

    # Calculate phase of each range bin center for correction
    # Brennan et al. (2014) eq. 17 measured at t=T/2
    self.phiref = 2.*np.pi*self.header.fc*tau - (self.header.chirp_grad*tau**2.)/2

    # --- Loop through for each chirp in burst --- #

    # pre-allocate
    spec = np.zeros((self.bnum, self.cnum, nf)).astype(np.cdouble)
    spec_cor = np.zeros((self.bnum, self.cnum, nf)).astype(np.cdouble)

    for ib in range(self.bnum):
        for ic in range(self.cnum):
            # isolate the chirp and preprocess before transform
            chirp = self.data[ib, ic, :].copy()
            chirp = chirp-np.mean(chirp)  # de-mean
            chirp *= win  # windowed

            # fourier transform
            fft_chirp = (np.sqrt(2.*p)/len(chirp))*np.fft.fft(chirp, p*self.snum)  # fft and scale for padding
            fft_chirp /= np.sqrt(np.mean(win**2.))  # scale with rms of window

            # output
            spec[ib, ic, :] = fft_chirp[:nf]  # positive frequency half of spectrum up to (nyquist minus deltaf)
            comp = np.exp(-1j*(self.phiref))  # unit phasor with conjugate of phiref phase
            spec_cor[ib, ic, :] = comp*fft_chirp[:nf] # positive frequency half of spectrum with ref phase subtracted

    self.data = spec_cor.copy()
    self.spec = spec.copy()
    self.data_dtype = self.data.dtype

    # precise range measurement
    self.Rfine = phase2range(self, np.angle(self.data), self.header.lambdac,
                             np.tile(self.Rcoarse, (self.bnum, self.cnum, 1)),
                             self.header.chirp_grad, self.header.ci)

    # Crop output variables to useful depth range only
    n = np.argmin(self.Rcoarse <= max_range)
    self.Rcoarse = self.Rcoarse[:n]
    self.Rfine = self.Rfine[:n]
    self.data = self.data[:, :, :n]
    self.spec = self.spec[:, :, :n]
    self.snum = n

    self.flags.range = max_range


def phase_uncertainty(self, bed_range):
    """
    Calculate the phase uncertainty using a noise phasor.

    Following Kingslake et al. (2014)

    Parameters
    ---------
    self: class
        ApresData object
    bed_range: float
        Range to the ice-bed interface, so noise floor can be calculated beneath it
    """

    if self.flags.range == 0:
        raise TypeError('The range filter has not been executed on this data class, do that before the uncertainty calculation.')

    # Get measured phasor from the data class, and use the median magnitude for noise phasor
    meas_phasor = np.squeeze(self.data)
    median_mag = np.nanmedian(abs(meas_phasor[np.argwhere(self.Rcoarse>bed_range)]))
    # Noise phasor with random phase and magnitude equal to median of measured phasor
    noise_phase = np.random.uniform(-np.pi,np.pi,np.shape(meas_phasor))
    noise_phasor = median_mag*(np.cos(noise_phase)+1j*np.sin(noise_phase))
    noise_orth = median_mag*np.sin(np.angle(meas_phasor)-np.angle(noise_phasor))
    # Phase uncertainty is the deviation in the phase introduced by the noise phasor when it is oriented perpendicular to the reflector phasor
    self.uncertainty = np.abs(np.arcsin(noise_orth/np.abs(meas_phasor)))

    self.flags.uncertainty = True


def phase2range(self, phi, lambdac=None, rc=None, K=None, ci=None):
    """
    Convert phase difference to range for FMCW radar

    Parameters
    ---------
    lambdac: float
        wavelength (m) at center frequency
    rc: float; optional
        coarse range of bin center (m)
    K:  float; optional
        chirp gradient (rad/s/s)
    ci: float; optional
        propagation velocity (m/s)


    ### Original Matlab File Notes ###
    Craig Stewart
    2014/6/10
    """

    if lambdac is None:
        lambdac = self.header.lambdac

    if not all([K,ci]) or rc is None:
        # First order method
        # Brennan et al. (2014) eq 15
        r = lambdac*phi/(4.*np.pi)
    else:
        # Precise
        r = phi/((4.*np.pi/lambdac) - (4.*rc*K/ci**2.))

    return r


def stacking(self, num_chirps=None):
    """
    Stack traces/chirps together to beat down the noise.
    Can stack across bursts if multiple are present.

    Parameters
    ---------
    self: class
        ApresData object
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

