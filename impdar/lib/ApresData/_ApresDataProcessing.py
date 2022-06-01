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

def apres_range(self,p,max_range=4000,winfun='blackman'):
    """
    Range conversion.

    Parameters
    ---------
    self: class
        ApresData object
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
    # window for fft
    if winfun not in ['blackman','bartlett','hamming','hanning','kaiser']:
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
    self.phiref = 2.*np.pi*self.header.fc*tau -(self.header.chirp_grad*tau**2.)/2

    # --- Loop through for each chirp in burst --- #

    # pre-allocate
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
    self.data_dtype = self.data.dtype

    # precise range measurement
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


def phase_uncertainty(self,bed_range):
    """
    Calculate the phase uncertainty using a noise phasor.

    Following Kingslake et al. (2014)

    Parameters
    ---------
    self: class
        ApresData object

    Returns
    --------
    phase_uncertainty: array
        uncertainty in the phase (rad)
    r_uncertainty: array
        uncertainty in the range (m) calculated from phase uncertainty
    """

    if self.flags.range == 0:
        raise TypeError('The range filter has not been executed on this data class, do that before the uncertainty calculation.')

    # Get measured phasor from the data class, and use the median magnitude for noise phasor
    meas_phasor = self.data
    median_mag = np.nanmedian(abs(meas_phasor[:,:,np.argwhere(self.Rcoarse>bed_range)]))
    # Noise phasor with random phase and magnitude equal to median of measured phasor
    noise_phase = np.random.uniform(-np.pi,np.pi,np.shape(meas_phasor))
    noise_phasor = median_mag*(np.cos(noise_phase)+1j*np.sin(noise_phase))
    noise_orth = median_mag*np.sin(np.angle(meas_phasor)-np.angle(noise_phasor))
    # Phase uncertainty is the deviation in the phase introduced by the noise phasor when it is oriented perpendicular to the reflector phasor
    phase_uncertainty = np.abs(np.arcsin(noise_orth/np.abs(meas_phasor)))
    # Convert phase to range
    r_uncertainty = phase2range(phase_uncertainty,
            self.header.lambdac,
            self.Rcoarse,
            self.header.chirp_grad,
            self.header.ci)

    return phase_uncertainty, r_uncertainty


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

    Returns
    --------
    r: float or array
        range (m)

    ### Original Matlab File Notes ###
    Craig Stewart
    2014/6/10
    """

    if not all([K,ci]) or rc is None:
        # First order method
        # Brennan et al. (2014) eq 15
        r = lambdac*phi/(4.*np.pi)
    else:
        # Precise
        r = phi/((4.*np.pi/lambdac) - (4.*rc*K/ci**2.))
    return r


def range_diff(self,acq1,acq2,win,step,Rcoarse=None,r_uncertainty=None,uncertainty='CR'):
    """
    Calculate the vertical motion using a correlation coefficient.

    Parameters
    ---------
    self: class
        ApresData object
    acq1: array
        first acquisition for comparison
    acq2: array
        second acquisition for comparison
    win: int
        window size over which to do the correlation coefficient calculation
    step: int
        step size for the window to move between calculations
    Rcoarse: array; optional
        if an external depth array is desired, input here
    r_uncertainty: array; optional
        if unceratinty based on the noise vector is desired input value here
        this should be the sum of uncertainty from both acquisitions.
    uncertainty: string;
        default 'CR' Cramer-Rao bound as in Jordan et al. (2020)

    Returns
    --------
    ds: array
        depths at which the correlation coefficient is calculated
    co: array
        correlation coefficient between acquisitions
        amplitude indicates how well reflection packets match between acquisitions
        phase is a measure of the vertical motion
    range_diff: array
        vertical motion in meters
    range_diff_unc: array
        uncertainty of vertical motion in meters
    """

    if np.shape(acq1) != np.shape(acq2):
        raise TypeError('Acquisition inputs must be of the same shape.')

    idxs = np.arange(win//2,(len(acq1)-win//2),step)
    if Rcoarse is not None:
        ds = Rcoarse[idxs]
    else:
        ds = self.Rcoarse[idxs]
    co = np.empty_like(ds).astype(np.complex)
    for i,idx in enumerate(idxs):
        # index two sub_arrays to compare
        arr1 = acq1[idx-win//2:idx+win//2]
        arr2 = acq2[idx-win//2:idx+win//2]
        # correlation coefficient to get the motion
        # the amplitude indicates how well the reflections match between acquisitions
        # the phase is a measure of the offset
        co[i] = np.corrcoef(arr1,arr2)[1,0]

    # convert the phase offset to a distance vector
    r_diff = phase2range(np.angle(co),
            self.header.lambdac,
            ds,
            self.header.chirp_grad,
            self.header.ci)

    if uncertainty == 'CR':
        # Error from Cramer-Rao bound, Jordan et al. (2020) Ann. Glac. eq. (5)
        sigma = (1./abs(co))*np.sqrt((1.-abs(co)**2.)/(2.*win))
        # convert the phase offset to a distance vector
        r_diff_unc = phase2range(sigma,
                self.header.lambdac,
                ds,
                self.header.chirp_grad,
                self.header.ci)

    elif uncertainty == 'noise_phasor':
        # Uncertainty from Noise Phasor as in Kingslake et al. (2014)
        # r_uncertainty should be calculated using the function phase_uncertainty defined in this script
        r_diff_unc = np.array([np.nanmean(r_uncertainty[i-win//2:i+win//2]) for i in idxs])

    return ds, co, r_diff, r_diff_unc


def stacking(self,num_chirps=None):
    """
    Stack traces/chirps together to beat down the noise.

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
