#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 Benjamin Hills <bhills@uw.edu>
#
# Distributed under terms of the GNU GPL3 license.

"""
Differencing between two ApRES data objects

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

May 20 2022
"""

import numpy as np

class ApresDiff():


def bed_picker(self,bed_window=np.array([1000,3000])):
    """
    Estimate the location of the ice-bed interface.

    Parameters
    ---------
    self: class
        ApresData object
    bed_window: 2x1 np.array
        The range over which the bed is predicted to be
    """
    self.bed_pick = np.array([bed_samp,bed_range,bed_power])

def range_diff(self,acq1,acq2,win,step,Rcoarse=None,r_uncertainty=None,uncertainty='CR'):
    """
    Calculate the vertical motion using a correlation coefficient.

    Parameters
    ---------
    self: class
        data object
    acq1: array
        first acquisition for comparison (dat1.data)
    acq2: array
        second acquisition for comparison (dat2.data)
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

