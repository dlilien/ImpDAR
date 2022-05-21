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

from . import ApresData
import numpy as np

class ApresDiff():
    """
    Class for differencing between two Apres acquisitions.
    """
    def __init__(self,dat1,dat2):
        """Initialize differencieng fields

        self: class
            Apres data object
        fn_secondary: string
            file name for secondary data object
        """

        if type(dat1) == str and type(dat2) == str:
            dat1 = ApresData(dat1)
            dat2 = ApresData(dat2)

        self.fn1 = dat1.fn
        self.fn2 = dat2.fn

        if not np.all([dat1.flags.range,dat2.flags.range]):
            raise ValueError('Need to do the range conversion on both data objects first.')

        # check that the two data objects are comparable
        if np.shape(dat1.data) != np.shape(dat2.data):
            raise TypeError('Acquisition inputs must be of the same shape.')
        if not np.all(abs(dat1.Rcoarse - dat2.Rcoarse) < 1e04):
            raise ValueError('Range vector should be the same for both acquisitions')

        # instantiate the differencing class
        self.data1 = dat1.data
        self.data2 = dat2.data
        self.Rcoarse = dat1.Rcoarse

        # take the header from the first data object
        self.header = dat1.header


    def phase_diff(self,win,step,Rcoarse=None):
        """
        Calculate the vertical motion using a correlation coefficient.

        Parameters
        ---------
        self: class
            Apres differencing object
        win: int
            window size over which to do the correlation coefficient calculation
        step: int
            step size for the window to move between calculations
        Rcoarse: array; optional
            if an external depth array is desired, input here
        """

        idxs = np.arange(win//2,(len(self.data1)-win//2),step)
        if Rcoarse is not None:
            self.ds = Rcoarse[idxs]
        else:
            self.ds = self.Rcoarse[idxs]
        self.co = np.empty_like(self.ds).astype(complex)
        acq1 = self.data1
        acq2 = self.data2
        for i,idx in enumerate(idxs):
            # index two sub_arrays to compare
            arr1 = acq1[idx-win//2:idx+win//2]
            arr2 = acq2[idx-win//2:idx+win//2]
            # correlation coefficient to get the motion
            # the amplitude indicates how well the reflections match between acquisitions
            # the phase is a measure of the offset
            self.co[i] = np.corrcoef(arr1,arr2)[1,0]


    def range_diff(self,r_uncertainty=None,uncertainty='CR'):
        """
        r_uncertainty: array; optional
            if unceratinty based on the noise vector is desired input value here
            this should be the sum of uncertainty from both acquisitions.
        uncertainty: string;
            default 'CR' Cramer-Rao bound as in Jordan et al. (2020)

        range_diff: array
            vertical motion in meters

        """

        # convert the phase offset to a distance vector
        self.w = phase2range(self.phi,
                self.header.lambdac,
                ds,
                self.header.chirp_grad,
                self.header.ci)

        if uncertainty == 'CR':
            # Error from Cramer-Rao bound, Jordan et al. (2020) Ann. Glac. eq. (5)
            sigma = (1./abs(self.co))*np.sqrt((1.-abs(self.co)**2.)/(2.*win))
            # convert the phase offset to a distance vector
            self.w_err = phase2range(sigma,
                    self.header.lambdac,
                    ds,
                    self.header.chirp_grad,
                    self.header.ci)

        elif uncertainty == 'noise_phasor':
            # Uncertainty from Noise Phasor as in Kingslake et al. (2014)
            # r_uncertainty should be calculated using the function phase_uncertainty defined in this script
            self.w_err = np.array([np.nanmean(r_uncertainty[i-win//2:i+win//2]) for i in idxs])


    def phase_unwrap(self,win=10,thresh=0.9):
        """
        Unwrap the phase profile to get a

        Parameters
        ---------
        self: class
            ApresData object
        win: int
            window; number of points to include in the wrap-finding window
        thresh: float (0-1)
            threshold; coherence threshold where to stop looking for phase wraps
        """

        # Check that the range difference has been done

        self.phi = np.angle(self.co).astype(float)

        for i in range(len(self.co)-1):
            idx = i+1
            if np.all(abs(self.co[idx-win:idx+win]) < thresh):
                continue
            if self.phi[idx]-self.phi[idx-1] > np.pi:
                self.phi[idx:] -= 2.*np.pi
            elif self.phi[idx]-self.phi[idx-1] < -np.pi:
                self.phi[idx:] += 2.*np.pi


    def strain_rate(self,strain_window=np.array([1000,3000])):
        """
        Estimate the location of the ice-bed interface.

        Parameters
        ---------
        self: class
            ApresData object
        strain_window: 2x1 np.array
            The range over which the bed is predicted to be
        """


    def bed_pick(self,bed_window=np.array([1000,3000])):
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
