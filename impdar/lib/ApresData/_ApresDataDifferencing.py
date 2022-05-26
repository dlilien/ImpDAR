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
from ._ApresDataProcessing import phase2range
import numpy as np
from scipy.stats import linregress
from scipy.signal import medfilt, find_peaks


class ApresDiff():
    """
    Class for differencing between two Apres acquisitions.
    """
    def __init__(self,dat1,dat2):
        """Initialize differencieng fields

        dat1: class or string
            Apres data object to load (if string should be an impdar apres file)
        dat2: class or string
            Apres data object to load (if string should be an impdar apres file)
        """

        if type(dat1) == str and type(dat2) == str:
            dat1 = ApresData(dat1)
            dat2 = ApresData(dat2)
        self.fn1 = dat1.fn
        self.fn2 = dat2.fn

        if dat1.flags.range ==0 or dat2.flags.range ==0:
            raise TypeError('The range filter has not been executed on this data class, do that before proceeding.')

        # check that the two data objects are comparable
        if np.shape(dat1.data) != np.shape(dat2.data):
            raise TypeError('Acquisition inputs must be of the same shape.')
        if not np.all(abs(dat1.Rcoarse - dat2.Rcoarse) < 1e04):
            raise ValueError('Range vector should be the same for both acquisitions')

        # instantiate the differencing class
        self.data1 = dat1.data
        self.unc1 = dat1.uncertainty
        self.data2 = dat2.data
        self.unc2 = dat2.uncertainty
        self.Rcoarse = dat1.Rcoarse

        # take the flags and header from the first data object
        self.flags = dat1.flags
        self.header = dat1.header


    def phase_diff(self,win,step,Rcoarse=None):
        """
        Calculate the phase offset using a correlation coefficient.

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

        # Fill a depth array which will be more sparse than the full Rcoarse vector
        idxs = np.arange(win//2,(len(self.data1)-win//2),step)
        if Rcoarse is not None:
            self.ds = Rcoarse[idxs]
        else:
            self.ds = self.Rcoarse[idxs]

        # Create data and coherence vectors
        acq1 = self.data1
        acq2 = self.data2
        self.co = np.empty_like(self.ds).astype(complex)
        for i,idx in enumerate(idxs):
            # index two sub_arrays to compare
            arr1 = acq1[idx-win//2:idx+win//2]
            arr2 = acq2[idx-win//2:idx+win//2]
            # correlation coefficient between acquisitions
            # amplitude is coherence between acquisitions and phase is the offset
            self.co[i] = np.corrcoef(arr2,arr1)[1,0]

        self.flags.phase_diff = np.array([win,step])


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
        if self.flags.phase_diff is None:
            raise ValueError('Need to do the phase difference calculation first.')

        self.phi = np.angle(self.co).astype(float)
        for i in range(len(self.co)-1):
            idx = i+1
            if np.all(abs(self.co[idx-win:idx+win]) < thresh):
                continue
            if self.phi[idx]-self.phi[idx-1] > np.pi:
                self.phi[idx:] -= 2.*np.pi
            elif self.phi[idx]-self.phi[idx-1] < -np.pi:
                self.phi[idx:] += 2.*np.pi


    def phase2range(self,phi=None,lambdac=None,rc=None,K=None,ci=None):
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

        Output
        --------
        r: float or array
            range (m)

        ### Original Matlab File Notes ###
        Craig Stewart
        2014/6/10
        """

        if phi is None:
            phi = self.phi

        if lambdac is None:
            lambdac = self.header.lambdac

        if not all([K,ci]) or rc is None:
            # First order method
            # Brennan et al. (2014) eq 15
            self.r = lambdac*phi/(4.*np.pi)
        else:
            # Precise
            self.r = phi/((4.*np.pi/lambdac) - (4.*rc*K/ci**2.))


    def range_diff(self,uncertainty='CR'):
        """
        Convert the phase profile to range offset (vertical velocity)

        Parameters
        ---------
        self: class
            ApresData object
        uncertainty: string;
            default 'CR' Cramer-Rao bound as in Jordan et al. (2020)
        """

        # Check for unwrap
        if not hasattr(self,'phi'):
            raise ValueError('Should unwrap the phase profile before converting to range')

        win, step = self.flags.phase_diff

        # convert the phase offset to a distance vector
        self.w = self.phase2range(self.phi,
                self.header.lambdac,
                self.ds,
                self.header.chirp_grad,
                self.header.ci)

        if uncertainty == 'CR':
            # Error from Cramer-Rao bound, Jordan et al. (2020) Ann. Glac. eq. (5)
            sigma = (1./abs(self.co))*np.sqrt((1.-abs(self.co)**2.)/(2.*win))
            # convert the phase offset to a distance vector
            self.w_err = phase2range(sigma,
                    self.header.lambdac,
                    self.ds,
                    self.header.chirp_grad,
                    self.header.ci)

        elif uncertainty == 'noise_phasor':
            # Uncertainty from Noise Phasor as in Kingslake et al. (2014)
            # r_uncertainty should be calculated using the function phase_uncertainty defined in this script
            r_uncertainty = phase2range(self.unc1,self.header.lambdac) + phase2range(self.unc2,self.header.lambdac)
            idxs = np.arange(win//2,(len(self.data1)-win//2),step)
            self.w_err = np.array([np.nanmean(r_uncertainty[i-win//2:i+win//2]) for i in idxs])


    def strain_rate(self,strain_window=(200,1200),w_surf=0.):
        """
        Estimate the location of the ice-bed interface.

        Parameters
        ---------
        self: class
            ApresData object
        strain_window: tuple
            The range over which the vertical velocity profile is quasi linear and the average strain rate should be calculated.
        w_surf: float
            Vertical velocity at the surface (use ice equivalent accumulation rate)
            This is just to line up the profile, it is not really needed, can leave as 0.
        """

        if not hasattr(self,'w'):
            raise ValueError("Get the vertical velocity profile first with 'range_diff()'.")

        print('Calculating vertical strain rate over range from %s to %s'%strain_window)
        idx = np.logical_and(self.ds>strain_window[0],self.ds<strain_window[1])
        slope, intercept, r_value, p_value, std_err = linregress(self.ds[idx],self.w[idx])
        self.eps_zz = slope
        self.w0 = intercept
        print('Vertical strain rate:',self.eps_zz)
        print('r_squared:',r_value**2.)

        self.w += w_surf - self.w0


    def bed_pick(self,sample_threshold=50,coherence_threshold=0.9,
                 filt_kernel=201,prominence=10,peak_width=300):
        """
        Estimate the location of the ice-bed interface.

        Parameters
        ---------
        self: class
            ApresData object
        sample_threshold: int
            Number of samples to tolerate for difference between the pick from acquisition 1 and 2
        coherence_threshold: int
            Minimum coherence allowed for the picked bed
        filt_kernel: int
            Kernel for median filter
        prominence: int
            How high the bed power needs to be above neighboring power profile
        peak_width: int
            Width of the power peak
        """

        P1 = 10.*np.log10(self.data1**2.)
        mfilt1 = medfilt(P1.astype(float),filt_kernel)
        peaks1 = find_peaks(mfilt1,prominence=prominence,width=peak_width)[0]
        bed_idx1 = max(peaks1)

        P2 = 10.*np.log10(self.data2**2.)
        mfilt2 = medfilt(P2.astype(float),filt_kernel)
        peaks2 = find_peaks(mfilt2,prominence=prominence,width=peak_width)[0]
        bed_idx2 = max(peaks2)

        if not abs(bed_idx1 - bed_idx2) < sample_threshold:
            raise ValueError('Bed pick from first and second acquisitions are too far apart.')

        bed_samp = (bed_idx1+bed_idx2)//2
        bed_power = (mfilt1[bed_idx1]+mfilt2[bed_idx2])/2.
        bed_range = self.Rcoarse[bed_samp]

        diff_idx = np.argmin(abs(self.ds-bed_range))
        bed_coherence= np.median(abs(self.co[diff_idx-10:diff_idx+10]))

        if not bed_coherence > coherence_threshold:
            raise ValueError('Bed pick has too low coherence.')

        self.bed = np.array([bed_samp,bed_range,bed_coherence,bed_power])
