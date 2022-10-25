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

from ._ApresDataProcessing import phase2range
import numpy as np
from scipy.stats import linregress
from scipy.signal import medfilt, find_peaks


def coherence(s1, s2):
    """
    Phase correlation between two elements of the scattering matrix
    Jodan et al. (2019) eq. 13

    Parameters
    ---------
    s1: array
        first acquisition
    s2:
        second acquisition

    Output
    ---------
    c:  array
        phase coherence
    """

    if hasattr(s1, '__len__') and hasattr(s2, '__len__'):
        top = np.sum(np.dot(s1, np.conj(s2)))
        bottom = np.sqrt(np.sum(np.abs(s1)**2.)*np.sum(np.abs(s2)**2.))
        c = top/bottom
    else:
        top = np.dot(s1, np.conj(s2))
        bottom = np.sqrt(np.abs(s1)**2.*np.abs(s2)**2.)
        c = top/bottom

    return c


def phase_diff(self, win, step, range_ext=None):
    """
    Calculate the phase offset along the full ApRES acquisition
    using a correlation coefficient within a moving window.

    Parameters
    ---------
    self: class
        Apres differencing object
    win: int
        window size over which to do the correlation coefficient calculation
    step: int
        step size for the window to move between calculations
    range_ext: array; optional
        if an external depth array is desired, input here
    """

    # Fill a depth array which will be more sparse than the full Range vector
    idxs = np.arange(win//2, len(self.data)-win//2, step).astype(int)
    if range_ext is not None:
        self.ds = range_ext[idxs]
    else:
        self.ds = self.range[idxs]

    # Create data and coherence vectors
    acq1 = self.data
    acq2 = self.data2
    self.co = np.empty_like(self.ds).astype(np.cdouble)
    for i, idx in enumerate(idxs):
        # index two sub_arrays to compare
        arr1 = acq1[idx-win//2:idx+win//2]
        arr2 = acq2[idx-win//2:idx+win//2]
        # correlation coefficient between acquisitions
        # amplitude is coherence between acquisitions and phase is the offset
        self.co[i] = coherence(arr1, arr2)

    self.flags.phase_diff = np.array([win, step])


def phase_unwrap(self, win=10, thresh=0.9):
    """
    Unwrap the phase profile to get one that is
    either monotonically increasing or monotonically decreasing.

    Parameters
    ---------
    self: class
        ApresTimeDiff object
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


def range_diff(self, uncertainty='noise_phasor'):
    """
    Convert the phase profile to range offset (vertical velocity).

    Parameters
    ---------
    self: class
        ApresTimeDiff object
    uncertainty: string;
        default 'noise_phasor' as in Kingslake et al. (2014)
    """

    # Check for unwrap
    if not hasattr(self, 'phi'):
        raise ValueError('Should unwrap the phase profile before converting to range')

    win, step = self.flags.phase_diff

    # convert the phase offset to a distance vector
    self.w = phase2range(self, self.phi,
                         self.header.lambdac,
                         self.ds,
                         self.header.chirp_grad,
                         self.header.ci)

    # If the individual acquisitions have had uncertainty calculations
    if self.unc1 is not None:

        if uncertainty == 'CR':
            # Error from Cramer-Rao bound, Jordan et al. (2020) Ann. Glac. eq. (5)
            sigma = (1./abs(self.co))*np.sqrt((1.-abs(self.co)**2.)/(2.*win))
            # convert the phase offset to a distance vector
            self.w_err = phase2range(self, sigma,
                                     self.header.lambdac,
                                     self.ds,
                                     self.header.chirp_grad,
                                     self.header.ci)

        elif uncertainty == 'noise_phasor':
            # Uncertainty from Noise Phasor as in Kingslake et al. (2014)
            # r_uncertainty should be calculated using the function phase_uncertainty defined in this script
            r_uncertainty = phase2range(self, self.unc1, self.header.lambdac) +\
                phase2range(self, self.unc2, self.header.lambdac)
            idxs = np.arange(win//2, len(self.data)-win//2, step)
            self.w_err = np.array([np.nanmean(r_uncertainty[i-win//2:i+win//2]) for i in idxs])


def strain_rate(self, strain_window=(200, 1200), w_surf=0.):
    """
    Estimate the average vertical strain rate within some
    range span provided.

    Parameters
    ---------
    self: class
        ApresTimeDiff object
    strain_window: tuple
        The range over which the vertical velocity profile is quasi linear and the average strain rate should be calculated.
    w_surf: float
        Vertical velocity at the surface (use ice equivalent accumulation rate)
        This is just to line up the profile, it is not really needed, can leave as 0.
    """

    if not hasattr(self, 'w'):
        raise ValueError("Get the vertical velocity profile first with 'range_diff()'.")

    print('Calculating vertical strain rate over range from %s to %s meters.' % strain_window)
    idx = np.logical_and(self.ds > strain_window[0], self.ds < strain_window[1])
    slope, intercept, r_value, p_value, std_err = linregress(self.ds[idx], self.w[idx])
    self.eps_zz = slope
    self.w0 = intercept
    print('Vertical strain rate (yr-1):', self.eps_zz)
    print('r_squared:', r_value**2.)

    self.w += w_surf - self.w0


def bed_pick(self, sample_threshold=50, coherence_threshold=0.9,
             filt_kernel=201, prominence=10, peak_width=300):
    """
    Estimate the location of the ice-bed interface.

    Parameters
    ---------
    self: class
        ApresTimeDiff object
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

    P1 = 10.*np.log10(self.data**2.)
    mfilt1 = medfilt(P1.real, filt_kernel)
    peaks1 = find_peaks(mfilt1, prominence=prominence, width=peak_width)[0]
    bed_idx1 = max(peaks1)

    P2 = 10.*np.log10(self.data2**2.)
    mfilt2 = medfilt(P2.real, filt_kernel)
    peaks2 = find_peaks(mfilt2, prominence=prominence, width=peak_width)[0]
    bed_idx2 = max(peaks2)

    if not abs(bed_idx1 - bed_idx2) < sample_threshold:
        raise ValueError('Bed pick from first and second acquisitions are too far apart.')

    bed_samp = (bed_idx1+bed_idx2)//2
    bed_power = (mfilt1[bed_idx1]+mfilt2[bed_idx2])/2.
    bed_range = self.range[bed_samp]

    diff_idx = np.argmin(abs(self.ds-bed_range))
    bed_coherence = np.median(abs(self.co[diff_idx-10:diff_idx+10]))

    if not bed_coherence > coherence_threshold:
        raise ValueError('Bed pick has too low coherence.')

    self.bed = np.array([bed_samp, bed_range, bed_coherence, bed_power])
