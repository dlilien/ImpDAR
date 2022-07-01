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

import os
import h5py
from . import ApresData, ApresFlags, ApresHeader
from ._ApresDataProcessing import phase2range
from ..ImpdarError import ImpdarError
import numpy as np
from scipy.io import loadmat
from scipy.stats import linregress
from scipy.signal import medfilt, find_peaks


class ApresDiff():
    """
    Class for differencing between two Apres acquisitions.

    We keep track of processing steps with the flags attribute.
    This base version's __init__ takes two filenames of .mat files in the old StODeep format to load.
    """
    #: Attributes that every ApresData object should have and should not be None.
    attrs_guaranteed = ['data',
                        'data2',
                        'Rcoarse',
                        'fn1',
                        'fn2',
                        'fn']

    #: Optional attributes that may be None without affecting processing.
    #: These may not have existed in old StoDeep files that we are compatible with,
    #: and they often cannot be set at the initial data load.
    #: If they exist, they all have units of meters.
    attrs_optional = ['unc1',
                      'unc2',
                      'ds',
                      'co',
                      'phi',
                      'w',
                      'w_err',
                      'w_0',
                      'u',
                      'eps_zz',
                      'bed']

    from ._ApresDataSaving import save

    def __init__(self,fn,dat2=None):
        """Initialize differencing fields

        fn: class or string
            Apres data object to load (if string should be an impdar apres file)
        dat2: class or string
            Apres data object to load (if string should be an impdar apres file)
        """

        if fn is None:
            # Write these out so we can document them
            # Very basics
            self.snum = None  #: int number of samples per chirp
            self.data = None  #: np.ndarray(snum x tnum) of the actual return power
            self.data2 = None  #: np.ndarray(snum x tnum) of the actual return power
            self.dt = None  #: float, The spacing between samples in travel time, in seconds

            # Per-trace attributes
            #: np.ndarray(tnum,) of the acquisition time of each trace
            #: note that this is referenced to Jan 1, 0 CE (matlabe datenum)
            #: for convenience, use the `datetime` attribute to access a python version of the day
            self.decday = None
            self.decday2 = None
            #: np.ndarray(tnum,) latitude along the profile. Generally not in projected coordinates
            self.lat = None
            self.lat2 = None
            #: np.ndarray(tnum,) longitude along the profile. Generally not in projected coords.
            self.long = None
            self.long2 = None

            # Sample-wise attributes
            #: np.ndarray(snum,) The range for each sample, in m
            self.Rcoarse = None

            #: np.ndarray(tnum,) Optional. Projected x-coordinate along the profile.
            self.x_coord = None
            self.x_coord2 = None
            #: np.ndarray(tnum,) Optional. Projected y-coordinate along the profile.
            self.y_coord = None
            self.y_coord2 = None
            #: np.ndarray(tnum,) Optional. Elevation along the profile.
            self.elev = None
            self.elev2 = None

            # Special attributes
            #: impdar.lib.RadarFlags object containing information about the processing steps done.
            self.flags = ApresFlags()
            self.header = ApresHeader()

            self.data_dtype = None
            return

        if dat2 is None:
            if os.path.splitext(fn)[1] == '.h5':
                with h5py.File(fn, 'r') as fin:
                    grp = fin['dat']
                    for attr in grp.keys():
                        if attr in ['ApresFlags', 'ApresHeader']:
                            continue
                        val = grp[attr][:]
                        if isinstance(val, h5py.Empty):
                            val = None
                        setattr(self, attr, val)
                    for attr in grp.attrs.keys():
                        val = grp.attrs[attr]
                        if isinstance(val, h5py.Empty):
                            val = None
                        setattr(self, attr, val)
                    self.flags = ApresFlags()
                    self.header = ApresHeader()
                    self.flags.read_h5(grp)
                    self.header.read_h5(grp)
                # Set the file name
                self.fn = fn

            elif os.path.splitext(fn)[1] == '.mat':
                mat = loadmat(fn)
                for attr in self.attrs_guaranteed:
                    if mat[attr].shape == (1, 1):
                        setattr(self, attr, mat[attr][0][0])
                    elif mat[attr].shape[0] == 1 or mat[attr].shape[1] == 1:
                        setattr(self, attr, mat[attr].flatten())
                    else:
                        setattr(self, attr, mat[attr])
                # We may have some additional variables
                for attr in self.attrs_optional:
                    if attr in mat:
                        if mat[attr].shape == (1, 1):
                            setattr(self, attr, mat[attr][0][0])
                        elif mat[attr].shape[0] == 1 or mat[attr].shape[1] == 1:
                            setattr(self, attr, mat[attr].flatten())
                        else:
                            setattr(self, attr, mat[attr])
                    else:
                        setattr(self, attr, None)
                self.data_dtype = self.data.dtype
                self.flags = ApresFlags()
                self.flags.from_matlab(mat['flags'])
                self.header = ApresHeader()
                self.header.from_matlab(mat['header'])
                # Set the file name
                self.fn = fn

            else:
                raise ImportError('ApresDiff can only load .h5, .mat, or 2 data objects in parallel.')

        else:
            dat1 = fn
            if type(dat1) == str and type(dat2) == str:
                dat1 = ApresData(dat1)
                dat2 = ApresData(dat2)
            self.fn1 = dat1.header.fn
            self.fn2 = dat2.header.fn
            self.fn = self.fn1+'_diff_'+self.fn2

            if dat1.flags.range == 0 or dat2.flags.range == 0:
                raise TypeError('The range filter has not been executed on this data class, do that before proceeding.')
            if dat1.flags.stack == 0 or dat2.flags.stack == 0:
                raise TypeError('The stacking filter has not been executed on this data class, do that before proceeding.')

            # check that the two data objects are comparable
            if np.shape(dat1.data) != np.shape(dat2.data):
                raise TypeError('Acquisition inputs must be of the same shape.')
            if not np.all(abs(dat1.Rcoarse - dat2.Rcoarse) < 1e04):
                raise ValueError('Range vector should be the same for both acquisitions')

            # instantiate the differencing class
            self.data = dat1.data[0][0]
            self.data2 = dat2.data[0][0]
            self.Rcoarse = dat1.Rcoarse

            # if uncertainties were calculated, bring those in too
            if hasattr(dat1,'uncertainty'):
                self.unc1 = dat1.uncertainty
            if hasattr(dat2,'uncertainty'):
                self.unc2 = dat2.uncertainty

            # take the flags and header from the first data object
            self.flags = dat1.flags
            self.header = dat1.header

        self.check_attrs()

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
        idxs = np.arange(win//2,(len(self.data)-win//2),step).astype(int)
        if Rcoarse is not None:
            self.ds = Rcoarse[idxs]
        else:
            self.ds = self.Rcoarse[idxs]

        # Create data and coherence vectors
        acq1 = self.data
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
        self.w = phase2range(self,self.phi,
                self.header.lambdac,
                self.ds,
                self.header.chirp_grad,
                self.header.ci)

        if uncertainty == 'CR':
            # Error from Cramer-Rao bound, Jordan et al. (2020) Ann. Glac. eq. (5)
            sigma = (1./abs(self.co))*np.sqrt((1.-abs(self.co)**2.)/(2.*win))
            # convert the phase offset to a distance vector
            self.w_err = phase2range(self,sigma,
                    self.header.lambdac,
                    self.ds,
                    self.header.chirp_grad,
                    self.header.ci)

        elif uncertainty == 'noise_phasor':
            # Uncertainty from Noise Phasor as in Kingslake et al. (2014)
            # r_uncertainty should be calculated using the function phase_uncertainty defined in this script
            r_uncertainty = phase2range(self,self.unc1,self.header.lambdac) + phase2range(self,self.unc2,self.header.lambdac)
            idxs = np.arange(win//2,(len(self.data)-win//2),step)
            self.w_err = np.array([np.nanmean(r_uncertainty[i-win//2:i+win//2]) for i in idxs])

        if self.x_coord and self.x_coord2:
            # Calculate the surface velocity
            self.u = (self.x_coord2-self.x_coord,
                      self.y_coord2-self.y_coord)


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

        print('Calculating vertical strain rate over range from %s to %s meters.'%strain_window)
        idx = np.logical_and(self.ds>strain_window[0],self.ds<strain_window[1])
        slope, intercept, r_value, p_value, std_err = linregress(self.ds[idx],self.w[idx])
        self.eps_zz = slope
        self.w0 = intercept
        print('Vertical strain rate (yr-1):',self.eps_zz)
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

        P1 = 10.*np.log10(self.data**2.)
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


    def check_attrs(self):
        """Check if required attributes exist.

        This is largely for development only; loaders should generally call this method last,
        so that they can confirm that they have defined the necessary attributes.

        Raises
        ------
        ImpdarError
            If any required attribute is None or any optional attribute is fully absent"""
        for attr in self.attrs_guaranteed:
            if not hasattr(self, attr):
                raise ImpdarError('{:s} is missing. \
                    It appears that this is an ill-defined ApresDiff object'.format(attr))
            if getattr(self, attr) is None:
                raise ImpdarError('{:s} is None. \
                    It appears that this is an ill-defined ApresDiff object'.format(attr))

        if not hasattr(self, 'data_dtype') or self.data_dtype is None:
            self.data_dtype = self.data.dtype
        return
