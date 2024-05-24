#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""The Picks structure tracks picks and picking parameters."""
import numpy as np
from scipy.signal import filtfilt, butter

from .ImpdarError import ImpdarError
from .LastTrace import LastTrace
from .LeaderTrailer import LeaderTrailer
from .PickParameters import PickParameters


class Picks():
    """Information about picks.

    This object holds all picks for a given radargram. The main containers are
    matrices holding information about the indices in the data matrix a given
    pick lies.

    Attributes
    ----------
    samp1: nsamp x tnum array
        Min/max above the center of each pick. A new row for each new pick.
    samp2: nsamp x tnum array
        Max/min at the center of each pick. A new row for each new pick.
    samp3: nsamp x tnum array
        Max/min below the center of each pick. A new row for each new pick.
    time: nsamp x tnum array
        In StoDeep used to contain TWTT samp2. Since this is redundant I'm
        deptrecating it, and it will be zeros for impdar processed data.
    power: nsamp x tnum array
        Power across a pick. To get this in decibels, you need to take
        10. * np.log10(power)
    picknums: list of length nsamp
        The number of each pick.
    lasttrace: impdar.lib.LastTrace.LastTrace
        Information about the end of the trace
    lt: impdar.lib.LeaderTrailer.LeaderTrailer
        StoDeep legacy for compatibility, unused
    pickparams: impdar.lib.PickParameters.PickParameters
        This structure contains important information used in picking,
        such as frequency for picks.
    """

    attrs = ['samp1', 'samp2', 'samp3', 'time', 'power', 'picknums']
    flatten = [False, False, False, False, False, True]
    spec_attrs = ['lasttrace', 'lt', 'pickparams']

    def __str__(self):
        try:
            if self.samp1 is not None:
                approx_indices = np.nanmean(self.samp1, axis=1).astype(int)
                approx_indices[approx_indices < 0] = 0  # for NaNs
                mean_twtts = self.radardata.travel_time[approx_indices]
                if self.radardata.nmo_depth is not None:
                    assume_depth = ''
                    mean_depths = self.radardata.nmo_depth[approx_indices]
                else:
                    assume_depth = ' assuming 1.68e8 m/s vel'
                    mean_depths = mean_twtts / 2.0 * 1.68e3
                string = 'Pick object with {:d} picks:'.format(len(self.picknums))
                for i in range(len(self.picknums)):
                    if approx_indices[i] != 0:
                        string += '\n    pick {:d} at ~{:4.2f} us (~{:4.2f} m{:s})'.format(int(self.picknums[i]), mean_twtts[i], mean_depths[i], assume_depth)
                    else:
                        string += '\n    empty pick {:d}'.format(int(self.picknums[i]))
            else:
                string = 'Empty pick object'
        except (ValueError, TypeError, IndexError):
            string = 'Picks Object'
        return string

    def __init__(self, radardata, pick_struct=None):
        if pick_struct is not None:
            # Loading from a file
            for attr, flat in zip(self.attrs, self.flatten):
                setattr(self, attr, pick_struct[attr][0][0])
                if flat:
                    setattr(self, attr, getattr(self, attr).flatten())

                # Convert matlab zeros to Nones
                if getattr(self, attr).shape == (1, 1) and (
                        getattr(self, attr)[0][0] == 0):
                    setattr(self, attr, None)
            self.lasttrace = LastTrace(pick_struct['lasttrace'])
            self.lt = LeaderTrailer(radardata, pick_struct['lt'])
            self.pickparams = PickParameters(radardata,
                                             pick_struct['pickparams'])
            self.picknums = self.picknums.tolist()
        else:
            # Blank initialization
            self.samp1 = None
            self.samp2 = None
            self.samp3 = None
            self.time = None
            self.power = None
            self.picknums = None
            self.lasttrace = LastTrace()
            self.lt = LeaderTrailer(radardata)
            self.pickparams = PickParameters(radardata)

        self.radardata = radardata
        self.lines = []

    def add_pick(self, picknum=0):
        """Add a new pick.

        This method handles any complexity in adding a new pick. If no picks
        exist, it creates the matrices. Otherwise, it just adds rows. If the
        last pick hasn't been used, this just recycles that pick.

        Parameters
        ----------
        picknum: int, optional
            Number to call the new pick. Default zero

        Returns
        -------
        index: int
            The index of the row number containing the new pick

        Raises
        ------
        ValueError if the picknum already exists--we do not deal with repeats

        """
        if self.samp1 is None:
            # We have no matrices yet
            self.samp1 = np.zeros((1, self.radardata.tnum))
            self.samp2 = np.zeros((1, self.radardata.tnum))
            self.samp3 = np.zeros((1, self.radardata.tnum))
            self.time = np.zeros((1, self.radardata.tnum))
            self.power = np.zeros((1, self.radardata.tnum))
            self.samp1[0, :] = np.nan
            self.samp2[0, :] = np.nan
            self.samp3[0, :] = np.nan
            self.time[0, :] = np.nan
            self.power[0, :] = np.nan
            self.picknums = [picknum]
            self.lasttrace.add_pick(-9999, 0)
        elif np.all(np.isnan(self.samp1[-1, :])):
            # If the last pick is blank, we just overwrite it. Zero the pick.
            self.samp1[-1, :] = np.nan
            self.samp2[-1, :] = np.nan
            self.samp3[-1, :] = np.nan
            self.time[-1, :] = np.nan
            self.power[-1, :] = np.nan
            self.picknums[-1] = picknum
        else:
            # If loading from matlab, need to cast picknums as a list
            if type(self.picknums) == np.ndarray:
                self.picknums = self.picknums.flatten().tolist()
            if picknum in self.picknums:
                raise ValueError('We already have that pick')

            # We are just adding a row to the existing matrices of samples etc.
            self.samp1 = np.vstack((self.samp1,
                                    np.zeros((1, self.radardata.tnum))))
            self.samp2 = np.vstack((self.samp2,
                                    np.zeros((1, self.radardata.tnum))))
            self.samp3 = np.vstack((self.samp3,
                                    np.zeros((1, self.radardata.tnum))))
            self.time = np.vstack((self.time,
                                   np.zeros((1, self.radardata.tnum))))
            self.power = np.vstack((self.power,
                                    np.zeros((1, self.radardata.tnum))))
            self.samp1[-1, :] = np.nan
            self.samp2[-1, :] = np.nan
            self.samp3[-1, :] = np.nan
            self.time[-1, :] = np.nan
            self.power[-1, :] = np.nan
            self.lasttrace.add_pick(-9999, 0)

            self.picknums.append(picknum)
        return self.samp1.shape[0]

    def update_pick(self, picknum, pick_info):
        """Update a pick with new information.

        Rather than go in and manually update a pick every time it is changed,
        we take in all information about an individual pick simultaneously and
        this method updates the pick's information in the Pick object.

        Parameters
        ----------
        picknum: int
            The pick number to update. Must be in picknums
        pick_info: 5xtnum np.ndarray
            Array where rows are upper pick bound, pick center, lower pick
            bound, time (deprecated, usually left as zeros or nans), and power
            across the pick.

        Raises
        ------
        ValueError if picknum is not in picknums or if the shape of the
        pick_info is bad.

        """
        try:
            ind = self.picknums.index(picknum)
        except ValueError:
            raise ValueError('picknum provided is not a pick; you must you \
                             use a picknum not an index')

        if pick_info.shape != (5, self.radardata.tnum):
            raise ValueError('pick_info must be a 5xtnum array')

        self.samp1[ind, :] = pick_info[0, :]
        self.samp2[ind, :] = pick_info[1, :]
        self.samp3[ind, :] = pick_info[2, :]
        self.time[ind, :] = pick_info[3, :]
        self.power[ind, :] = pick_info[4, :]

    def smooth(self, lowpass, units='tnum'):
        """Smooth the picks.

        For now there are no choices on the filter--it is 3rd order Butterworth.
        Picks that have NaNs in the middle are left alone--this avoids edge effects.
        Power is not recalculated--too much risk of bias. Do manually at own risk.

        Parameters
        ----------
        lowpass: float
            The cutoff value for filtering, in units determined by 'units'
        units: str, optional
            The units in which lowpass are provided. Choices are tnum or dist, default tnum.

        Raises
        ------
        ValueError
            If the wavelength is less than 1 or greater than tnum.
            If the units are not in [dist, tnum].
        ImpDARError
            If units are dist but the data are not constant spaced.
            If the data have been elevation corrected.
        """
        # Quickly do nothing if we don't have picks
        if self.samp1 is None:
            return

        if (self.radardata.flags.interp is None or
                not self.radardata.flags.interp[0]) and units == 'dist':
            raise ImpdarError('Use units=tnum for non-respaced data')

        if self.radardata.flags.elev:
            raise ImpdarError('This will not work with elevation corrected data')

        tracespace = self.radardata.flags.interp[1]

        # Calculate the number of samples per wavelength.
        if units == 'dist':
            nsamp = lowpass / tracespace
        elif units == 'tnum':
            nsamp = lowpass
        else:
            raise ValueError('Units must be dist or tnum')
        if nsamp <= 2:
            raise ValueError('wavelength is too small, causing no samples per wavelength')
        if nsamp > self.radardata.tnum:
            raise ValueError('wavelength is too large, bigger than the whole radargram')

        high_corner_freq = 1. / float(nsamp)
        corner_freq = high_corner_freq * 2
        b, a = butter(3, corner_freq, 'low')
        padlen = 12

        for attr in ['samp1', 'samp2', 'samp3']:
            dat = getattr(self, attr)
            for row in range(dat.shape[0]):
                # We cannot smooth if there are gaps in the middle
                # But we do want to smooth everything, so iterate through non-nan chunks
                nn = np.where(~np.isnan(dat[row, :]))[0]
                isn = np.where(np.isnan(dat[row, :]))[0]
                if len(nn) == 0:
                    continue
                else:
                    start_ind = nn[0]
                while start_ind < self.radardata.tnum:
                    nans_remaining = isn[isn > start_ind]
                    if len(nans_remaining) > 0:
                        # Smooth a chunk that does not reach the right
                        end_ind = isn[isn > start_ind][0]

                        # need to check that the chunk is long enough to smooth
                        if end_ind - start_ind < padlen:
                            # If we still have remaining non-nans, update start_ind
                            # otherwise, we are done
                            if len(nn[nn > end_ind]) > 0:
                                start_ind = nn[nn > end_ind][0]
                                continue
                            else:
                                break

                        dat[row, start_ind:end_ind] = np.around(
                            filtfilt(b, a, dat[row, start_ind:end_ind], padlen=padlen))

                        # If we still have remaining non-nans, update start_ind
                        # otherwise, we are done
                        if len(nn[nn > end_ind]) > 0:
                            start_ind = nn[nn > end_ind][0]
                        else:
                            break
                    else:
                        # Everything left is not nan
                        if self.radardata.tnum - start_ind < nsamp:
                            break
                        dat[row, start_ind:] = np.around(filtfilt(b, a, dat[row, start_ind:], padlen=padlen))
                        break
            setattr(self, attr, dat)

    def reverse(self):
        """Flip left-right.

        Called by the overall RadarData.reverse
        """
        if self.samp1 is not None:
            self.samp1 = np.flip(self.samp1, 1)
        if self.samp2 is not None:
            self.samp2 = np.flip(self.samp2, 1)
        if self.samp3 is not None:
            self.samp3 = np.flip(self.samp3, 1)
        if self.power is not None:
            self.power = np.flip(self.power, 1)
        if self.time is not None:
            self.time = np.flip(self.time, 1)

    def hcrop(self, limits):
        """Crop to limits.

        Called by the overall RadarData.hcrop
        """
        attrs = ['samp1', 'samp2', 'samp3', 'time', 'power']
        for attr in attrs:
            val = getattr(self, attr)
            if val is not None:
                setattr(self, attr, val[:, limits[0]:limits[1]])

    def to_struct(self):
        """Convert to a format writable to a .mat file.

        Returns
        -------
        mat: dict
            Dictionary of attributes for export with scipy.io.savemat
        """
        mat = {}
        for attr in self.attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                mat[attr] = 0
        for attr in self.spec_attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr).to_struct()
            else:
                mat[attr] = 0
        return mat

    def crop(self, ind):
        """Crop the picks.

        This is just subtraction with some bounds checks.
        It is easy since we assume that the work to convert to indices has already been done.
        Not needed for bottom crops.

        Parameters
        ----------
        ind: int or ndarray(tnum, ):
            How much we need to shift by
        """
        for attr in ['samp1', 'samp2', 'samp3']:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                val = getattr(self, attr)

                # Be sure to preserve nans
                nanmask = np.isnan(val)
                val -= ind
                val[nanmask] = np.nan

                # And allow cropping such that picks are no longer within the window
                val[val < 0] = np.nan
                val[val >= self.radardata.snum] = np.nan

                setattr(self, attr, val)

    def restack(self, traces):
        for attr, nptype in zip(['samp1', 'samp2', 'samp3', 'time', 'power'], [int, int, int, float, float]):
            if hasattr(self, attr) and getattr(self, attr) is not None:
                val = getattr(self, attr)
                tnum = int(np.floor(val.shape[1] / traces))
                new_vals = np.zeros((val.shape[0], tnum))
                new_vals[:] = np.nan

                for j in range(tnum):
                    # It is not totally clear if this should be a mean or nanmean
                    new_vals[:, j] = np.nanmean(val[:, j * traces:min((j + 1) * traces, val.shape[1])], axis=1).astype(nptype)
                    new_vals[new_vals < 0] = np.nan
                    new_vals[new_vals >= self.radardata.snum] = np.nan
                setattr(self, attr, new_vals)
