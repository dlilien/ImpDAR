#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.
"""
The class methods for filtering.
"""
import numpy as np
from scipy.signal import filtfilt, butter, tukey, cheby1, bessel, firwin, lfilter, wiener
from scipy.ndimage import median_filter
from .. import migrationlib
from ..ImpdarError import ImpdarError


def adaptivehfilt(self, window_size, *args, **kwargs):
    """Adaptively filter to reduce noise in upper layers

    This subtracts the average of traces around an individual trace in order to filter it.
    You can call this method directly, or it can be called by sending the
    'adaptive' option to :func:`RadarData.hfilt() <impdar.lib.RadarData.RadarData.hfilt>`

    Parameters
    ----------
    window_size: int
        number of traces to include in the moving average to be removed

    Original StoDeep Documentation:
       HFILTDEEP-This StoDeep subroutine processes bandpass filtered
       or NMO data to reduce the horizontal noise in the upper layers.
       The user need not specify any frequencies.  This program simply
       takes the average of all of the traces and subtracts it from the
       bandpassed data.  It will remove most horizontally-oriented
       features in a radar profile: ringing, horizontal reflectors.  It
       will also create artifacts at travel-times corresponding to any
       bright horizontal reflectors included in the average trace.

       You will want to experiment with the creation of the average trace.
       Ideally choose an area where all reflectors are sloped and
       relatively dim so that they average out while the horizontal noise
       is amplified.  Note that generally there is no perfect horizontal
       filter that will work at all depths.  You will have to experiment
       to get the best results for your area of interest.


       WARNING: Do not use hfiltdeep on elevation-corrected data!!!

       Created: Logan Smith - 6/12/02
       Modified by: L. Smith, 5/27/03. K. Dwyer, 6/3/03.
       B. Welch, 5/2/06. B. Youngblood, 6/13/08.  J. Olson, 7/10/08
    """

    print('Adaptive filtering')

    # taper average trace so it mostly affects only the upper layers in the data
    avg_trace_scale = (np.exp(-self.travel_time.flatten() * 0.05) / np.exp(-self.travel_time[0] * 0.05))

    # preallocate array to fill with processed data
    ahfilt_data = np.zeros_like(self.data, dtype=self.data.dtype)

    # begin looping through data trace-by-trace
    for i in range(int(self.tnum)):
        # build a packet of window_size # of traces around the trace in question
        if i <= window_size // 2:
            scpacket = self.data[:, 0:window_size // 2 + i].copy()
        elif i >= self.tnum - window_size // 2:
            scpacket = self.data[:, int(self.tnum) - window_size:int(self.tnum)].copy()
        else:
            scpacket = self.data[:, i - window_size // 2 + 1:i + window_size // 2].copy()

        # average the packet horizontally and double filter it (allows the
        # program to maintain small horizontal artifacts that are likely real)
        avg_trace_scan_low = filtfilt([.25, .25, .25, .25],
                                      1,
                                      np.mean(scpacket, axis=-1)
                                      ).flatten() * avg_trace_scale.flatten()

        # subtract the average trace off the data trace
        # hfiltdata_scan_low[:, i] = hfiltdata_mass[:, i] - avg_trace_scan_low
        ahfilt_data[:, i] = self.data[:, i].copy() - avg_trace_scan_low

    self.data = ahfilt_data.astype(self.data.dtype)
    print('Adaptive filtering complete')

    # set flags structure components
    self.flags.hfilt[0] = 1
    self.flags.hfilt[1] = 4


def horizontalfilt(self, ntr1, ntr2, *args, **kwargs):
    """Remove the average trace.

    Parameters
    ----------
    ntr1: int
        Leftmost trace for averaging
    ntr2: int
        Rightmost trace for averaging


    Original StoDeep Documentation:
           HFILTDEEP - This StoDeep subroutine processes bandpass filtered
           or NMO data to reduce the horizontal noise in the upper layers.
           The user need not specify any frequencies.  This program simply
           takes the average of all of the traces and subtracts it from the
           bandpassed data.  It will remove most horizontally-oriented
           features in a radar profile: ringing, horizontal reflectors.  It
           will also create artifacts at travel-times corresponding to any
           bright horizontal reflectors included in the average trace.

           You will want to experiment with the creation of the average trace.
           Ideally choose an area where all reflectors are sloped and
           relatively dim so that they average out while the horizontal noise
           is amplified.  Note that generally there is no perfect horizontal
           filter that will work at all depths.  You will have to experiment
           to get the best results for your area of interest.
    """
    # avoid a value less than 1 or greater than tnum
    htr1 = int(max(0, min(ntr1, self.tnum - 1)))
    htrn = int(max(htr1 + 1, min(ntr2, self.tnum)))
    print('Subtracting mean trace found between {:d} and {:d}'.format(htr1, htrn))

    avg_trace = np.mean(self.data[:, htr1:htrn], axis=-1)

    # taper average trace so it mostly affects only the upper layers in the data
    avg_trace = avg_trace * (
        np.exp(-self.travel_time.flatten() * 0.05) / np.exp(-self.travel_time[0] * 0.05))
    self.data = self.data - np.atleast_2d(avg_trace).transpose().astype(self.data.dtype)
    print('Horizontal filter complete.')

    # set flags structure components
    self.flags.hfilt = np.ones((2,))


def highpass(self, wavelength):
    """High pass in the horizontal for a given wavelength.

    This only works if the data have constant trace spacing;
    we check the processing flags to enforce this.

    Parameters
    ----------
    wavelength: int
        The wavelength to pass, in meters.


    Original StoDeep Documentation:
        HIGHPASSDEEP - This is NOT a highpass frequency filter--rather it is a
        horizontal filter to be used after interpolation because our data now has
        constant spacing and a constant time.  Note that this horizontal filter
        requires constant trace-spacing in order to be effective.

        You will want to experiment with the creation of the average trace.
        Ideally choose an area where all reflectors are sloped and relatively
        dim so that they average out while the horizontal noise is amplified.
        Note that generally there is no perfect horizontal filter that will
        work at all depths.  You will have to experiment to get the best
        results for your area of interest.


        WARNING: Do not use highpassdeep on elevation-corrected data!!!


        Created by L. Smith and modified by
        A. Hagen, 6/15/04. B. Welch, 5/3/06. J. Werner, 6/30/08. J. Olson, 7/10/08
    """
    if self.flags.interp is None or not self.flags.interp[0]:
        raise ImpdarError('This method can only be used on constantly spaced data')

    if self.flags.elev:
        raise ImpdarError('This will not work with elevation corrected data')

    tracespace = self.flags.interp[1]

    # Convert wavelength to meters.
    wavelength = int(wavelength)
    # Set an approximate sampling frequency (10ns ~ 10m --> 100MHz).
    fsamp = 100.
    # Calculate the number of samples per wavelength.
    nsamp = int(wavelength / tracespace)
    if nsamp < 1:
        raise ValueError('wavelength is too small, causing no samples per wavelength')
    if nsamp > self.tnum:
        raise ValueError('wavelength is too large, bigger than the whole radargram')
    print('Sample resolution = {:d}'.format(nsamp))
    # The high corner frequency is the ratio of the sampling frequency (in MHz)
    # and the number of samples per wavelength (unitless).
    high_corner_freq = fsamp / float(nsamp)
    print('High cutoff at {:4.2f} MHz...'.format(high_corner_freq))

    sample_freq = 1. / self.dt

    nyquist_freq = sample_freq / 2.0
    # Convert High_Corner_Freq to Hz.
    high_corner_freq = high_corner_freq * 1.0e6

    # Corner_Freq is used in the olaf_butter routine.
    corner_freq = high_corner_freq / nyquist_freq

    b, a = butter(5, corner_freq, 'high')

    self.data = filtfilt(b, a, self.data)

    # set flags structure components
    self.flags.hfilt = np.ones((2,))
    self.flags.hfilt[1] = 3

    print('Highpass filter complete.')


def lowpass(self, wavelength):
    """Low pass in the horizontal for a given wavelength.

    This only works if the data have constant trace spacing;
    we check the processing flags to enforce this.

    Parameters
    ----------
    wavelength: int
        The wavelength to pass, in meters.


    Original StoDeep Documentation:
        HIGHPASSDEEP - This is NOT a highpass frequency filter--rather it is a
        horizontal filter to be used after interpolation because our data now has
        constant spacing and a constant time.  Note that this horizontal filter
        requires constant trace-spacing in order to be effective.

        You will want to experiment with the creation of the average trace.
        Ideally choose an area where all reflectors are sloped and relatively
        dim so that they average out while the horizontal noise is amplified.
        Note that generally there is no perfect horizontal filter that will
        work at all depths.  You will have to experiment to get the best
        results for your area of interest.


        WARNING: Do not use highpassdeep on elevation-corrected data!!!


        Created by L. Smith and modified by
        A. Hagen, 6/15/04. B. Welch, 5/3/06. J. Werner, 6/30/08. J. Olson, 7/10/08
    """
    if self.flags.interp is None or not self.flags.interp[0]:
        raise ImpdarError('This method can only be used on constantly spaced data')

    if self.flags.elev:
        raise ImpdarError('This will not work with elevation corrected data')

    tracespace = self.flags.interp[1]

    # Convert wavelength to meters.
    wavelength = int(wavelength)
    # Set an approximate sampling frequency (10ns ~ 10m --> 100MHz).
    fsamp = 100.
    # Calculate the number of samples per wavelength.
    nsamp = int(wavelength / tracespace)
    if nsamp < 1:
        raise ValueError('wavelength is too small, causing no samples per wavelength')
    if nsamp > self.tnum:
        raise ValueError('wavelength is too large, bigger than the whole radargram')
    print('Sample resolution = {:d}'.format(nsamp))
    # The high corner frequency is the ratio of the sampling frequency (in MHz)
    # and the number of samples per wavelength (unitless).
    high_corner_freq = fsamp / float(nsamp)
    print('Low cutoff at {:4.2f} MHz...'.format(high_corner_freq))

    sample_freq = 1. / self.dt

    nyquist_freq = sample_freq / 2.0
    # Convert High_Corner_Freq to Hz.
    high_corner_freq = high_corner_freq * 1.0e6

    # Corner_Freq is used in the olaf_butter routine.
    corner_freq = high_corner_freq / nyquist_freq

    b, a = butter(3, corner_freq, 'low')

    self.data = filtfilt(b, a, self.data)

    # set flags structure components
    self.flags.hfilt = np.ones((2,))
    self.flags.hfilt[1] = 3

    print('Lowpass filter complete.')


def horizontal_band_pass(self, low, high):
    """Bandpass in the horizontal for a given pair of wavelengths

    This only works if the data have constant trace spacing;
    we check the processing flags to enforce this.

    Parameters
    ----------
    low: float
        The minimum wavelength to pass, in meters.
    high: float
        The maximum wavelength to pass, in meters.

    """
    if self.flags.interp is None or not self.flags.interp[0]:
        raise ImpdarError('This method can only be used on constantly spaced data')
    if self.flags.elev:
        raise ImpdarError('This will not work with elevation corrected data')
    if low >= high:
        raise ValueError('Low must be less than high')

    tracespace = self.flags.interp[1]

    # Convert wavelength to meters.
    wavelength_high = high
    wavelength_low = low
    # Set an approximate sampling frequency (10ns ~ 10m --> 100MHz).
    fsamp = 100.
    # Calculate the number of samples per wavelength.
    nsamp_high = int(wavelength_high / tracespace)
    nsamp_low = int(wavelength_low / tracespace)
    if nsamp_low < 1:
        raise ValueError('Minimum wavelength is too small, causing no samples per wavelength')
    if nsamp_high > self.tnum:
        raise ValueError('Maximum wavelength is too long, causing more samples per wavelength than tnum, use lowpass instead?')
    print('Sample resolution high = {:d}'.format(nsamp_high))
    print('Sample resolution low = {:d}'.format(nsamp_low))
    # The high corner frequency is the ratio of the sampling frequency (in MHz)
    # and the number of samples per wavelength (unitless).
    high_corner_freq = fsamp / float(nsamp_high)
    low_corner_freq = fsamp / float(nsamp_low)

    nyquist_freq = fsamp / 2.0
    high_corner_freq = high_corner_freq
    low_corner_freq = low_corner_freq

    corner_freq = np.zeros((2,))
    corner_freq[0] = low_corner_freq / nyquist_freq
    corner_freq[1] = high_corner_freq / nyquist_freq

    b, a = butter(5, corner_freq, 'bandpass')

    self.data = filtfilt(b, a, self.data, axis=1)

    # set flags structure components
    self.flags.hfilt = np.ones((2,))
    self.flags.hfilt[1] = 3

    print('Highpass filter complete.')


def winavg_hfilt(self, avg_win, taper='full', filtdepth=100):
    """Uses a moving window to find the average trace, then subtracts this from the data.

    Parameters
    ----------
    avg_win: int
        The size of moving window. Must be odd and less than tnum.
        We will correct these rather than raise an exception.
    taper: str
        How the filter varies with depth. Options are full or pexp. For full,
        the power tapers exponentially. For pexp, the filter stops after the sample
        number given by filtdepth.


    Original StoDeep Documentation:
        WINAVG_HFILTDEEP - This StoDeep subroutine performs a horizontal filter on
        the data to reduce the ringing in the upper layers. It uses a moving
        window to create an average trace for each individual trace, applies an
        exponential taper to it and then subtracts it from the trace. This is
        based off of the original hfiltdeep and uses the moving window designed
        in cosine_win_hfiltdeep.  The extent of the taper can help to minimize
        artifacts that are often produced near regions of bright bed
        reflectors.

        You will want to experiment with the creation of the average trace.
        Ideally choose an area where all reflectors are sloped and relatively
        dim so that they average out while the horizontal noise is amplified.
        Note that generally there is no perfect horizontal filter that will
        work at all depths.  You will have to experiment to get the best
        results for your area of interest.

       WARNING: Do not use winavg_hfiltdeep on elevation-corrected data!!!

       Created: Kieran Dwyer 6/18/03
       Modified by B. Welch, 5/2/06. J. Olson, 7/10/08.
    """

    if avg_win > self.tnum:
        print('Cannot average over more than the whole data matrix. Reducing avg_win to tnum')
        avg_win = self.tnum
    if avg_win % 2 == 0:
        avg_win = avg_win + 1
        print('The averaging window must be an odd number of traces.')
        print('The averaging window has been changed to {:d}'.format(avg_win))
    exptaper = np.exp(-self.travel_time.flatten() * 0.05) / np.exp(-self.travel_time[0] * 0.05)
    if taper == 'full':
        # Don't modify the taper
        pass
    elif taper == 'pexp':
        # Commented out this next line since it seems wrong
        # filtdepth = filtdepth * 100 + self.trig + 1

        # set taper to effect only the initial times
        exptaper[:filtdepth] = exptaper[:filtdepth] - exptaper[filtdepth]
        exptaper[filtdepth:self.snum] = 0
        exptaper = exptaper / np.max(exptaper)
    elif taper == 'tukey':  # pragma: no cover
        # This taper is weirdly hard-coded in StoDeep, and I'm assuming it is unused
        exptaper[1:30] = np.ones((30,))
        tukey_win = tukey(60, 0.5)
        exptaper[31:45] = tukey_win[46:60]
        print('I am not sure this Tukey filter works for most data, use with caution')
    else:
        raise ValueError('Unrecognized taper. Options are full, pexp, or tukey')

    hfiltdata = np.zeros_like(self.data, dtype=self.data.dtype)
    # set up ranges, create average, taper average, subtract average

    for i in range(int(self.tnum)):
        range_start = i - ((avg_win - 1) // 2)
        range_end = i + ((avg_win - 1) // 2)

        # As opposed to StoDeep, don't wrap just cutoff for simplicity
        if range_start < 0:
            range_start = 0
        if range_end > self.tnum:
            range_end = self.tnum
        # Create an average trace
        avg_trace = np.mean(self.data[:, range_start:range_end], axis=-1) * exptaper

        # Subtract avg_trace from each trace
        hfiltdata[:, i] = self.data[:, i] - avg_trace

    self.data = hfiltdata.copy()
    self.flags.hfilt = np.zeros((2,))
    self.flags.hfilt[1] = 2

    print('Horizontal filter complete.')


def hfilt(self, ftype='hfilt', bounds=None, window_size=None):
    """Horizontally filter the data.

    This is a wrapper around other filter types.
    Horizontal filters are implemented (and documented) in the
    :mod:`impdar.lib.horizontal_filters` module.

    Parameters
    ----------
    ftype: str, optional
        The filter type. Options are :func:`hfilt <impdar.lib.horizontal_filters.hfilt>`
        and :func:`adaptive <impdar.lib.horizontal_filters.adaptivehfilt>`. Default hfilt
    bounds: tuple, optional
        Bounds for the hfilt. Default is None, but required if ftype is hfilt.
    window_size: int, optional
        number of traces in the moving average. Default is None, but required if ftype is adaptive.

    """
    if ftype == 'hfilt':
        self.horizontalfilt(bounds[0], bounds[1])
    elif ftype == 'adaptive':
        self.adaptivehfilt(window_size=window_size)
    else:
        raise ValueError('Unrecognized filter type')


def vertical_band_pass(self,
                       low,
                       high,
                       order=5,
                       filttype='butter',
                       cheb_rp=5,
                       fir_window='hamming',
                       *args,
                       **kwargs):
    """Vertically bandpass the data

    This function uses a forward-backward filter on the data.
    Returns power that is not near the wavelength of the transmitter is assumed to be noise,
    so the limits for the filtering should generally surround the radar frequency.
    Some experimentation may be needed to see what gives the clearest results for any given data

    There are a number of options for the filter type. Depending on the type of filter chosen,
    some other p

    Parameters
    ----------
    low: float
        Lowest frequency passed, in MHz
    high: float
        Highest frequency passed, in MHz
    order: int
        Filter order (default 5)
    filttype: str, optional
        The filter type to use. Options are butter(worth), cheb(yshev type I),
        bessel, or FIR (finite impulse response). Default is butter.
    cheb_rp: float, optional
        Maximum ripple, in decibels, of Chebyshev filter.
        Only used if filttype=='cheb'. Default is 5.
    fir_window: str, optional
        The window type passed to scipy.signal.firwin.
        Only used if filttype=='fir'. Default is hamming'
    """

    # first determine the cut-off corner frequencies - expressed as a
    # fraction of the Nyquist frequency (half the sample freq).
    # Note: all of this is in Hz
    sample_freq = 1.0 / self.dt  	# dt=time/sample (seconds)

    # Calculate the Nyquist frequency
    nyquist_freq = 0.5 * sample_freq

    low_corner_freq = low * 1.0e6
    high_corner_freq = high * 1.0e6

    corner_freq = np.zeros((2,))
    corner_freq[0] = low_corner_freq / nyquist_freq
    corner_freq[1] = high_corner_freq / nyquist_freq

    # provide feedback to the user
    print('Bandpassing from {:4.1f} to {:4.1f} MHz...'.format(low, high))

    # FIR operates a little differently, and cheb has rp arg,
    # so we need to do each case separately
    if filttype.lower() in ['butter', 'butterworth']:
        b, a = butter(order, corner_freq, 'bandpass')
        self.data = filtfilt(b, a, self.data, axis=0).astype(self.data.dtype)
    elif filttype.lower() in ['cheb', 'chebyshev']:
        b, a = cheby1(order, cheb_rp, corner_freq, 'bandpass')
        self.data = filtfilt(b, a, self.data, axis=0).astype(self.data.dtype)
    elif filttype.lower() == 'bessel':
        b, a = bessel(order, corner_freq, 'bandpass')
        self.data = filtfilt(b, a, self.data, axis=0).astype(self.data.dtype)
    elif filttype.lower() == 'fir':
        taps = firwin(order + 1, corner_freq, pass_zero=False)
        # I'm leaving the data past the filter--this is not filtfilt so we have a delay
        self.data[:-order, :] = lfilter(taps, 1.0, self.data, axis=0
                                        ).astype(self.data.dtype)[order:, :]
    else:
        raise ValueError('Filter type {:s} is not recognized'.format(filttype))

    print('Bandpass filter complete.')

    # set flags structure components
    self.flags.bpass[0] = 1
    self.flags.bpass[1] = low
    self.flags.bpass[2] = high


def denoise(self, vert_win=1, hor_win=10, noise=None, ftype='wiener'):
    """
    Denoising filter

    For now this just uses the scipy wiener or median filter,
    We could experiment with other options though

    Parameters
    ---------
    vert_win: int; optional
        vertical window size
    hor_win: int; optional
        horizontal window size
    noise: float; optional
        power of noise reduction for weiner, default is the average of the local variance
    ftype: string; optional
        filter type

    Raises
    ------

    """
    if ftype == 'wiener':
        if noise is None:
            # We want an error if there is no variance
            with np.errstate(divide='raise'):
                try:
                    self.data = wiener(self.data, mysize=(vert_win, hor_win))
                except FloatingPointError:
                    raise ValueError('Could not compute variance, specify noise for denoise')
        else:
            self.data = wiener(self.data, mysize=(vert_win, hor_win), noise=noise)
    elif ftype == 'median':
        self.data = median_filter(self.data, size=(vert_win, hor_win))
    else:
        raise ValueError('Only the wiener filter has been implemented for denoising.')


def migrate(self,
            mtype='stolt',
            vtaper=10,
            htaper=10,
            tmig=0,
            vel_fn=None,
            vel=1.68e8,
            nxpad=10,
            nearfield=False,
            verbose=0):
    """Migrate the data.

    This is a wrapper around all the migration routines in migration_routines.py.

    Parameters
    ----------
    mtype: str, optional
        The chosen migration routine. Options are: kirch, stolt, phsh.
        Default: stolt
    """
    if mtype == 'kirch':
        migrationlib.migrationKirchhoff(self, vel=vel, nearfield=nearfield)
    elif mtype == 'stolt':
        migrationlib.migrationStolt(self, vel=vel, htaper=htaper, vtaper=vtaper)
    elif mtype == 'phsh':
        migrationlib.migrationPhaseShift(self, vel=vel, vel_fn=vel_fn, htaper=htaper, vtaper=vtaper)
    elif mtype == 'tk':
        migrationlib.migrationTimeWavenumber(self,
                                             vel=vel,
                                             vel_fn=vel_fn,
                                             htaper=htaper,
                                             vtaper=vtaper)
    elif mtype[:2] == 'su':
        migrationlib.migrationSeisUnix(self,
                                       mtype=mtype,
                                       vel=vel,
                                       vel_fn=vel_fn,
                                       tmig=tmig,
                                       verbose=verbose,
                                       nxpad=nxpad,
                                       htaper=htaper,
                                       vtaper=vtaper)
    else:
        raise ValueError('Unrecognized migration routine')

    # change migration flag
    mflag = mtype
    self.flags.mig = mflag
