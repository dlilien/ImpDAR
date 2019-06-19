import numpy as np
from scipy.signal import filtfilt, butter, tukey, cheby1, bessel, firwin, lfilter
from .migration_routines import *


class RadarDataFiltering:

    def adaptivehfilt(self, *args, **kwargs):
        """Adaptively filter to reduce noise in upper layers

        This subtracts the average of traces around an individual trace in order to filter it. You can call this method directly, or it can be called by sending the 'adaptive' option to :func:`RadarData.hfilt() <impdar.lib.RadarData.RadarData.hfilt>`"""
        # v3.1
        #	HFILTDEEP - This StoDeep subroutine processes bandpass filtered
        #		or NMO data to reduce the horizontal noise in the upper layers.
        #		The user need not specify any frequencies.  This program simply
        #		takes the average of all of the traces and subtracts it from the
        #		bandpassed data.  It will remove most horizontally-oriented
        #		features in a radar profile: ringing, horizontal reflectors.  It
        #		will also create artifacts at travel-times corresponding to any
        #		bright horizontal reflectors included in the average trace.
        #
        #       You will want to experiment with the creation of the average trace.
        #       Ideally choose an area where all reflectors are sloped and
        #       relatively dim so that they average out while the horizontal noise
        #       is amplified.  Note that generally there is no perfect horizontal
        #       filter that will work at all depths.  You will have to experiment
        #       to get the best results for your area of interest.
        #
        #	WARNING: Do not use hfiltdeep on elevation-corrected data!!!
        #
        #	Created: Logan Smith - 6/12/02
        #
        #   Modifications:
        #       1) Now has option to horizontally filter any data that exist in
        #       memory - L. Smith, 5/27/03
        #       2)Now uses menudeep instead of individual menu. K. Dwyer 6/3/03
        #       3)No longer main platform for horizontal filtering, acts as one
        #       option in a set of horizontal filters. K. Dwyer 6/3/03
        #       4)Coverted to function and added documentation. B. Welch 5/2/06
        #       5)Added double layer filtering and optimized for low gain data.
        #       Added smoothing of average trace--allows retention of more real
        #       data with an adaptive filter.  B. Youngblood 6/13/08
        #		6) Added flags structure to function input and output.  Also added
        #		code to set hfilt components of flags structure.  J. Olson 7/10/08
        #

        print('Adaptive filtering')
        #create average trace for first (rough) scan of data
        avg_trace = np.mean(self.data, axis=1)
        hfiltdata_mass = self.data - np.atleast_2d(avg_trace).transpose()

        #preallocate array
        avg_trace_scale = np.zeros_like(self.travel_time)

        # create a piecewise scaling function (insures that the filter only affects
        # the top layers of data)
        mask = self.travel_time <= 1.25
        avg_trace_scale[mask] = -0.1 * (self.travel_time[mask] - 0.25) * (self.travel_time[mask] - 0.25) + 1
        avg_trace_scale[~mask] = np.exp(-30. * ((self.travel_time[~mask] - 0.25) - 0.9) * ((self.travel_time[~mask] - 0.25) - 0.9))

        #preallocate array
        hfiltdata_scan_low = np.zeros_like(hfiltdata_mass, dtype=self.data.dtype)

        #begin looping through data
        for i in range(int(self.tnum)):
            # build a packet of 100 traces around the trace in question
            if i <= 50:
                scpacket = hfiltdata_mass[:, 0:100 - i]
            elif i >= self.tnum - 50:
                scpacket = hfiltdata_mass[:, int(self.tnum) - 100:int(self.tnum)]
            else:
                scpacket = hfiltdata_mass[:, i - 49:i + 50]

            # average the packet horizontally and double filter it (allows the
            # program to maintain small horizontal artifacts that are likely real)
            avg_trace_scan_low = filtfilt([.25, .25, .25, .25], 1, np.mean(scpacket, axis=-1)).flatten() * avg_trace_scale.flatten()

            # subtract the average trace off the data trace
            hfiltdata_scan_low[:, i] = hfiltdata_mass[:, i] - avg_trace_scan_low

        self.data = hfiltdata_scan_low.astype(self.data.dtype)
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
        """
        # v2.1
        #	HFILTDEEP - This StoDeep subroutine processes bandpass filtered
        #		or NMO data to reduce the horizontal noise in the upper layers.
        #		The user need not specify any frequencies.  This program simply
        #		takes the average of all of the traces and subtracts it from the
        #		bandpassed data.  It will remove most horizontally-oriented
        #		features in a radar profile: ringing, horizontal reflectors.  It
        #		will also create artifacts at travel-times corresponding to any
        #		bright horizontal reflectors included in the average trace.
        #
        #       You will want to experiment with the creation of the average trace.
        #       Ideally choose an area where all reflectors are sloped and
        #       relatively dim so that they average out while the horizontal noise
        #       is amplified.  Note that generally there is no perfect horizontal
        #       filter that will work at all depths.  You will have to experiment
        #       to get the best results for your area of interest.

        htr1 = int(max(0, min(ntr1, self.tnum - 1)))  # avoid a value less than 1 or greater than tnum
        htrn = int(max(htr1 + 1, min(ntr2, self.tnum)))
        print('Subtracting mean trace found between {:d} and {:d}'.format(htr1, htrn))

        avg_trace = np.mean(self.data[:, htr1:htrn], axis=-1)

        #taper average trace so it mostly affects only the upper layers in the
        #data
        avg_trace = avg_trace * (np.exp(-self.travel_time.flatten() * 0.05) / np.exp(-self.travel_time[0] * 0.05))
        self.data = self.data - np.atleast_2d(avg_trace).transpose().astype(self.data.dtype)
        print('Horizontal filter complete.')

        # set flags structure components
        self.flags.hfilt = np.ones((2,))

    def highpass(self, wavelength, tracespace):
        #HIGHPASSDEEP - This is NOT a highpass frequency filter--rather it is a
        # horizontal filter to be used after interpolation because our data now has
        # constant spacing and a constant time.  Note that this horizontal filter
        # requires constant trace-spacing in order to be effective.
        #
        #   You will want to experiment with the creation of the average trace.
        #   Ideally choose an area where all reflectors are sloped and relatively
        #   dim so that they average out while the horizontal noise is amplified.
        #   Note that generally there is no perfect horizontal filter that will
        #   work at all depths.  You will have to experiment to get the best
        #   results for your area of interest.
        #:func:`hfilt <impdar.lib.horizontal_filters.hfilt>`
        #	WARNING: Do not use highpassdeep on elevation-corrected data!!!
        #
        #
        # Created by L. Smith and modified by A. Hagen, 6/15/04.
        #
        #Modifications:
        #    1) Converted to function and added documentation. B. Welch 5/3/06
        #    2) Added ability to bypass user input when running batch files,
        #    increased function inputs and outputs. - J. Werner, 6/30/08.
        #	 3) Added flags structure to function input and output.  Also added
        #		code to set hfilt components of flags structure.  J. Olson 7/10/08
        #

        #Convert wavelength to meters.
        wavelength = int(wavelength * 1000)
        #Set an approximate sampling frequency (10ns ~ 10m --> 100MHz).
        fsamp = 100.
        #Calculate the number of samples per wavelength.
        nsamp = (wavelength - (wavelength % tracespace)) // tracespace
        if nsamp < 1:
            raise ValueError('wavelength is too small, causing no samples per wavelength')
        print('Sample resolution = {:d}'.format(nsamp))
        #The high corner frequency is the ratio of the sampling frequency (in MHz)
        #and the number of samples per wavelength (unitless).
        High_Corner_Freq = fsamp / float(nsamp)
        print('High cutoff at {:4.2f} MHz...'.format(High_Corner_Freq))

        Sample_Freq = 1. / self.dt

        Nyquist_Freq = Sample_Freq / 2.0
        #Convert High_Corner_Freq to Hz.
        High_Corner_Freq = High_Corner_Freq

        #Corner_Freq is used in the olaf_butter routine.
        Corner_Freq = High_Corner_Freq / Nyquist_Freq

        b, a = butter(5, Corner_Freq, 'high')

        self.data = filtfilt(b, a, self.data)

        # set flags structure components
        self.flags.hfilt = np.ones((2,))
        self.flags.hfilt[1] = 3

        print('Highpass filter complete.')

    def winavg_hfilt(self, avg_win, taper='full', filtdepth=100):
        #WINAVG_HFILTDEEP - This StoDeep subroutine performs a horizontal filter on
        #   the data to reduce the ringing in the upper layers. It uses a moving
        #   window to creat an average trace for each individual trace, applies an
        #   exponential taper to it and then subtracts it from the trace. This is
        #   based off of the original hfiltdeep and uses the moving window designed
        #   in cosine_win_hfiltdeep.  The extent of the taper can help to minimize
        #   artifacts that are often produced near regions of bright bed
        #   reflectors.
        #
        #   You will want to experiment with the creation of the average trace.
        #   Ideally choose an area where all reflectors are sloped and relatively
        #   dim so that they average out while the horizontal noise is amplified.
        #   Note that generally there is no perfect horizontal filter that will
        #   work at all depths.  You will have to experiment to get the best
        #   results for your area of interest.
        #
        #	WARNING: Do not use winavg_hfiltdeep on elevation-corrected data!!!
        #
        #Created: Kieran Dwyer 6/18/03
        #
        #Modifications:
        #    1) Coverted to function and added documentation. B. Welch 5/2/06
        #	 2) Added flags structure to function input and output.  Also added
        #		code to set hfilt components of flags structure.  J. Olson 7/10/08
        #

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

            #set taper to effect only the initial times
            exptaper[:filtdepth] = exptaper[:filtdepth] - exptaper[filtdepth]
            exptaper[filtdepth:self.snum] = 0
            exptaper = exptaper / np.max(exptaper)
        elif taper == 'tukey':  # pragma: no cover
            # This taper is weirdly hard-coded in StoDeep, and I'm assuming it is unused
            exptaper[1:30] = np.ones((30,))
            tukey_win = tukey(60, 0.5)
            exptaper[31:45] = tukey_win[46:60]
        else:
            raise ValueError('Unrecognized taper. Options are full, pexp, or tukey')

        hfiltdata = np.zeros_like(self.data, dtype=self.data.dtype)
        #set up ranges, create average, taper average, subtract average

        for i in range(int(self.tnum)):
            range_start = i - ((avg_win - 1) // 2)
            range_end = i + ((avg_win - 1) // 2)

            # As opposed to StoDeep, don't wrap just cutoff for simplicity
            if range_start < 0:
                range_start = 0
            if range_end > self.tnum:
                range_end = self.tnum
            #create an average trace
            avg_trace = np.mean(self.data[:, range_start:range_end], axis=-1) * exptaper

            #subtract avg_trace from each trace
            hfiltdata[:, i] = self.data[:, i] - avg_trace

        self.data = hfiltdata.copy()
        self.flags.hfilt = np.zeros((2,))
        self.flags.hfilt[1] = 2

        print('Horizontal filter complete.')

    def hfilt(self, ftype='hfilt', bounds=None):
        """Horizontally filter the data.

        This is a wrapper around other filter types. Horizontal filters are implemented (and documented) in the :mod:`impdar.lib.horizontal_filters` module.

        Parameters
        ----------
        ftype: str, optional
            The filter type. Options are :func:`hfilt <impdar.lib.horizontal_filters.hfilt>` and :func:`adaptive <impdar.lib.horizontal_filters.adaptivehfilt>`. Default hfilt
        bounds: tuple, optional
            Bounds for the hfilt. Default is None, but required if ftype is hfilt.
        """
        if ftype == 'hfilt':
            self.horizontalfilt(bounds[0], bounds[1])
        elif ftype == 'adaptive':
            self.adaptivehfilt()
        else:
            raise ValueError('Unrecognized filter type')

    def vertical_band_pass(self, low, high, order=5, filttype='butter', cheb_rp=5, fir_window='hamming', *args, **kwargs):
        """Vertically bandpass the data

        This function uses a forward-backward Butterworth filter to filter the data. Returns power that is not near the wavelength of the transmitter is assumed to be noise, so the limits for the filtering should generally surround the radar frequency. Some experimentation to see what provides the clearest results is probably needed for any given dataset.

        Parameters
        ----------
        low: float
            Lowest frequency passed, in MHz
        high: float
            Highest frequency passed, in MHz
        order: int
            Filter order (default 5)
        filttype: str, optional
            The filter type to use. Options are butter(worth), cheb(yshev type I), bessel, or FIR (finite impulse response). Default is butter.
        cheb_rp: float, optional
            Maximum ripple, in decibels, of Chebyshev filter. Only used if filttype=='cheb'. Default is 5.
        fir_window: str, optional
            The window type passed to scipy.signal.firwin. Only used if filttype=='fir'. Default is hamming'
        """
        # adapted from bandpass.m  v3.1 - this function performs a banspass filter in the
        # time-domain of the radar data to remove environmental noise

        # first determine the cut-off corner frequencies - expressed as a
        #	fraction of the Nyquist frequency (half the sample freq).
        # 	Note: all of this is in Hz

        Sample_Freq = 1.0 / self.dt  	# dt=time/sample (seconds)

        #calculate the Nyquist frequency
        Nyquist_Freq = 0.5 * Sample_Freq

        Low_Corner_Freq = low * 1.0e6
        High_Corner_Freq = high * 1.0e6

        corner_freq = np.zeros((2,))
        corner_freq[0] = Low_Corner_Freq / Nyquist_Freq
        corner_freq[1] = High_Corner_Freq / Nyquist_Freq

        # provide feedback to the user
        print('Bandpassing from {:4.1f} to {:4.1f} MHz...'.format(low, high))

        # FIR operates a little differently, and cheb has rp arg, so we need to do each case separately
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
            self.data[:-order, :] = lfilter(taps, 1.0, self.data, axis=0).astype(self.data.dtype)[order:, :]
        else:
            raise ValueError('Filter type {:s} is not recognized'.format(filttype))
        
        print('Bandpass filter complete.')

        # set flags structure components
        self.flags.bpass[0] = 1
        self.flags.bpass[1] = low
        self.flags.bpass[2] = high

    def migrate(self, mtype='stolt', **kwargs):
        """Migrate the data.

        This is a wrapper around all the migration routines in migration_routines.py.

        Parameters
        ----------
        mtype: str, optional
            The chosen migration routine. Options are: kirch, stolt, phsh.
            Default: stolt
        """
        if mtype == 'kirch':
            migrationKirchhoff(self, **kwargs)
        elif mtype == 'stolt':
            migrationStolt(self, **kwargs)
        elif mtype == 'phsh':
            migrationPhaseShift(self, **kwargs)
        elif mtype == 'su':
            migrationSeisUnix(self, **kwargs)
        else:
            raise ValueError('Unrecognized migration routine')

        # change migration flag
        self.flags.mig = mtype
