#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""

import numpy as np
from scipy.signal import butter, filtfilt


def vertical_band_pass(dat, low, high, *args, **kwargs):
    # function [filtdata,flags] = bandpassdeep(data_in,dt,tnum,snum,Low_Corner_Freq,High_Corner_Freq,flags)
    # bandpass.m  v2.1 - this function performs a banspass filter in the
    # time-domain of the radar data to remove environmental noise.  The routine
    # currently uses a 5th order Butterworth filter.
    #
    #Required inputs: 
    #   data_in - standard StoDeep radar data variable
    #   dt      - time per sample in radar data (in sec.)
    #   tnum    - number of traces in profile array
    #   snum    - number of samples in profile array
    #   Low_Corner_Freq - lower frequency limit for the bandpass filter (in MHz)
    #   High_Corner_Freq - upper frequency limit for the bandpass filter (in MHz)
    #	flags   - structure containing record of processes & parameters used on
    #			given file
    #
    #Outputs:
    #   filtdata - filtered radar data with dimensions (snum,tnum)
    #	flags   - structure containing record of processes & parameters used on
    #			given file
    #	
    #Created as stand alone script bandpass.m prior to 1997
    #  Modification history:
    #   1) Input changes made by A. Weitzel 7/10/97
    #   2) Switched to 5th-order Butterworth filter - P. Pearson, 7/2001
    #   3) Coverted for use in StoDeep and added pre-allocation of filtdata
    #       variable - B. Welch 10/2001
    #   4) Filters "stackdata" by default (if it exists),
    #		otherwise filters "data" - Peter Pearson, 2/13/02
    #   5) Now user can filter any standard StoDeep data variable that exists in
    #       memory using a menu displayed for the user - L. Smith, 5/27/03
    #   6) Converted to function and data variable is now passed to the function
    #       rather than selected within the script. - B. Welch, 5/1/06
    #	7) Updated input and outputs to include flags structure. Also added
    #		code to update flags structure - J. Olson 7/10/08
    #	

    # first determine the cut-off corner frequencies - expressed as a
    #	fraction of the Nyquist frequency (half the sample freq).
    # 	Note: all of this is in Hz
    Sample_Freq = 1.0 / dat.dt  	# dt=time/sample (seconds)

    #calculate the Nyquist frequency
    Nyquist_Freq = 0.5 * Sample_Freq

    Low_Corner_Freq = low * 1.0e6
    High_Corner_Freq = high * 1.0e6

    #now put the scaled value in the variable Corner_Freq
    corner_freq = np.zeros((2,))
    corner_freq[0] = Low_Corner_Freq / Nyquist_Freq
    corner_freq[1] = High_Corner_Freq / Nyquist_Freq

    b, a = butter(5, corner_freq, 'bandpass')

    # provide feedback to the user
    print('Bandpassing from {:4.1f} to {:4.1f} MHz...'.format(low, high))

    #pre-allocate the filtdata variable to make the script run faster
    filtdata = np.zeros((dat.snum, dat.tnum))

    #loop through the profile data trace by trace and filter the data
    for i in range(dat.tnum):
        #olaf_filtfilt.m is a St. Olaf version of the Matlab signal processing
        #version of the script.
        filtdata[:, i] = filtfilt(b, a, dat.data[:, i])
        #provide user feedback on the progress of the script
        if (dat.tnum - i) % 500 == 0:
            print('	Traces remaining: {:d}'.format(dat.tnum - i))

    print('Bandpass filter complete.')

    # places zeros in pre-trigger bins and re-generate time
    #BW 10/31/05 - commented out the following line that zeros out the
    #pre-trigger data.  We want to keep the pre-trigger data as a measure
    #of background noise levels.
    #filtdata(1:(trig-1),1:sz(2)) = zeros((trig-1),sz(2));

    # set flags structure components
    dat.flags.bpass = np.zeros((3,))
    dat.flags.bpass[0] = 1
    dat.flags.bpass[1] = Low_Corner_Freq / 1e6
    dat.flags.bpass[2] = High_Corner_Freq / 1e6


def horizontal_band_pass(*args, **kwargs):
    pass
