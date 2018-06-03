import numpy as np
from scipy.signal import filtfilt, butter, tukey


def adaptivehfilt(dat, *args, **kwargs):
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
    avg_trace = np.mean(dat.data, axis=1)
    hfiltdata_mass = dat.data - np.atleast_2d(avg_trace).transpose()

    #preallocate array
    avg_trace_scale = np.zeros_like(dat.travel_time)

    # create a piecewise scaling function (insures that the filter only affects
    # the top layers of data)
    mask = dat.travel_time <= 1.25
    avg_trace_scale[mask] = -0.1 * (dat.travel_time[mask] - 0.25) * (dat.travel_time[mask] - 0.25) + 1
    avg_trace_scale[~mask] = np.exp(-30. * ((dat.travel_time[~mask] - 0.25) - 0.9) * ((dat.travel_time[~mask] - 0.25) - 0.9))

    #preallocate array
    hfiltdata_scan_low = np.zeros_like(hfiltdata_mass, dtype=dat.data.dtype)

    #begin looping through data
    for i in range(int(dat.tnum)):
        # build a packet of 100 traces around the trace in question
        if i <= 50:
            scpacket = hfiltdata_mass[:, 0:100 - i]
        elif i >= dat.tnum - 50:
            scpacket = hfiltdata_mass[:, int(dat.tnum) - 100:int(dat.tnum)]
        else:
            scpacket = hfiltdata_mass[:, i - 49:i + 50]

        # average the packet horizontally and double filter it (allows the
        # program to maintain small horizontal artifacts that are likely real)
        avg_trace_scan_low = filtfilt([.25, .25, .25, .25], 1, np.mean(scpacket, axis=-1)).flatten() * avg_trace_scale.flatten()
        
        # subtract the average trace off the data trace
        hfiltdata_scan_low[:, i] = hfiltdata_mass[:, i] - avg_trace_scan_low
        
    dat.data = hfiltdata_scan_low.astype(dat.data.dtype)
    print('Adaptive filtering complete')

    # set flags structure components
    dat.flags.hfilt[0] = 1
    dat.flags.hfilt[1] = 4


def hfilt(dat, ntr1, ntr2, *args, **kwargs):
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

    htr1 = int(max(0, min(ntr1, dat.tnum - 1)))  # avoid a value less than 1 or greater than tnum
    htrn = int(max(htr1 + 1, min(ntr2, dat.tnum)))
    print('Subtracting mean trace found between {:d} and {:d}'.format(htr1, htrn))

    avg_trace = np.mean(dat.data[:, htr1:htrn], axis=-1)

    #taper average trace so it mostly affects only the upper layers in the
    #data
    avg_trace = avg_trace * (np.exp(-dat.travel_time.flatten() * 0.05) / np.exp(-dat.travel_time[0] * 0.05))
    dat.data = dat.data - np.atleast_2d(avg_trace).transpose().astype(dat.data.dtype)
    print('Horizontal filter complete.')

    # set flags structure components
    dat.flags.hfilt = np.ones((2,))


def highpass(dat, wavelength, tracespace):
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
    #
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
    wavelength = wavelength * 1000
    #Set an approximate sampling frequency (10ns ~ 10m --> 100MHz).
    fsamp = 100
    #Calculate the number of samples per wavelength.
    nsamp = (wavelength - (wavelength % tracespace)) // tracespace
    print('Sample resolution = {:d}'.format(nsamp))
    #The high corner frequency is the ratio of the sampling frequency (in MHz)
    #and the number of samples per wavelength (unitless).
    High_Corner_Freq = fsamp / nsamp
    print('	Lowpass cutoff at {:4.2f} MHz...'.format(High_Corner_Freq))

    Sample_Freq = 1 / dat.dt

    Nyquist_Freq = Sample_Freq / 2.0
    #Convert High_Corner_Freq to Hz.
    High_Corner_Freq = High_Corner_Freq * 1e6

    #Corner_Freq is used in the olaf_butter routine.
    Corner_Freq = High_Corner_Freq / Nyquist_Freq

    b, a = butter(5, Corner_Freq, 'high')

    dat.data = filtfilt(b, a, dat.data)

    # set flags structure components
    dat.flags.hfilt[1] = 1
    dat.flags.hfilt[2] = 3

    print('Highpass filter complete.')


def winavg_hfiltdeep(dat, avg_win, taper='full', filtdepth=100):
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

    if avg_win % 2 == 0:
        avg_win = avg_win + 1
        print('The averaging window must be an odd number of traces.')
        print('The averaging window has been changed to {:d}'.format(avg_win))

    if taper == 'full':
        exptaper = np.exp(-dat.travel_time * 0.05) / np.exp(-dat.travel_time[1] * 0.05)
    elif taper == 'pexp':
        filtdepth = filtdepth * 100 + dat.trig + 1
        exptaper = np.exp(-dat.travel_time * 0.05) / np.exp(-dat.travel_time[1] * 0.05)
        #set taper to effect only the initial times
        exptaper[1:filtdepth] = exptaper[1:filtdepth] - exptaper[filtdepth]
        exptaper[filtdepth:dat.snum] = 0
        exptaper = exptaper / np.max(exptaper)
    elif taper == 'tukey':
        exptaper = np.zeros((950,))
        exptaper[1:30] = np.ones((30,))
        tukey_win = tukey(60, 0.5)
        exptaper[31:45] = tukey_win[46:60]
    else:
        raise ValueError('Unrecognized taper. Options are full, pexp, or tukey')

    hfiltdata = np.zeros_like(dat.data)
    #set up ranges, create average, taper average, subtract average

    for i in range(int(dat.tnum)):
        range_start = i - ((avg_win - 1) / 2)
        range_end = i + ((avg_win - 1) / 2)

        # wrap range_start around to otherside of window to maintain
        # constant range size when i is less than 101 and range_start would
        # be zero or less
        if range_start < 0:
            range_end = range_end + range_start
            range_start = 0

        #wrap range_end around to beginning of window to maintain constant 
        #range size when range would be greater than tnum, this
        #perserves size of range.
        if range_end > dat.tnum:
            range_start = range_start - (range_end - dat.tnum)
            range_end = dat.tnum

        #create an average trace
        avg_trace = (np.sum(dat.data[:, range_start:range_end]) / (range_end - range_start)).flatten() * exptaper

        #subtract avg_trace from each trace
        hfiltdata[:, i] = dat.data[:, i] - avg_trace

    dat.data = hfiltdata.astype(dat.data.dtype)
    dat.flags.hfilt[1] = 1
    dat.flags.hfilt[2] = 2

    print('Horizontal filter complete.')
