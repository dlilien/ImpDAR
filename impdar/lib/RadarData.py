#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Define a class that just has the necessary attributes for a StODeep file--this should be subclassed per filetype
"""

import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import butter, filtfilt
from scipy.interpolate import interp1d
from .horizontal_filters import hfilt, adaptivehfilt


class RadarData():
    """A class that holds the relevant information for a radar profile.
    
    We keep track of processing steps with the flags attribute. This thing gets subclassed per input filetype to override the init method, do any necessary initial processing, etc. This version's __init__ takes a filename of a .mat file in the old StODeep format to load.
    
    Attributes
    ----------
    data: :class:`numpy.ndarray`
        The data matrix.
    """
    chan = None
    data = None
    decday = None
    dist = None
    dt = None
    elev = None
    flags = None
    lat = None
    long = None
    pressure = None
    snum = None
    tnum = None
    trace_int = None
    trace_num = None
    travel_time = None
    trig = None
    trig_level = None
    x_coord = None
    y_coord = None
    nmo_depth = None
    picks = None
    attrs_guaranteed = ['chan', 'data', 'decday', 'dist', 'dt', 'elev', 'lat', 'long', 'pressure', 'snum', 'tnum', 'trace_int', 'trace_num', 'travel_time', 'trig', 'trig_level', 'x_coord', 'y_coord']
    attrs_optional = ['nmo_depth', 'picks']

    # Now make some load/save methods that will work with the matlab format
    def __init__(self, fn):
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
                self.attr = None
        self.flags = RadarFlags()
        self.flags.from_matlab(mat['flags'])

    def save(self, fn):
        """Save the radar data

        Parameters
        ----------
        fn: str
            Filename. Should have a .mat extension
        """
        mat = {}
        for attr in self.attrs_guaranteed:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                # this guards against error in matlab format
                mat[attr] = 0
        for attr in self.attrs_optional:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
        if self.flags is not None:
            mat['flags'] = self.flags
        else:
            # We want the structure available to prevent read errors from corrupt files
            mat['flags'] = RadarFlags()
        savemat(fn, mat)

    def vertical_band_pass(self, low, high, *args, **kwargs):
        """Vertically bandpass the data

        Parameters
        ----------
        low: float
            Lowest frequency passed, in MHz
        high: float
            Highest frequency passed, in MHz
        """
        # adapted from bandpass.m  v3.1 - this function performs a banspass filter in the
        # time-domain of the radar data to remove environmental noise.  The routine
        # currently uses a 5th order Butterworth filter.
        # We have a lot of power to mess with this because scipy. Keeping the butter for now.
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
        Sample_Freq = 1.0 / self.dt  	# dt=time/sample (seconds)

        #calculate the Nyquist frequency
        Nyquist_Freq = 0.5 * Sample_Freq

        Low_Corner_Freq = low * 1.0e6
        High_Corner_Freq = high * 1.0e6

        corner_freq = np.zeros((2,))
        corner_freq[0] = Low_Corner_Freq / Nyquist_Freq
        corner_freq[1] = High_Corner_Freq / Nyquist_Freq

        b, a = butter(5, corner_freq, 'bandpass')

        # provide feedback to the user
        print('Bandpassing from {:4.1f} to {:4.1f} MHz...'.format(low, high))
        self.data = filtfilt(b, a, self.data, axis=0).astype(self.data.dtype)
        print('Bandpass filter complete.')

        # set flags structure components
        self.flags.bpass[0] = 1
        self.flags.bpass[1] = low
        self.flags.bpass[2] = high

    def hfilt(self, ftype='hfilt', bounds=None):
        """Horizontally filter the data.

        This is a wrapper around other filter types.

        Parameters
        ----------
        ftype: str, optional
            The filter type. Options are 'hfilt' and 'adaptive'. Default hfilt
        bounds: tuple, optional
            Bounds for the hfilt. Default is None, but required for hfilt.
        """
        if ftype == 'hfilt':
            hfilt(self, bounds[0], bounds[1])
        elif ftype == 'adaptive':
            adaptivehfilt(self)
        else:
            raise ValueError('Unrecognized filter type')

    def reverse(self):
        """Reverse radar data

        The St Olaf version of this function had a bunch of options. They seemed irrelevant to me. We just try to flip everything that might need flipping.
        """
        self.data = np.fliplr(self.data)

        self.x_coord = np.fliplr(self.x_coord)
        self.y_coord = np.fliplr(self.y_coord)
        self.decday = np.fliplr(self.decday)
        self.lat = np.fliplr(self.lat)
        self.long = np.fliplr(self.long)
        self.elev = np.fliplr(self.elev)

        # allow for re-reverasl
        if self.flags.reverse:
            print('Back to original direction')
            self.flags.reverse = False
        else:
            print('Profile direction reversed')
            self.flags.reverse = True

    def nmo(self, ant_sep, uice=1.69e8, uair=3.0e8):
        """Normal move-out correction.

        Converts travel time to distance accounting for antenna separation. This potentially affects data and snum. It also defines nmo_depth, the moveout-corrected depth

        Parameters
        ----------
        ant_sep: float
            Antenna separation in meters.
        uice: float, optional
            Speed of light in ice, in m/s. (different from StoDeep!!). Default 1.69e8
        uar: float, optional
            Speed of light in air. Default 3.0e8
        """
        # Conversion of StO nmodeep v2.4
        #   Modifications:
        #       1) Adapted for Stodeep - L. Smith, 6/15/03
        #       2) Fixed rounding bug - L. Smith, 6/16/03
        #       3) Converted to function and added documentation - B. Welch, 5/8/06
        #		4) Rescaled uice at end of script to original input value for use
        #		in other programs.  - J. Olson, 7/3/08
        #       5) Increased function outputs - J. Werner, 7/7/08.
        #       6) Converted flag variables to structure, and updated function I/O
        #       - J. Werner, 7/10/08
        #

        #calculate travel time of air wave
        tair = ant_sep / uair

        if np.round(tair / self.dt) > self.trig:
            self.trig = int(np.round(1.1 * np.round(tair / self.dt)))
            nmodata = np.vstack((np.zeros((self.trig, self.data.shape[1])), self.data))
            self.snum = nmodata.shape[0]
        else:
            nmodata = self.data.copy()

        #calculate direct ice wave time
        tice = ant_sep / uice

        #calculate sample location for time=0 when radar pulse is transmitted
        #	(Note: trig is the air wave arrival)
        # switched from rounding whole right side to just rounding (tair/dt)
        #    -L. Smith, 6/16/03
        nair = int((self.trig + 1) - np.round(tair / self.dt))

        #calculate sample location for direct ice wave arrival
        # switched from rounding whole right side to just rounding (tice/dt)
        #    -L. Smith, 6/16/03
        nice = int(nair + np.round(tice / self.dt))

        #calculate vector of recorded travel times starting at time=0
        #calculate new time vector where t=0 is the emitted pulse
        pulse_time = np.arange(-self.dt * (nair), (self.snum - nair) * self.dt, self.dt)
        pulse_time[nice] = tice

        #create an empty vector
        tadj = np.zeros((self.snum, ))

        #calculate legs of travel path triangle (in one-way travel time)
        hyp = pulse_time[nice:self.snum] / 2.

        #calculate the vertical one-way travel time
        tadj[nice:self.snum] = (np.sqrt((hyp**2.) - ((tice / 2.)**2.)))
        tadj[np.imag(tadj) != 0] = 0

        #convert the vertical one-way travel time to vector indices
        nadj = np.zeros((self.snum, ))
        nadj[nice:self.snum] = pulse_time[nice: self.snum] / self.dt - tadj[nice: self.snum] * 2. / self.dt

        #loop through samples for all traces in profile
        for n in range(nice, self.snum):
            #correct the data by shifting the samples earlier in time by "nadj"
            nmodata[n - np.round(nadj[n]).astype(int), :] = nmodata[np.round(n).astype(int), :]
            #Since we are stretching the data near the ice surface we need to interpolate samples in the gaps. B. Welch June 2003
            if (n - np.round(nadj[n])) - ((n - 1) - np.round(nadj[n - 1])) > 1 and n != nice:
                interper = interp1d(np.array([((n - 1) - int(round(nadj[n - 1]))), (n - int(round(nadj[n])))]), nmodata[[((n - 1) - int(round(nadj[n - 1]))), (n - int(round(nadj[n])))], :].transpose())
                nmodata[((n - 1) - int(round(nadj[n - 1]))): (n - int(round(nadj[n]))), :] = interper(np.arange(((n - 1) - int(round(nadj[n - 1]))), (n - int(round(nadj[n]))))).transpose()

        #define the new pre-trigger value based on the nmo-adjustment calculations above:
        self.trig = nice - np.round(nadj[nice])

        #calculate the new variables for the y-axis after NMO is complete
        self.travel_time = np.arange((-self.trig) * self.dt, (nmodata.shape[0] - nair) * self.dt, self.dt) * 1.0e6
        self.nmo_depth = self.travel_time / 2. * uice * 1.0e-6
        
        self.data = nmodata

        self.flags.nmo[0] = 1
        try:
            self.flags.nmo[1] = ant_sep
        except IndexError:
            self.flags.nmo = np.ones((2, ))
            self.flags.nmo[1] = ant_sep

    def crop(self, lim, top_or_bottom='top', dimension='snum', uice=1.69e8):
        """Crop the radar data in the vertical. We can take off the top or bottom.

        This will affect data, travel_time, and snum.

        Parameters
        ----------
        lim: float (int if dimension=='snum')
            The value at which to crop.
        top_or_bottom: str, optional
            Crop off the top (lim is the first remaining return) or the bottom (lim is the last remaining return).
        dimension: str, optional
            Evaluate in terms of sample (snum), travel time (twtt), or depth (depth). If depth, uses nmo_depth if present and use uice with no transmit/receive separation.
        uice: float, optional
            Speed of light in ice. Used if nmo_depth is None and dimension=='depth'
        """
        if top_or_bottom not in ['top', 'bottom']:
            raise ValueError('top_or_bottom must be "top" or "bottom" not {:s}'.format(top_or_bottom))
        if dimension not in ['snum', 'twtt', 'depth']:
            raise ValueError('Dimension must be in [\'snum\', \'twtt\', \'depth\']')

        if dimension == 'twtt':
            ind = np.min(np.argwhere(self.travel_time >= lim))
        elif dimension == 'depth':
            if self.nmo_depth is not None:
                depth = self.nmo_depth
            else:
                depth = self.travel_time / 2. * uice * 1.0e-6
            ind = np.min(np.argwhere(depth >= lim))
        else:
            ind = int(lim)

        if top_or_bottom == 'top':
            lims = [ind, self.data.shape[0]]
        else:
            lims = [0, ind]

        self.data = self.data[lims[0]:lims[1], :]
        self.travel_time = self.travel_time[lims[0]:lims[1]] - self.travel_time[lims[0]]
        self.snum = self.data.shape[0]
        self.flags.crop[0] = 1
        try:
            self.flags.crop[2] = self.flags.crop[1] + lims[1]
        except IndexError:
            self.flags.crop = np.zeros((3,))
            self.flags.crop[0] = 1
            self.flags.crop[2] = self.flags.crop[1] + lims[1]
        self.flags.crop[1] = self.flags.crop[1] + lims[0]
        print('Vertical samples reduced to subset [{:d}:{:d}] of original'.format(int(self.flags.crop[1]), int(self.flags.crop[2])))

    def restack(self, traces):
        """Restack all relevant data to the given number of traces.

        There are fancier ways to do this--if you can, you probably want to restack to constant trace spacing instead

        Parameters
        ----------
        traces: int
            The (odd) number of traces to stack
        """
        traces = int(traces)
        if traces % 2 == 0:
            print('Only will stack odd numbers of traces. Using {:d}'.format(int(traces + 1)))
            traces = traces + 1
        tnum = int(np.floor(self.tnum / traces))
        stack = np.zeros((self.snum, tnum))
        trace_int = np.zeros((tnum, ))
        oned_restack_vars = ['dist', 'pressure', 'trig_level', 'lat', 'long', 'x_coord', 'y_coord', 'elev', 'decday']
        oned_newdata = {key: np.zeros((tnum, )) for key in oned_restack_vars}
        for j in range(tnum):
            stack[:, j] = np.mean(self.data[:, j * traces:min((j + 1) * traces, self.data.shape[1])], axis=1)
            trace_int[j] = np.sum(self.trace_int[j * traces:min((j + 1) * traces, self.data.shape[1])])
            for var, val in oned_newdata.items():
                val[j] = np.mean(getattr(self, var)[j * traces:min((j + 1) * traces, self.data.shape[1])])
        self.tnum = tnum
        self.data = stack
        self.trace_num = np.arange(self.tnum).astype(int) + 1
        self.trace_int = trace_int
        for var, val in oned_newdata.items():
            setattr(self, var, val)
        self.flags.restack = True

    def rangegain(self, slope):
        """Apply a range gain.

        Parameters
        ----------
        slope: float
            The slope of the linear range gain to be applied. Maybe try 1.0e-2?
        """
        gain = self.travel_time[self.trig + 1:] * slope
        self.data[self.trig + 1:, :] *= np.atleast_2d(gain).transpose()
        self.flags.rgain = True

    def agc(self, window=50, scaling_factor=50):
        """Try to do some automatic gain control

        This is from StoDeep--I'm not sure it is useful but it was easy to roll over so I'm going to keep it. I think you should have most of this gone with a bandpass, but whatever.

        Parameters
        ----------
        window: int, optional
            The size of window we use in number of samples (default 50)
        scaling_factor: int, optional
            The scaling factor. This gets divided by the max amplitude when we rescale the input. Default 50.
        """
        maxamp = np.zeros((self.snum,))
        for i in range(window // 2):
            maxamp[i] = np.max(np.abs(self.data[max(0, i - window // 2):min(i + window // 2, self.snum)]))
        maxamp[maxamp == 0] = 1.0e-6
        self.data *= (scaling_factor / np.atleast_2d(maxamp).transpose()).astype(self.data.dtype)
        self.flags.agc = True

    def constant_space(self, spacing, min_movement=1.0e-2):
        """Restack the radar data to a constant spacing.

        This comprises the second half of what was done by StoDeep's interpdeep. If you have good GPS data, you would first want to call `impdar.lib.gpslib.kinematic_gps_control` so that the GPS-related variables are all improved, then you would want to call this method. `impdar.lib.gpslib` provides some wrappings for doing both steps and for loading in the external GPS data

        Parameters
        ----------
        spacing: float
            Target trace spacing, in meters
        min_movement: float, optional
            Minimum trace spacing. If there is not this much separation, toss the next shot. Set high to keep everything. Default 1.0e-2.
        """
        # eliminate an interpolation error by masking out little movement
        good_vals = np.hstack((np.array([True]), np.diff(self.dist * 1000) >= min_movement))
        new_dists = np.arange(np.min(self.dist[good_vals]), np.max(self.dist[good_vals]), step=spacing / 1000.0)
        self.data = interp1d(self.dist[good_vals], self.data[:, good_vals])(new_dists)

        for attr in ['lat', 'long', 'elev', 'x_coord', 'y_coord', 'decday']:
            setattr(self, attr, interp1d(self.dist[good_vals], getattr(self, attr)[good_vals])(new_dists))

        self.tnum = self.data.shape[1]
        self.trace_num = np.arange(self.tnum).astype(int) + 1
        self.dist = new_dists
        try:
            self.flag.interp[0] = 1
            self.flag.interp[1] = spacing
        except IndexError:
            self.flag.interp = np.ones((2,))
            self.flag.interp[1] = spacing


class RadarFlags():
    """Flags that indicate the processing that has been used on the data.

    These are used for figuring out whether different processing steps have been performed. They also contain some information about the input arguments for some (but not all) of the processing steps.

    Attributes
    ----------
    batch: bool
        Legacy indication of whether we are batch processing. Always False.
    agc: bool
        Automatic gain control has been applied.
    reverse: bool
        Data have been reversed.
    restack: bool
        Data have been restacked.
    rgain: bool
        Data have a linear range gain applied.
    bpass: 3x1 :class:`numpy.ndarray`
        Elements: (1) 1 if bandpassed; (2) Low; and (3) High (MHz) bounds
    hfilt: 2x1 :class:`numpy.ndarray`
        Elements: (1) 1 if horizontally filtered; (2) Filter type
    interp: 2x1 :class:`numpy.ndarray`
        Elements: (1) 1 if constant distance spacing applied (2) The constant spacing (m)
    """

    def __init__(self):
        self.batch = False
        self.bpass = np.zeros((3,))
        self.hfilt = np.zeros((2,))
        self.rgain = False
        self.agc = False
        self.restack = False
        self.reverse = False
        self.crop = np.zeros((3,))
        self.nmo = np.zeros((2,))
        self.interp = np.zeros((2,))
        self.mig = False
        self.elev = 0
        self.attrs = ['batch', 'bpass', 'hfilt', 'rgain', 'agc', 'restack', 'reverse', 'crop', 'nmo', 'interp', 'mig', 'elev']
        self.bool_attrs = ['agc', 'batch', 'restack', 'reverse', 'rgain', 'mig']

    def to_matlab(self):
        """Convert all associated attributes into a dictionary formatted for use with :func:`scipy.io.savemat`
        """
        outmat = {att: getattr(self, att) for att in self.attrs}
        for attr in self.bool_attrs:
            outmat[attr] = 1 if outmat[attr] else 0
        return outmat

    def from_matlab(self, matlab_struct):
        """Associate all values from an incoming .mat file (i.e. a dictionary from :func:`scipy.io.loadmat`) with appropriate attributes
        """
        for attr in self.attrs:
            setattr(self, attr, matlab_struct[attr][0][0][0])
        for attr in self.bool_attrs:
            setattr(self, attr, True if matlab_struct[attr][0][0][0] == 1 else 0)
