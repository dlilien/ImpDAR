#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
A class that just has the necessary attributes for a StODeep file.

This is the basic object around which ImpDAR is written. It contains most of the methods used for processing and saving radar data. The source code is split across several files containing superclasses to keep file size down, but it all the methods should be documented here.

This base object is then subclassed for each type of data that ImpDAR can load, so that each of those subclasses can have easy initialization.
"""

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp1d
from ._RadarDataSaving import RadarDataSaving
from ._RadarDataFiltering import RadarDataFiltering
from .RadarFlags import RadarFlags
from .Picks import Picks

from .migration_routines import *
try:
    from osgeo import osr, ogr
    conversions_enabled = True
except ImportError:
    conversions_enabled = False


class RadarData(RadarDataSaving, RadarDataFiltering):
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
    fn = None
    attrs_guaranteed = ['chan', 'data', 'decday', 'dist', 'dt', 'elev', 'lat', 'long', 'pressure', 'snum', 'tnum', 'trace_int', 'trace_num', 'travel_time', 'trig', 'trig_level']
    attrs_optional = ['nmo_depth', 'elevation', 'x_coord', 'y_coord']

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
        self.fn = fn
        self.flags = RadarFlags()
        self.flags.from_matlab(mat['flags'])
        if ('picks' not in mat):
            self.picks = Picks(self)
        else:
            self.picks = Picks(self, mat['picks'])

    def reverse(self):
        """Reverse radar data

        Essentially flip the profile left-right. This is desirable in a number of instances, but is particularly useful for concatenating profiles that were acquired in different directions. The St Olaf version of this function had a bunch of options. They seemed irrelevant to me. We just try to flip everything that might need flipping.
        """
        self.data = np.fliplr(self.data)

        self.x_coord = np.flip(self.x_coord, 0)
        self.y_coord = np.flip(self.y_coord, 0)
        self.decday = np.flip(self.decday, 0)
        self.lat = np.flip(self.lat, 0)
        self.long = np.flip(self.long, 0)
        self.elev = np.flip(self.elev, 0)

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
        uair: float, optional
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
        nair = int((self.trig) - np.round(tair / self.dt))

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

            #Since we're stretching the data near the ice surface we need to interpolate samples in the gaps. B. Welch June 2003 I cannot hit this with a test. D. Lilien 4/2019
            if (n - np.round(nadj[n])) - ((n - 1) - np.round(nadj[n - 1])) > 1 and n != nice:
                interper = interp1d(np.array([((n - 1) - int(round(nadj[n - 1]))), (n - int(round(nadj[n])))]), nmodata[[((n - 1) - int(round(nadj[n - 1]))), (n - int(round(nadj[n])))], :].transpose())
                nmodata[((n - 1) - int(round(nadj[n - 1]))): (n - int(round(nadj[n]))), :] = interper(np.arange(((n - 1) - int(round(nadj[n - 1]))), (n - int(round(nadj[n]))))).transpose()

        #define the new pre-trigger value based on the nmo-adjustment calculations above:
        self.trig = nice - np.round(nadj[nice])

        #calculate the new variables for the y-axis after NMO is complete
        self.travel_time = np.arange((-self.trig) * self.dt, (nmodata.shape[0] - nair) * self.dt, self.dt) * 1.0e6
        self.nmo_depth = self.travel_time / 2. * uice * 1.0e-6

        self.data = nmodata

        try:
            self.flags.nmo[0] = 1
            self.flags.nmo[1] = ant_sep
        except (IndexError, TypeError):
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
            If pretrig, uses the recorded trigger sample to crop.
        uice: float, optional
            Speed of light in ice. Used if nmo_depth is None and dimension=='depth'

        Modified by Joshua Driscol, 03/04/2019
        Include an option to use the pretrigger value to clip
        """
        if top_or_bottom not in ['top', 'bottom']:
            raise ValueError('top_or_bottom must be "top" or "bottom" not {:s}'.format(top_or_bottom))
        if dimension not in ['snum', 'twtt', 'depth', 'pretrig']:
            raise ValueError('Dimension must be in [\'snum\', \'twtt\', \'depth\']')
        if top_or_bottom == 'bottom' and dimension == 'pretrig':
            raise ValueError('Only use pretrig to crop from the top')

        if dimension == 'twtt':
            ind = np.min(np.argwhere(self.travel_time >= lim))
        elif dimension == 'depth':
            if self.nmo_depth is not None:
                depth = self.nmo_depth
            else:
                depth = self.travel_time / 2. * uice * 1.0e-6
            ind = np.min(np.argwhere(depth >= lim))
        elif dimension == 'pretrig':
            if not (type(self.trig) == np.ndarray):
                ind = int(self.trig)
            else:
                ind = self.trig.astype(int)
        else:
            ind = int(lim)

        if not (type(ind) == np.ndarray) or (dimension != 'pretrig'):
            if top_or_bottom == 'top':
                lims = [ind, self.data.shape[0]]
            else:
                lims = [0, ind]
            self.data = self.data[lims[0]:lims[1], :]
            self.travel_time = self.travel_time[lims[0]:lims[1]] - self.travel_time[lims[0]]
            self.snum = self.data.shape[0]
        else:
            # pretrig, vector input
            # Need to figure out if we need to do any shifting
            # The extra shift compared to the smallest
            mintrig = np.min(ind)
            trig_ends = self.data.shape[0] - (ind - mintrig) - 1
            data_old = self.data.copy()
            self.data = np.zeros((data_old.shape[0] - mintrig, data_old.shape[1]))
            self.data[:, :] = np.nan
            for i in range(self.data.shape[1]):
                self.data[:trig_ends[i], i] = data_old[ind[i]:, i]
            lims = [0, mintrig]

        try:
            self.flags.crop[0] = 1
            self.flags.crop[2] = self.flags.crop[1] + lims[1]
        except (IndexError, TypeError):
            self.flags.crop = np.zeros((3,))
            self.flags.crop[0] = 1
            self.flags.crop[2] = self.flags.crop[1] + lims[1]
        self.flags.crop[1] = self.flags.crop[1] + lims[0]
        print('Vertical samples reduced to subset [{:d}:{:d}] of original'.format(int(self.flags.crop[1]), int(self.flags.crop[2])))

    def restack(self, traces):
        """Restack all relevant data to the given number of traces.

        This function just takes the average of the given number of traces. This reduces file size and can get rid of noise. There are fancier ways to do this--if you can, you probably want to restack to constant trace spacing instead.

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
        if type(self.trig) in [float, int]:
            gain = self.travel_time[int(self.trig) + 1:] * slope
            self.data[int(self.trig + 1):, :] *= np.atleast_2d(gain).transpose()
        else:
            for i, trig in enumerate(self.trig):
                gain = self.travel_time[int(trig) + 1:] * slope
                self.data[int(trig) + 1:, i] *= gain
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
        # In the for loop, old code indexed used range(window // 2). This did not make sense to me.
        for i in range(self.snum):
            print(i, max(0, i - window // 2), min(i + window // 2, self.snum))
            maxamp[i] = np.max(np.abs(self.data[max(0, i - window // 2):min(i + window // 2, self.snum), :]))
        maxamp[maxamp == 0] = 1.0e-6
        self.data *= (scaling_factor / np.atleast_2d(maxamp).transpose()).astype(self.data.dtype)
        self.flags.agc = True

    def constant_space(self, spacing, min_movement=1.0e-2):
        """Restack the radar data to a constant spacing.

        This method uses the GPS information (i.e. the distance, x, y, lat, and lon), to do a 1-d interpolation to get new values in the radargram. It also updates related variables like lat, long, elevation, and coordinates. To avoid retaining sections of the radargram when the antenna was in fact stationary, some minimum movement between traces is enforced. This value is in meters, and should change to be commensurate with the collection strategy (e.g. skiing a radar is slower than towing it with a snowmobile).

        This function comprises the second half of what was done by StoDeep's interpdeep. If you have GPS data from an external, high-precision GPS, you would first want to call `impdar.lib.gpslib.kinematic_gps_control` so that the GPS-related variables are all improved, then you would want to call this method. `impdar.lib.gpslib` provides some wrappings for doing both steps and for loading in the external GPS data.

        Parameters
        ----------
        spacing: float
            Target trace spacing, in meters
        min_movement: float, optional
            Minimum trace spacing. If there is not this much separation, toss the next shot. Set high to keep everything. Default 1.0e-2.
        """
        # eliminate an interpolation error by masking out little movement
        good_vals = np.hstack((np.array([True]), np.diff(self.dist * 1000.) >= min_movement))
        new_dists = np.arange(np.min(self.dist[good_vals]), np.max(self.dist[good_vals]), step=spacing / 1000.0)
        self.data = interp1d(self.dist[good_vals], self.data[:, good_vals])(new_dists)

        for attr in ['lat', 'long', 'elev', 'x_coord', 'y_coord', 'decday']:
            setattr(self, attr, interp1d(self.dist[good_vals], getattr(self, attr)[good_vals])(new_dists))

        self.tnum = self.data.shape[1]
        self.trace_num = np.arange(self.tnum).astype(int) + 1
        self.dist = new_dists
        try:
            self.flags.interp[0] = 1
            self.flags.interp[1] = spacing
        except (IndexError, TypeError):
            self.flags.interp = np.ones((2,))
            self.flags.interp[1] = spacing

    def elev_correct(self, v=1.69e8):
        if self.nmo_depth is None:
            raise ValueError('nmo must have been run before elev_correct so that we have depth scale')
        # calculate number of rows that must be added
        elev_diffs = np.max(self.elev) - self.elev
        max_diff = np.max(elev_diffs)

        dz = self.dt * (v / 2.)
        max_samp = int(np.floor(max_diff / dz))

        data_old = self.data.copy()
        self.data = np.zeros((data_old.shape[0] + max_samp, data_old.shape[1]))
        self.data[:, :] = np.nan

        for i in range(self.data.shape[1]):
            self.data[int(elev_diffs[i] // dz): int(elev_diffs[i] // dz) + data_old.shape[0], i] = data_old[:, i]

        self.elevation = np.hstack((np.arange(np.max(self.elev), np.min(self.elev), -dz), np.min(self.elev) - self.nmo_depth))
        self.flags.elev = 1
