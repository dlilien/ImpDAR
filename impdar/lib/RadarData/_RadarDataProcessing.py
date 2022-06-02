#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU-GPL3.0 license.

"""
Define processing steps for Radar Data. These are all instance methods.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from ..permittivity_models import firn_permittivity
from ..ImpdarError import ImpdarError


def reverse(self):
    """Reverse radar data

    Essentially flip the profile left-right.
    This is desirable in a number of instances,
    but is particularly useful for concatenating profiles acquired in opposite directions.
    The St Olaf version of this function had a bunch of options.
    They seemed irrelevant to me. We just try to flip everything that might need flipping.
    """
    self.data = np.fliplr(self.data)

    self.x_coord = np.flip(self.x_coord, 0)
    self.y_coord = np.flip(self.y_coord, 0)
    self.decday = np.flip(self.decday, 0)
    self.lat = np.flip(self.lat, 0)
    self.long = np.flip(self.long, 0)
    self.elev = np.flip(self.elev, 0)
    if self.picks is not None:
        self.picks.reverse()

    # allow for re-reverasl
    if self.flags.reverse:
        print('Back to original direction')
        self.flags.reverse = False
    else:
        print('Profile direction reversed')
        self.flags.reverse = True


def constant_sample_depth_spacing(self):
    # First, we make a new array of depths
    if self.nmo_depth is None:
        raise AttributeError('Call nmo first...')
    if np.allclose(np.diff(self.nmo_depth), np.ones((self.snum - 1,)) * (self.nmo_depth[1] - self.nmo_depth[0])):
        print('No constant sampling when you already have constant sampling...')
        return 1

    depths = np.linspace(np.min(self.nmo_depth[0], 0), self.nmo_depth[-1], len(self.nmo_depth))
    self.data = interp1d(self.nmo_depth, self.data.transpose())(depths).transpose()
    self.travel_time = interp1d(self.nmo_depth, self.travel_time)(depths)
    self.nmo_depth = depths


def nmo(self, ant_sep, uice=1.69e8, uair=3.0e8, const_firn_offset=None, rho_profile=None, permittivity_model=firn_permittivity, const_sample=True):
    """Normal move-out correction.

    Converts travel time to distance accounting for antenna separation.
    This potentially affects data and snum.
    It also defines nmo_depth, the moveout-corrected depth

    Updated to accomodate vertical velocity profile. B Hills 8/2019

    Parameters
    ----------
    ant_sep: float
        Antenna separation in meters.
    uice: float or np.ndarray (2 x m), optional
        Speed of light in ice, in m/s. (different from StoDeep!!). Default 1.69e8.
    uair: float, optional
        Speed of light in air. Default 3.0e8
    const_firn_offset: float, optional
        Offset all depths by this much to account for firn-air.
        THIS IS THE TWO-WAY THICKNESS (not the firn-air).
        Useful if data start below thesurface so you can skip
        complicated, depth-dependent corrections. Default None.
    rho_profile: str,optional
        Filename for a csv file with density profile (depths in first column and densities in second)
        Units should be meters for depth, kgs per meter cubed for density.
        Note that using a variable uice will break the linear scaling between depth and time,
        so we are forced to choose whether the y-axis is linear in speed or time.
        I chose time, since this eases calculations like migration.
        For plotting vs. depth, however, the functions just use the bounds,
        so the depth variations are averaged out.
        You can use the helper function constant_sample_depth_spacing() in order to correct this,
        but you should call that after migration.
    permittivity_model: fun, optional
        density to permittivity model from the literature
    const_sample: bool, optional
        interpolate to constant sample spacing after the nmo correction
    """
    # Conversion of StO nmodeep v2.4
    #   Modifications:
    #       1) Adapted for Stodeep - L. Smith, 6/15/03
    #       2) Fixed rounding bug - L. Smith, 6/16/03
    #       3) Converted to function and added documentation - B. Welch, 5/8/06
    #       4) Rescaled uice at end of script to original input value for use
    #                   other programs.  - J. Olson, 7/3/08
    #       5) Increased function outputs - J. Werner, 7/7/08.
    #       6) Converted flag variables to structure, and updated function I/O
    #       - J. Werner, 7/10/08

    # need to crop out the pretrigger before doing the move-out
    if np.any(self.trig > 0):
        raise ImpdarError('Crop out the pretrigger before doing the nmo correction.')

    # --- Load the velocity profile --- #
    if rho_profile is not None:
        try:
            # load the density profile
            rho_profile_data = np.genfromtxt(rho_profile, delimiter=',')
            profile_depth = rho_profile_data[:, 0]
            profile_rho = rho_profile_data[:, 1]
        except IndexError:
            raise IndexError('Cannot load the depth-density profile')
        # Density to velocity
        eps = np.real(permittivity_model(profile_rho))
        profile_u = uair / np.sqrt(eps)
        # Interpolate velocity profile onto constant depth spacing
        d_interp = np.linspace(np.min(profile_depth, 0), max(profile_depth), self.snum)
        u_interp = interp1d(profile_depth, profile_u)(d_interp)

    # --- Do the move-out correction --- #

    # empty array to fill with new times
    nmotime = np.zeros((len(self.travel_time)))

    for i,t in enumerate(self.travel_time):
        if rho_profile is None:
            u_rms = uice
        else:
            # get RMS velocity used for correction
            d = minimize(optimize_moveout_depth,.5*t*uice,args=(t,ant_sep,d_interp,u_interp),
                         tol=1e-8,bounds=[(0,0.5*t*uair)])['x'][0]
            u_rms = np.sqrt(np.mean(u_interp[d_interp<d]**2.))
        # get the upper leg of the trave_path triangle (direct arrival) from the antenna separation and the rms velocity
        tsep_ice = 1e6*(ant_sep / u_rms)
        # hypotenuese, adjust to 'transmit time' by adding the separation time
        thyp = t+tsep_ice
        # calculate the vertical two-way travel time
        nmotime[i] = np.sqrt((thyp)**2. - tsep_ice**2.)

    # --- Cleanup --- #

    # save the updated time vector
    self.travel_time = nmotime
    # time to depth conversion
    if rho_profile is None:
        self.nmo_depth = self.travel_time / 2. * uice * 1.0e-6
    else:
        self.nmo_depth = traveltime_to_depth(self, profile_depth, profile_rho, c=uair, permittivity_model=permittivity_model)
    if const_sample:
        constant_sample_depth_spacing(self)

    if const_firn_offset is not None:
        self.nmo_depth = self.nmo_depth + const_firn_offset

    # Set flags
    try:
        self.flags.nmo[0] = 1
        self.flags.nmo[1] = ant_sep
    except (IndexError, TypeError):
        self.flags.nmo = np.ones((2, ))
        self.flags.nmo[1] = ant_sep


def optimize_moveout_depth(d_in, t, ant_sep, profile_depth, profile_u):
    """Optimize depth in the nmo filter.

    In the case of variable velocity, we need to iterate on the depth
    and rms velocity within this function until it converges.

    Parameters
    ----------
    d_in: float
        initial guess at the depth
    t: float
        time
    ant_sep: float
        antennae separation
    profile_depth: array
        depths corresponding to input velocity profile
    profile_u: array
        velocity
    """
    args = np.argwhere(profile_depth<=d_in)
    if len(args) == 0:
        raise ValueError('Profile not shallow enough. Extend to cover top')
    vels = profile_u[args][:,0]
    u_rms = np.sqrt(np.mean(vels**2.))
    return abs(np.sqrt((t/2.*u_rms)**2.-ant_sep**2.)-d_in)


def traveltime_to_depth(self, profile_depth, profile_rho, c=3.0e8, permittivity_model=firn_permittivity):
    """
    Convert travel_time to depth based on density profile

    This is called from within the nmo processing function
    It returns the depth for the moveout-corrected depth, in the case of a variable velocity

    Parameters
    ----------
    profile_depth: array
        input depths for measured densities
    profile_rho: array
        input densities to match depth locations above
    c: float, optional
        speed of light in vacuum
    permittivity_model: function
        specific density-to-permittivity model to use

    Returns
    -------
    depth: np.ndarray (self.snum x 1)
    """
    # get the input velocity-depth profile
    eps = np.real(permittivity_model(profile_rho))
    profile_u = c / np.sqrt(eps)
    # iterate over time, moving down according to the velocity at each step
    z = 0.
    depth = self.travel_time / 2. * c / np.sqrt(np.real(permittivity_model(917.))) * 1.0e-6
    for i, t in enumerate(self.travel_time):
        if t < 0.:
            continue
        elif t < self.dt * 1.0e6:
            step_u = profile_u[0]
            z += t / 2. * step_u * 1.0e-6
            depth[i] = z
        else:
            step_u = profile_u[np.nanargmin(abs(profile_depth - z))]
            z += self.dt / 2. * step_u
            depth[i] = z
    return depth


def crop(self, lim, top_or_bottom='top', dimension='snum', uice=1.69e8, rezero=True, zero_trig=True):
    """Crop the radar data in the vertical. We can take off the top or bottom.

    This will affect data, travel_time, and snum.

    Parameters
    ----------
    lim: float (int if dimension=='snum')
        The value at which to crop.
    top_or_bottom: str, optional
        Crop off the top (lim is the first remaining return) or the bottom
        (lim is the last remaining return).
    dimension: str, optional
        Evaluate in terms of sample (snum), travel time (twtt), or depth (depth).
        If depth, uses nmo_depth if present and use uice with no transmit/receive separation.
        If pretrig, uses the recorded trigger sample to crop.
    rezero: bool, optional
        Set the zero on the y axis to the cropped value (if cropping off the top). Default True.
        This is desirable if you are zeroing to the surface.
    zero_trig: bool, optional
        Reset the trigger to zero. Effectively asserts that the crop was to the surface. Default True.
    uice: float, optional
        Speed of light in ice. Used if nmo_depth is None and dimension=='depth'
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
        if not isinstance(self.trig, np.ndarray):
            ind = int(self.trig)
        else:
            ind = self.trig.astype(int)
    else:
        ind = int(lim)

    if not isinstance(ind, np.ndarray) or (dimension != 'pretrig'):
        if top_or_bottom == 'top':
            lims = [ind, self.data.shape[0]]
            self.trig = self.trig - ind
            if zero_trig:
                self.trig = np.zeros_like(self.trig)
        else:
            lims = [0, ind]
        self.data = self.data[lims[0]:lims[1], :]
        self.travel_time = self.travel_time[lims[0]:lims[1]]
        if rezero:
            self.travel_time = self.travel_time - self.travel_time[0]
        self.snum = self.data.shape[0]
        mintrig = 0
    else:
        # pretrig, vector input
        # Need to figure out if we need to do any shifting
        # The extra shift compared to the smallest
        mintrig = np.nanmin(ind)
        lims = [mintrig, self.data.shape[0]]
        self.trig = self.trig - ind
        data_old = self.data.copy()
        self.data = np.zeros((data_old.shape[0] - mintrig, data_old.shape[1]))
        self.data[:, :] = np.nan
        trig_ends = self.data.shape[0] - (ind - mintrig)
        for i in range(self.data.shape[1]):
            self.data[:trig_ends[i], i] = data_old[ind[i]:, i]
        self.travel_time = self.travel_time[lims[0]:lims[1]]
        if rezero:
            self.travel_time = self.travel_time - self.travel_time[0]
        self.snum = self.data.shape[0]

    if top_or_bottom == 'top':
        if self.picks is not None:
            self.picks.crop(ind)

    try:
        self.flags.crop[0] = 1
        self.flags.crop[2] = self.flags.crop[1] + lims[1]
    except (IndexError, TypeError):
        self.flags.crop = np.zeros((3,))
        self.flags.crop[0] = 1
        self.flags.crop[2] = self.flags.crop[1] + lims[1]
    self.flags.crop[1] = self.flags.crop[1] + lims[0]
    print('Vertical samples reduced to subset [{:d}:{:d}] of original'.format(
        int(self.flags.crop[1]), int(self.flags.crop[2])))


def hcrop(self, lim, left_or_right='left', dimension='tnum'):
    """Crop the radar data in the horizontal. We can take off the left or right.

    This will affect all trace-wise variables:
    data, lat, long, dist, x_coord, y_coord, elev, trace_num, trace_int, tnum,
    decday, pressure, and trig.

    Parameters
    ----------
    lim: float (int if dimension=='snum')
        The value at which to crop.
    top_or_bottom: str, optional
        Crop off the left (lim is the first remaining trace/dist) or the right
        (lim is the first trace deleted).
        Default is left
    dimension: str, optional
        Evaluate in terms of trace_num (tnum) or distance (dist).
        Note that trace_num is 1-indexed.
        This defines the units for left_or_right.
        Default is tnum.
    """
    if left_or_right not in ['left', 'right']:
        raise ValueError('left_or_right must be left or right, not {:s}'.format(left_or_right))
    if dimension not in ['tnum', 'dist']:
        raise ValueError('Dimension must be in ["tnum", "dist"]')

    if dimension == 'dist':
        if lim > np.max(self.dist):
            raise ValueError('lim is larger than largest distance')
        if lim <= 0:
            raise ValueError('Distance should be strictly positive')
        ind = np.min(np.argwhere(self.dist >= lim))
    else:
        if int(lim) in (0, 1):
            raise ValueError('lim should be at least two to preserve some data')
        if lim > self.tnum:
            raise ValueError('lim should be less than tnum+1 {:d} in order to do anything'.format(self.tnum + 1))
        if lim == -1 or lim < -int(self.tnum):
            raise ValueError('If negative, lim should be in [-self.tnum; -1)')
        ind = int(lim) - 1

    if left_or_right == 'left':
        lims = [ind, self.data.shape[1]]
    else:
        lims = [0, ind]

    # Most variable just need subsetting
    self.data = self.data[:, lims[0]:lims[1]]
    for var in ['lat', 'long', 'pressure', 'trace_int', 'trig', 'elev', 'x_coord', 'y_coord', 'decday']:
        # some of these are optional, and trig may be set to a float rather than an array
        if getattr(self, var) is not None and isinstance(getattr(self, var), np.ndarray):
            setattr(self, var, getattr(self, var)[lims[0]:lims[1]])

    if self.picks is not None:
        self.picks.hcrop(lims)

    # More complex modifications for these two
    if self.dist is not None:
        self.dist = self.dist[lims[0]:lims[1]] - self.dist[lims[0]]
    self.trace_num = self.trace_num[lims[0]:lims[1]] - lims[0] + 1

    # Finally tnum
    self.tnum = self.data.shape[1]


def restack(self, traces):
    """Restack all relevant data to the given number of traces.

    This function just takes the average of the given number of traces.
    This reduces file size and can get rid of noise.
    There are fancier ways to do this---
    if you have GPS, you probably want to restack to constant trace spacing instead.

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
    oned_restack_vars = ['dist',
                         'pressure',
                         'lat',
                         'long',
                         'x_coord',
                         'y_coord',
                         'elev',
                         'decday',
                         'trig']
    oned_newdata = {key: np.zeros((tnum, )) if getattr(self, key) is not None else None for key in oned_restack_vars}
    for j in range(tnum):
        stack[:, j] = np.mean(self.data[:, j * traces:min((j + 1) * traces, self.data.shape[1])],
                              axis=1)
        # trace_int[j] = np.sum(self.trace_int[j * traces:min((j + 1) * traces, self.data.shape[1])])
        for var, val in oned_newdata.items():
            if val is not None:
                val[j] = np.mean(getattr(self, var)[j * traces:
                                                    min((j + 1) * traces, self.data.shape[1])])
    self.tnum = tnum
    self.data = stack
    self.trace_num = np.arange(self.tnum).astype(int) + 1
    self.trace_int = trace_int

    if hasattr(self, 'picks') and self.picks is not None:
        self.picks.restack(traces)

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
    if isinstance(self.trig, (float, int, np.float, np.int64)):
        gain = self.travel_time[int(self.trig) + 1:] * slope
        self.data[int(self.trig + 1):, :] *= np.atleast_2d(gain).transpose()
    else:
        for i, trig in enumerate(self.trig):
            gain = self.travel_time[int(trig) + 1:] * slope
            self.data[int(trig) + 1:, i] *= gain
    self.flags.rgain = True


def agc(self, window=50, scaling_factor=50):
    """Try to do some automatic gain control

    This is from StoDeep--I'm not sure it is useful but it was easy to roll over so
    I'm going to keep it. I think you should have most of this gone with a bandpass,
    but whatever.

    Parameters
    ----------
    window: int, optional
        The size of window we use in number of samples (default 50)
    scaling_factor: int, optional
        The scaling factor. This gets divided by the max amplitude when we rescale the input.
        Default 50.
    """
    maxamp = np.zeros((self.snum,))
    # In the for loop, old code indexed used range(window // 2). This did not make sense to me.
    for i in range(self.snum):
        maxamp[i] = np.max(np.abs(self.data[max(0, i - window // 2):
                                            min(i + window // 2, self.snum), :]))
    maxamp[maxamp == 0] = 1.0e-6
    self.data *= (scaling_factor / np.atleast_2d(maxamp).transpose()).astype(self.data.dtype)
    self.flags.agc = True


def constant_space(self, spacing, min_movement=1.0e-2, show_nomove=False):
    """Restack the radar data to a constant spacing.

    This method uses the GPS information (i.e. the distance, x, y, lat, and lon),
    to do a 1-d interpolation to get new values in the radargram.
    It also updates related variables like lat, long, elevation, and coordinates.
    To avoid retaining sections of the radargram when the antenna was in fact stationary,
    some minimum movement between traces is enforced.
    This value is in meters, and should change to be commensurate with the collection strategy
    (e.g. skiing a radar is slower than towing it with a snowmobile).

    This function comprises the second half of what was done by StoDeep's interpdeep.
    If you have GPS data from an external, high-precision GPS, you would first want to call
    `impdar.lib.gpslib.kinematic_gps_control` so that the GPS-related variables are all
    improved, then you would want to call this method.
    `impdar.lib.gpslib` provides some wrappings for doing both steps
    and for loading in the external GPS data.

    Parameters
    ----------
    spacing: float
        Target trace spacing, in meters
    min_movement: float, optional
        Minimum trace spacing. If there is not this much separation, toss the next shot.
        Set high to keep everything. Default 1.0e-2.
    show_nomove: bool, optional
        If True, make a plot shading the areas where we think there is no movement.
        This can be really helpful for diagnosing what is wrong if you have lingering stationary traces.
        Untested.
    """
    # eliminate an interpolation error by masking out little movement
    good_vals = np.hstack((np.array([True]), np.diff(self.dist * 1000.) >= min_movement))

    # Correct the distances to reduce noise
    for i in range(len(self.dist)):
        if not good_vals[i]:
            self.dist[i:] = self.dist[i:] - (self.dist[i] - self.dist[i - 1])
    temp_dist = self.dist[good_vals]

    if show_nomove:  # pragma: no cover
        from ..plot import plot_radargram
        fig, ax = plot_radargram(self)
        ax.fill_between(self.trace_num, np.ones_like(self.dist) * np.max(self.travel_time), np.ones_like(self.dist) * np.min(self.travel_time), where=~good_vals, alpha=0.5)
        fig.show()

    new_dists = np.arange(np.min(temp_dist),
                          np.max(temp_dist),
                          step=spacing / 1000.0)

    # interp1d can only handle real values
    if self.data.dtype in [np.complex128]:
        self.data = interp1d(temp_dist, np.real(self.data[:, good_vals]))(new_dists) + 1.j * interp1d(temp_dist, np.imag(self.data[:, good_vals]))(new_dists)
    else:
        self.data = interp1d(temp_dist, self.data[:, good_vals])(new_dists)

    for attr in ['lat', 'long', 'x_coord', 'y_coord', 'decday', 'pressure', 'trig']:
        setattr(self,
                attr,
                interp1d(temp_dist, getattr(self, attr)[good_vals])(new_dists))
    for attr in ['elev']:
        if getattr(self, attr) is not None:
            setattr(self,
                    attr,
                    interp1d(temp_dist, getattr(self, attr)[good_vals])(new_dists))

    if self.picks is not None:
        for attr in ['samp1', 'samp2', 'samp3']:
            if getattr(self.picks, attr) is not None:
                setattr(self.picks, attr, np.round(interp1d(temp_dist, getattr(self.picks, attr)[:, good_vals])(new_dists)))
        for attr in ['power', 'time']:
            if getattr(self.picks, attr) is not None:
                setattr(self.picks, attr, interp1d(temp_dist, getattr(self.picks, attr)[:, good_vals])(new_dists))

    self.tnum = self.data.shape[1]
    self.trace_num = np.arange(self.tnum).astype(int) + 1
    self.dist = new_dists
    self.trace_int = np.hstack((np.array(np.nanmean(np.diff(self.dist))),
                                np.diff(self.dist))) * 1000.
    try:
        self.flags.interp[0] = 1
        self.flags.interp[1] = spacing
    except (IndexError, TypeError):
        self.flags.interp = np.ones((2,))
        self.flags.interp[1] = spacing


def elev_correct(self, v_avg=1.69e8):
    """Move the surface down in the data array to account for surface elevation.

    NMO depth attribute must have been created before elev_correct is called.
    This method should generally be called after you have interpolated precision GPS onto
    the data, otherwise the noise a handheld GPS will make the results look pretty bad.

    Be aware that this is not a precise conversion if your nmo correction had antenna
    separation, or if you used a depth-variable velocity structure. This is because we have
    a single vector describing depths in the data array, so only one value for each depth
    step.

    Parameters
    ----------
    v_avg: float, optional
        Average velocity. This is what will define the depth slices in the new data array.
        Default 1.69e8.

    Raises
    ------
    ValueError:
        If there is no nmo_depth since this is used for calculating depths

    """
    if self.nmo_depth is None:
        raise ValueError('Run nmo before elev_correct so that we have depth scale')
    # calculate number of rows that must be added
    elev_diffs = np.max(self.elev) - self.elev
    max_diff = np.max(elev_diffs)

    dz_avg = self.dt * (v_avg / 2.)
    max_samp = int(np.floor(max_diff / dz_avg))

    data_old = self.data.copy()
    self.data = np.zeros((data_old.shape[0] + max_samp, data_old.shape[1]))
    self.data[:, :] = np.nan

    top_inds = (elev_diffs / dz_avg).astype(int)
    for i in range(self.data.shape[1]):
        top_ind = top_inds[i]
        self.data[top_ind: top_ind + data_old.shape[0], i] = data_old[:, i]

    if hasattr(self, 'picks') and self.picks is not None:
        self.picks.crop(-top_inds)

    self.elevation = np.hstack((np.arange(np.max(self.elev), np.min(self.elev), -dz_avg),
                                np.min(self.elev) - self.nmo_depth))
    self.flags.elev = 1
