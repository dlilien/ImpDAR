#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Benjamin Hills <bhills@uw.edu>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
An alternative ImpDAR class for ApRES data.
This should be considered separate from impulse data.
This class has a different set of loading and filtering scripts.

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 24 2019
"""

import os
import datetime

import numpy as np
from scipy.io import loadmat
import h5py

from .ApresFlags import ApresFlags, TimeDiffFlags, QuadPolFlags
from .ApresHeader import ApresHeader
from ..ImpdarError import ImpdarError

FILETYPE_OPTIONS = ['DAT', 'dat', 'mat', 'h5', 'nc']

class ApresData(object):
    """A class that holds the relevant information for an ApRES acquisition.

    Data are bundled across bursts into the same object if they are loaded together.
    Some information such as the instrument temperature, battery voltage, etc., are recorded for each burst.
    Some information is written to the header subclass instead of each burst individually.
    We track the processing steps with the flags attribute.
    """
    #: Attributes that every ApresData object should have and should not be None.
    attrs_guaranteed = ['data',
                        'decday',
                        'dt',
                        'snum',
                        'cnum',
                        'bnum',
                        'chirp_num',
                        'chirp_att',
                        'chirp_time',
                        'travel_time',
                        'frequencies']

    #: Optional attributes that may be None without affecting processing.
    attrs_optional = ['lat',
                      'long',
                      'x_coord',
                      'y_coord',
                      'elev',
                      'temperature1',
                      'temperature2',
                      'battery_voltage',
                      'Rcoarse',
                      'uncertainty']

    from ._ApresDataProcessing import apres_range, phase_uncertainty, stacking
    from ._ApresDataSaving import save

    # Now make some load/save methods that will work with the matlab format
    def __init__(self, fn):
        if fn is None:
            # Write these out so we can document them
            # Very basics
            self.snum = None  #: int number of samples per chirp
            self.cnum = None  #: int, the number of chirps in a burst
            self.bnum = None  #: int, the number of bursts
            self.data = None  #: np.ndarray(bnum, cnum, snum) of the actual return amplitude
            self.dt = None  #: float, The spacing between samples in travel time, in seconds
            self.uncertainty = None  #: float, phase uncertainty

            # Per-chirp attributes
            #: np.ndarray(bnum, cnum) of the acquisition time of each trace
            #: note that this is referenced to Jan 1, 0 CE (matlabe datenum)
            #: for convenience, use the `datetime` attribute to access a python version of the day
            self.decday = None
            #: np.ndarray(bnum, cnum) latitude/longitude along the profile. Generally not in projected coordinates
            self.lat = None
            self.long = None

            # chirp
            self.chirp_num = None  #: np.ndarray(bnum, cnum,) The 1-indexed number of the chirp
            self.chirp_att = None  #: np.ndarray(bnum, cnum,) Chirp attenuation settings
            self.chirp_time = None  #: np.ndarray(bnum, cnum,) Time at beginning of chirp (serial day)

            # Sample-wise attributes
            #: np.ndarray(snum,) The two way travel time to each sample, in us
            self.travel_time = None
            self.Rcoarse = None
            self.frequencies = None

            #: float Optional. Projected coordinates of the acquisition location
            self.x_coord = None
            self.y_coord = None
            self.elev = None

            #: float Optional. Temperatures
            self.temperature1 = None
            self.temperature2 = None
            self.battery_voltage = None

            # Special attributes
            self.flags = ApresFlags()
            self.header = ApresHeader()

            self.data_dtype = None
            return

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
        else:
            mat = loadmat(fn)
            for attr in self.attrs_guaranteed:
                if mat[attr].shape == (1, 1):
                    setattr(self, attr, mat[attr][0][0])
                elif (mat[attr].shape[0] == 1 or mat[attr].shape[1] == 1) and attr != 'data':
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

        self.fn = fn
        self.check_attrs()

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
                    It appears that this is an ill-defined ApresData object'.format(attr))
            if getattr(self, attr) is None:
                raise ImpdarError('{:s} is None. \
                    It appears that this is an ill-defined ApresData object'.format(attr))

        for attr in self.attrs_optional:
            if not hasattr(self, attr):
                raise ImpdarError('{:s} is missing. \
                    It appears that this is an ill-defined ApresData object'.format(attr))

        if not hasattr(self, 'data_dtype') or self.data_dtype is None:
            self.data_dtype = self.data.dtype
        return

    @property
    def datetime(self):
        """A python operable version of the time of acquisition of each trace"""
        return np.array([datetime.datetime.fromordinal(int(dd)) + datetime.timedelta(days=dd % 1) - datetime.timedelta(days=366)
                         for dd in self.decday], dtype=np.datetime64)


# ------------------------------------------------


class ApresTimeDiff(object):
    """
    Class for differencing between two Apres acquisitions.
    """
    #: Attributes that every ApresTimeDiff object should have and should not be None.
    attrs_guaranteed = ['data',
                        'data2',
                        'decday',
                        'decday2',
                        'dt',
                        'snum',
                        'range',
                        'fn1',
                        'fn2',
                        'fn']

    #: Optional attributes that may be None without affecting processing.
    #: These may not have existed in old StoDeep files that we are compatible with,
    #: and they often cannot be set at the initial data load.
    #: If they exist, they all have units of meters.
    attrs_optional = ['lat',
                      'lat2',
                      'long',
                      'long2',
                      'x_coord',
                      'x_coord2',
                      'y_coord',
                      'y_coord2',
                      'elev',
                      'elev2',
                      'unc1',
                      'unc2',
                      'ds',
                      'co',
                      'phi',
                      'w',
                      'w_err',
                      'w_0',
                      'eps_zz',
                      'bed']

    from ._TimeDiffProcessing import phase_diff, phase_unwrap, range_diff, strain_rate, bed_pick
    from ._ApresDataSaving import save

    def __init__(self,fn):
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
            self.data = None  #: np.ndarray(bnum, cnum, snum) of the actual return amplitude
            self.data2 = None  #: np.ndarray(bnum, cnum, snum) of the actual return amplitude
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
            self.range = None

            #: np.ndarray(tnum,) Optional. Projected x-coordinate along the profile.
            self.x_coord = None
            self.x_coord2 = None
            #: np.ndarray(tnum,) Optional. Projected y-coordinate along the profile.
            self.y_coord = None
            self.y_coord2 = None
            #: np.ndarray(tnum,) Optional. Elevation along the profile.
            self.elev = None
            self.elev2 = None

            # Differencing attributes
            self.ds = None
            self.co = None
            self.w = None

            # Special attributes
            #: impdar.lib.RadarFlags object containing information about the processing steps done.
            self.flags = TimeDiffFlags()
            self.header = ApresHeader()

            self.data_dtype = None
            return

        if os.path.splitext(fn)[1] == '.h5':
            with h5py.File(fn, 'r') as fin:
                grp = fin['dat']
                for attr in grp.keys():
                    if attr in ['TimeDiffFlags', 'ApresHeader']:
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
                self.flags = TimeDiffFlags()
                self.header = ApresHeader()
                self.flags.read_h5(grp)
                self.header.read_h5(grp)

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
            self.flags = TimeDiffFlags()
            self.flags.from_matlab(mat['flags'])
            self.header = ApresHeader()
            self.header.from_matlab(mat['header'])

        else:
            raise ImportError('ApresTimeDiff() is looking for an .h5 or .mat file \
                              saved as an Apdar object.')

        self.fn = fn
        self.check_attrs()

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
                    It appears that this is an ill-defined ApresTimeDiff object'.format(attr))
            if getattr(self, attr) is None:
                raise ImpdarError('{:s} is None. \
                    It appears that this is an ill-defined ApresTimeDiff object'.format(attr))

        if not hasattr(self, 'data_dtype') or self.data_dtype is None:
            self.data_dtype = self.data.dtype
        return


# ------------------------------------------------


class ApresQuadPol(object):
    """
    A class that holds the relevant information for a quad-polarized ApRES acquisition.
    """

    #: Attributes that every ApresQuadPol object should have and should not be None.
    attrs_guaranteed = ['data',
                        'shh',
                        'shv',
                        'svh',
                        'svv',
                        'range',
                        'decday',
                        'dt',
                        'snum',
                        'travel_time']

    #: Optional attributes that may be None without affecting processing.
    #: These may not have existed in old StoDeep files that we are compatible with,
    #: and they often cannot be set at the initial data load.
    #: If they exist, they all have units of meters.
    attrs_optional = ['lat',
                      'long',
                      'x_coord',
                      'y_coord',
                      'elev',
                      'ant_sep',
                      'ant_azi',
                      'thetas',
                      'HH',
                      'HV',
                      'VH',
                      'VV',
                      'chhvv',
                      'dphi_dz',
                      'cpe',
                      'cpe_idxs',
                      'chhvv_cpe',
                      'dphi_dz_cpe',
                      'phi']

    from ._QuadPolProcessing import rotational_transform, coherence2d, phase_gradient2d, find_cpe
    from ._ApresDataSaving import save

    # Now make some load/save methods that will work with the matlab format
    def __init__(self, fn):
        if fn is None:
            # Write these out so we can document them
            # Very basics
            self.data = None  #: np.ndarray(snum) of the actual return amplitude; copy of shh so that this class has a data array
            self.snum = None  #: int, number of samples per chirp
            self.dt = None  #: float, The spacing between samples in travel time, in seconds

            # Sample-wise attributes
            #: np.ndarray(snum,)
            self.shh = None  #: returned amplitude for hh polarization
            self.shv = None  #: returned amplitude for hv polarization
            self.svh = None  #: returned amplitude for vh polarization
            self.svv = None  #: returned amplitude for vv polarization
            self.travel_time = None  #: The two way travel time to each sample, in us

            # Float attributes relative to the time and location of the acquisition
            #: note that decimal days are referenced to Jan 1, 0 CE (matlabe datenum)
            #: for convenience, use the `datetime` attribute to access a python version of the day
            self.decday = None  #: acquisition time in days
            self.lat = None  #: latitude along the profile. Generally not in projected coordinates
            self.long = None  #: longitude along the profile. Generally not in projected coordinates
            self.x_coord = None  #: Optional. Projected x-coordinate along the profile.
            self.y_coord = None  #: Optional. Projected y-coordinate along the profile.
            self.elev = None  #: Optional. Elevation along the profile.

            # Special attributes
            #: impdar.lib.RadarFlags object containing information about the processing steps done.
            self.flags = QuadPolFlags()

            self.data_dtype = None
            return

        # Load from a file that has already been initialized in ImpDAR
        if os.path.splitext(fn)[1] == '.h5':
            with h5py.File(fn, 'r') as fin:
                grp = fin['dat']
                for attr in grp.keys():
                    if attr in ['QuadPolFlags']:
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
                self.flags = QuadPolFlags()
                self.flags.read_h5(grp)
        else:
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

            self.data_dtype = self.shh.dtype

            self.flags = QuadPolFlags()
            self.flags.from_matlab(mat['flags'])
            self.header = ApresHeader()
            self.header.from_matlab(mat['header'])

        self.fn = fn
        self.check_attrs()

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
                    It appears that this is an ill-defined ApresQuadPol object'.format(attr))
            if getattr(self, attr) is None:
                raise ImpdarError('{:s} is None. \
                    It appears that this is an ill-defined ApresQuadPol object'.format(attr))

        if not hasattr(self, 'data_dtype') or self.data_dtype is None:
            self.data_dtype = self.shh.dtype
        return

    @property
    def datetime(self):
        """A python operable version of the time of acquisition of each trace"""
        return np.array([datetime.datetime.fromordinal(int(dd)) + datetime.timedelta(days=dd % 1) - datetime.timedelta(days=366)
                         for dd in self.decday], dtype=np.datetime64)
