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

from .ApresFlags import ApresFlags
from .ApresHeader import ApresHeader
from ..ImpdarError import ImpdarError

class ApresData(object):
    """A class that holds the relevant information for an ApRES acquisition.

    We keep track of processing steps with the flags attribute.
    This base version's __init__ takes a filename of a .mat file in the old StODeep format to load.
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
    #: These may not have existed in old StoDeep files that we are compatible with,
    #: and they often cannot be set at the initial data load.
    #: If they exist, they all have units of meters.
    attrs_optional = ['lat',
                      'long',
                      'x_coord',
                      'y_coord',
                      'elev',
                      'temperature1',
                      'temperature2',
                      'battery_voltage',
                      'Rcoarse']

    from ._ApresDataProcessing import apres_range, phase_uncertainty, phase2range, range_diff, stacking
    from ._ApresDataSaving import save

    # Now make some load/save methods that will work with the matlab format
    def __init__(self, fn):
        if fn is None:
            # Write these out so we can document them
            # Very basics
            self.snum = None  #: int number of samples per chirp
            self.cnum = None  #: int, the number of chirps in a burst
            self.bnum = None  #: int, the number of bursts
            self.data = None  #: np.ndarray(snum x tnum) of the actual return power
            self.dt = None  #: float, The spacing between samples in travel time, in seconds

            # Per-trace attributes
            #: np.ndarray(tnum,) of the acquisition time of each trace
            #: note that this is referenced to Jan 1, 0 CE (matlabe datenum)
            #: for convenience, use the `datetime` attribute to access a python version of the day
            self.decday = None
            #: np.ndarray(tnum,) latitude along the profile. Generally not in projected coordinates
            self.lat = None
            #: np.ndarray(tnum,) longitude along the profile. Generally not in projected coords.
            self.long = None

            # chirp
            self.chirp_num = None  #: np.ndarray(cnum,) The 1-indexed number of the chirp
            self.chirp_att = None  #: np.ndarray(cnum,) Chirp attenuation settings
            self.chirp_time = None  #: np.ndarray(cnum,) Time at beginning of chirp (serial day)

            # Sample-wise attributes
            #: np.ndarray(snum,) The two way travel time to each sample, in us
            self.travel_time = None
            self.Rcoarse = None

            #: np.ndarray(tnum,) Optional. Projected x-coordinate along the profile.
            self.x_coord = None
            #: np.ndarray(tnum,) Optional. Projected y-coordinate along the profile.
            self.y_coord = None
            #: np.ndarray(tnum,) Optional. Elevation along the profile.
            self.elev = None

            # Special attributes
            #: impdar.lib.RadarFlags object containing information about the processing steps done.
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

        self.fn = fn
        self.header = ApresHeader()
        self.flags.from_matlab(mat['flags'])
        self.header.from_matlab(mat['header'])
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
                    It appears that this is an ill-defined RadarData object'.format(attr))
            if getattr(self, attr) is None:
                raise ImpdarError('{:s} is None. \
                    It appears that this is an ill-defined RadarData object'.format(attr))

        if not hasattr(self, 'data_dtype') or self.data_dtype is None:
            self.data_dtype = self.shh.dtype
        return

    @property
    def datetime(self):
        """A python operable version of the time of acquisition of each trace"""
        return np.array([datetime.datetime.fromordinal(int(dd)) + datetime.timedelta(days=dd % 1) - datetime.timedelta(days=366)
                         for dd in self.decday], dtype=np.datetime64)
