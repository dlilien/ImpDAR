#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
The basic class of ImpDAR. Methods for saving, processing, and filtering are defined externally.
"""


import datetime
import numpy as np
from scipy.io import loadmat
from ..RadarFlags import RadarFlags
from ..ImpdarError import ImpdarError
from ..Picks import Picks


class RadarData(object):
    """A class that holds the relevant information for a radar profile.

    We keep track of processing steps with the flags attribute.
    This base version's __init__ takes a filename of a .mat file in the old StODeep format to load.
    """
    #: Attributes that every RadarData object should have and should not be None.
    attrs_guaranteed = ['chan',
                        'data',
                        'decday',
                        'dt',
                        'lat',
                        'long',
                        'pressure',
                        'snum',
                        'tnum',
                        'trace_int',
                        'trace_num',
                        'travel_time',
                        'trig',
                        'trig_level']

    #: Optional attributes that may be None without affecting processing.
    #: These may not have existed in old StoDeep files that we are compatible with,
    #: and they often cannot be set at the initial data load.
    #: If they exist, they all have units of meters.
    attrs_optional = ['nmo_depth',
                      'elev',
                      'dist',
                      'x_coord',
                      'y_coord',
                      'fn']

    from ._RadarDataProcessing import reverse, nmo, crop, hcrop, restack, \
        rangegain, agc, constant_space, elev_correct, constant_sample_depth_spacing, \
        traveltime_to_depth
    from ._RadarDataSaving import save, save_as_segy, output_shp, output_csv, _get_pick_targ_info
    from ._RadarDataFiltering import adaptivehfilt, horizontalfilt, highpass, \
        winavg_hfilt, hfilt, vertical_band_pass, migrate

    # Now make some load/save methods that will work with the matlab format
    def __init__(self, fn_mat):
        if fn_mat is None:
            # Write these out so we can document them
            # Very basics
            self.snum = None  #: int number of samples per trace
            self.tnum = None  #: int, the number of traces in the file
            self.data = None  #: np.ndarray(snum x tnum) of the actual return power
            self.trace_int = None  #: float, the time between traces
            self.chan = None  #: The Channel number of the data
            self.dt = None  #: float, The spacing between samples in travel time, in seconds
            self.trig_level = None  #: float, The value on which the radar was triggering

            # Per-trace attributes
            #: np.ndarray(tnum,) of the acquisition time of each trace
            #: note that this is referenced to Jan 1, 0 CE (matlabe datenum)
            #: for convenience, use the `datetime` attribute to access a python version of the day
            self.decday = None
            #: np.ndarray(tnum,) latitude along the profile. Generally not in projected coordinates
            self.lat = None
            #: np.ndarray(tnum,) longitude along the profile. Generally not in projected coords.
            self.long = None
            #: np.ndarray(tnum,) of the distances along the profile.
            #: units will depend on whether geographic coordinate transforms,
            #: as well as GPS data, are available.
            self.dist = None
            #: np.ndarray(tnum,) The pressure at acquisition. ImpDAR does not use this at present.
            self.pressure = None
            self.trace_num = None  #: np.ndarray(tnum,) The 1-indexed number of the trace
            #: np.ndarray(tnum,) the index in each trace where the trigger was met.
            self.trig = None

            # Sample-wise attributes
            #: np.ndarray(snum,) The two way travel time to each sample, in us
            self.travel_time = None

            # Optional attributes
            #: str, the input filename. May be left as None.
            self.fn = None
            #: np.ndarray(tnum,) Optional. Projected x-coordinate along the profile.
            self.x_coord = None
            #: np.ndarray(tnum,) Optional. Projected y-coordinate along the profile.
            self.y_coord = None
            #: np.ndarray(tnum,) Optional. Elevation along the profile.
            self.elev = None
            #: np.ndarray(tnum,) Optional. Depth of each trace below the surface
            self.nmo_depth = None

            # Special attributes
            #: impdar.lib.RadarFlags object containing information about the processing steps done.
            self.flags = RadarFlags()
            #: impdar.lib.Picks object with information about picks in/interpretation of the data.
            #: The init method of picks needs some basic data to calculate frequencies, etc, so it
            #: is not created until it is needed (maybe after some modifications to the data).
            self.picks = None

            self.data_dtype = None
            return

        mat = loadmat(fn_mat)
        for attr in self.attrs_guaranteed:
            if attr not in mat:
                raise KeyError('.mat file does not appear to be in the StoDeep/ImpDAR format')
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

        self.fn = fn_mat
        self.flags = RadarFlags()
        self.flags.from_matlab(mat['flags'])
        if 'picks' not in mat:
            self.picks = Picks(self)
        else:
            self.picks = Picks(self, mat['picks'])
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

        for attr in self.attrs_optional:
            if not hasattr(self, attr):
                raise ImpdarError('{:s} is missing. \
                    It appears that this is an ill-defined RadarData object'.format(attr))

        if (self.data.shape != (self.snum, self.tnum)) and (self.elev is None):
            print(self.data.shape, (self.snum, self.tnum))
            raise ImpdarError('The data shape does not match the snum and tnum values!!!')

        if not hasattr(self, 'data_dtype') or self.data_dtype is None:
            self.data_dtype = self.data.dtype
        return

    @property
    def datetime(self):
        """A python operable version of the time of acquisition of each trace"""
        return np.array([datetime.datetime.fromordinal(int(dd)) + datetime.timedelta(days=dd % 1) - datetime.timedelta(days=366)
                         for dd in self.decday], dtype=np.datetime64)
