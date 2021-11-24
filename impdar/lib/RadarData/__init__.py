#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""The basic class of ImpDAR.

Methods for saving, processing, and filtering are defined externally.
"""


import datetime
import numpy as np
from scipy.io import loadmat
from ..RadarFlags import RadarFlags
from ..ImpdarError import ImpdarError
from ..Picks import Picks
from .. import gpslib

STODEEP_ATTRS = ['data', 'migdata', 'interp_data', 'nmo_data', 'filtdata', 'hfilt_data']


class RadarData(object):
    """A class that holds the relevant information for a radar profile.

    We keep track of processing steps with the flags attribute.
    This base version's __init__ takes a filename of a .mat file in the old
    StODeep format to load.
    """

    #: Attributes that every RadarData object should have.
    #: These should not be None.
    attrs_guaranteed = ['chan',
                        'data',
                        'decday',
                        'dt',
                        'pressure',
                        'snum',
                        'tnum',
                        'trace_int',
                        'trace_num',
                        'travel_time',
                        'trig',
                        'trig_level']

    #: Optional attributes that may be None without affecting processing.
    #: These may not have existed in old StoDeep files,
    #: and they often cannot be set at the initial data load.
    #: If they exist, they all have units of meters.
    attrs_optional = ['nmo_depth',
                      'lat',
                      'long',
                      'elev',
                      'dist',
                      'x_coord',
                      'y_coord',
                      'fn',
                      't_srs']

    #: The names of the optional StoDeep data matrices
    stodeep_attrs = STODEEP_ATTRS

    def __str__(self):
        try:
            if (self.snum is not None) and (self.tnum is not None):
                string = '{:d}x{:d} RadarData object'.format(self.snum, self.tnum)
                proc = False
                if (self.flags.bpass is not None) and (self.flags.bpass[0]):
                    proc = True
                    string += ', vertically bandpassed {:4.1f}:{:4.1f} Mhz'.format(self.flags.bpass[0], self.flags.bpass[1])
                if (self.flags.hfilt is not None) and (self.flags.hfilt[0]):
                    proc = True
                    string += ', horizontally filtered'
                if (self.flags.interp is not None) and (self.flags.interp[0]):
                    proc = True
                    string += ', re-interpolated to {:4.2f}-m spacing'.format(self.flags.interp[1])
                if (self.flags.crop is not None) and (self.flags.crop[0]):
                    proc = True
                    string += ', cropped to {:d}:{:d}'.format(int(self.flags.crop[1]), int(self.flags.crop[2]))
                if self.nmo_depth is not None:
                    string += ', moveout-corrected'
                if (self.flags.restack is not None) and self.flags.restack > 0:
                    proc = True
                    string += ', restacked by {:d}'.format(int(self.flags.restack))
                if (self.flags.mig is not None) and (self.flags.mig != 'none'):
                    proc = True
                    string += ', migrated'
                if not proc:
                    string += ', unprocessed'
                string += '.\n'

                if self.fn is not None:
                    string += '\n    from file {:s}'.format(self.fn)
                if self.x_coord is not None:
                    string += '\n    Projected geographic coordinates'
                    if self.t_srs is not None:
                        string += (': ' + self.t_srs)
                elif self.lat is not None:
                    string += '\n    Unprojected geographic coordinates'
                if (self.picks is not None) and (self.picks.samp1 is not None):
                    string += ('\nAssociate picks are: ' + str(self.picks))
                else:
                    string += '\nno picks'
            else:
                string = 'RadarData object, undefined dimensions'
        except (ValueError, TypeError, IndexError):
            string = 'RadarData Object'
        return string

    from ._RadarDataProcessing import reverse, nmo, crop, hcrop, restack, \
        rangegain, agc, constant_space, elev_correct, \
        constant_sample_depth_spacing, traveltime_to_depth
    from ._RadarDataSaving import save, save_as_segy, output_shp, output_csv, \
        _get_pick_targ_info
    from ._RadarDataFiltering import adaptivehfilt, horizontalfilt, highpass, \
        winavg_hfilt, hfilt, vertical_band_pass, denoise, migrate, \
        horizontal_band_pass, lowpass

    # Now make some load/save methods that will work with the matlab format
    def __init__(self, fn_mat):
        if fn_mat is None:
            # Store this for possible later filename modification
            self.fn = fn_mat

            # Write these out so we can document them
            # Very basics

            #: int number of samples per trace
            self.snum = None
            #: int, the number of traces in the file
            self.tnum = None
            #: np.ndarray(snum x tnum) of the actual return power.
            self.data = None
            #: float, the time between traces.
            self.trace_int = None
            #: The Channel number of the data.
            self.chan = None
            #: float, The spacing between samples in travel time, in seconds.
            self.dt = None
            #: float, The value on which the radar was triggering.
            self.trig_level = None

            # Per-trace attributes
            #: np.ndarray(tnum,) of the acquisition time of each trace
            #: note that this is referenced to Jan 1, 0 CE (matlabe datenum)
            #: for convenience, use the `datetime` attribute to access a python
            #: version of the day
            self.decday = None
            #: np.ndarray(tnum,) latitude along the profile. Generally not in
            #: projected coordinates
            self.lat = None
            #: np.ndarray(tnum,) longitude along the profile.
            #: Generally not in projected coords.
            self.long = None
            #: np.ndarray(tnum,) of the distances along the profile.
            #: units will depend on whether geographic coordinate transforms,
            #: as well as GPS data, are available.
            self.dist = None
            #: np.ndarray(tnum,) The pressure at acquisition.
            #: ImpDAR does not use this at present.
            self.pressure = None
            #: np.ndarray(tnum,) The 1-indexed number of the trace
            self.trace_num = None

            #: np.ndarray(tnum,) the index in each trace where the triggered.
            self.trig = None

            # Sample-wise attributes
            #: np.ndarray(snum,) The two way travel time to each sample, in us
            self.travel_time = None

            # Optional attributes
            #: str, the input filename. May be left as None.
            self.fn = None
            #: str, the projected coordinate system of the data
            self.t_srs = None
            #: np.ndarray(tnum,) Optional.
            #: Projected x-coordinate along the profile.
            self.x_coord = None
            #: np.ndarray(tnum,) Optional.
            #: Projected y-coordinate along the profile.
            self.y_coord = None
            #: np.ndarray(tnum,) Optional. Elevation along the profile.
            self.elev = None
            #: np.ndarray(tnum,) Optional.
            #: Depth of each trace below the surface
            self.nmo_depth = None

            # Special attributes
            #: impdar.lib.RadarFlags object containing information about the
            #: processing steps done.
            self.flags = RadarFlags()
            #: impdar.lib.Picks object with information about picks
            #: in/interpretation of the data.
            #: The init method of picks needs some basic data to calculate
            #: frequencies, etc, so it is not created until it is needed (maybe
            #: after some modifications to the data).
            self.picks = None

            self.data_dtype = None
            return

        mat = loadmat(fn_mat)
        for attr in self.attrs_guaranteed:
            # Exceptional case for 'data' variable because there are alternative names
            if attr == 'data':
                self._parse_stodeepdata(mat)
            elif attr not in mat:
                raise KeyError('.mat file does not appear to be in the StoDeep/ImpDAR format')
            else:
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
                elif mat[attr].shape[0] == 1 or (len(mat[attr].shape) > 1 and
                                                 mat[attr].shape[1] == 1):
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

    def _parse_stodeepdata(self, mat, data_attrs=STODEEP_ATTRS):
        """Set data attribute in a prioritized order."""
        data_dict = {}
        for data_attr in data_attrs:
            if data_attr in mat:
                if len(mat[data_attr].dtype) > 0:
                    print('Warning: Multiple arrays stored in {:s}, taking the first.'.format(data_attr))
                    data_dict[data_attr] = mat[data_attr][0][0][0]
                else:
                    data_dict[data_attr] = mat[data_attr]
        for i, attr in enumerate(data_attrs):
            if attr in data_dict:
                data_dict['data'] = data_dict[attr]
                if attr != 'data':
                    del data_dict[attr]
                if i > 0:
                    print('First priority data {:s} not in structure, using {:s}'.format(data_attrs[0], attr))
                    print('(caused a rename of {:s}'.format(attr))
                break
        else:
            raise KeyError('Data do not appear to be in StoDeep format')
        for attr, val in data_dict.items():
            setattr(self, attr, val)

    def check_attrs(self):
        """Check if required attributes exist.

        This is largely for development only; loaders should generally call
        this method last, so that they can confirm that they have defined the
        necessary attributes.

        Raises
        ------
        ImpdarError
            If any required attribute is None,
            or any optional attribute is fully absent
        """
        # fn is required but defined separately
        for attr in self.attrs_guaranteed + ['fn']:
            if not hasattr(self, attr):
                raise ImpdarError('{:s} is missing. \
                    It appears that this is an ill-defined \
                        RadarData object'.format(attr))
            if getattr(self, attr) is None:
                raise ImpdarError('{:s} is None. \
                    It appears that this is an ill-defined \
                        RadarData object'.format(attr))

        for attr in self.attrs_optional:
            if not hasattr(self, attr):
                raise ImpdarError('{:s} is missing. \
                    It appears that this is an ill-defined \
                        RadarData object'.format(attr))

        # Do some shape checks, but we need to be careful since
        # variable-surface will screw this up
        if (self.data.shape != (self.snum, self.tnum)) and (self.elev is None):
            raise ImpdarError('The data shape does not match \
                              the snum and tnum values!!!')

        if hasattr(self, 'nmo_depth') and (self.nmo_depth is not None):
            if (self.nmo_depth.shape[0] != self.snum) and (self.elev is None):
                raise ImpdarError('The nmo_depth shape does not match \
                                  the tnum value!!!')

        # We checked for existence, so we can just confirm that these have
        # the right length if they exist
        for attr in ['lat', 'long', 'pressure', 'trig', 'elev', 'dist',
                     'x_coord', 'y_coord', 'decday']:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                if (not hasattr(getattr(self, attr), 'shape')) or (
                        len(getattr(self, attr).shape) < 1):
                    if getattr(self, attr) == 0:
                        # This is just caused by None being weird with matlab
                        setattr(self, attr, None)
                    else:
                        if attr == 'trig':
                            self.trig = np.ones((self.tnum,), dtype=int) * int(self.trig)
                        else:
                            raise ImpdarError('{:s} needs to be a vector'.format(attr))
                elif getattr(self, attr).shape[0] != self.tnum:
                    raise ImpdarError('{:s} needs length tnum {:d}'.format(attr, self.tnum))

        if not hasattr(self, 'data_dtype') or self.data_dtype is None:
            self.data_dtype = self.data.dtype
        return

    def get_projected_coords(self, t_srs=None):
        """Convert to projected coordinates

        Parameters
        ----------
        t_srs: str, optional
            A text string accepted by GDAL (e.g. EPSG:3031)
            If None (default) use UTM.
        """
        if t_srs is not None:
            transform, self.t_srs = gpslib.get_conversion(t_srs=t_srs)
        else:
            transform, self.t_srs = gpslib.get_utm_conversion(np.nanmean(self.lat), np.nanmean(self.long))

        pts = np.array(transform(np.vstack((self.long, self.lat)).transpose()))

        self.x_coord, self.y_coord = pts[:, 0], pts[:, 1]
        self.dist = np.zeros((len(self.y_coord), ))
        self.dist[1:] = np.cumsum(np.sqrt(np.diff(self.x_coord) ** 2.0
            + np.diff(self.y_coord) ** 2.0)) / 1000.0

    def get_ll(self, s_srs):
        """Convert to projected coordinates

        Parameters
        ----------
        t_srs: str, optional
            A text string accepted by GDAL (e.g. EPSG:3031)
            If None (default) use UTM.
        """
        transform, self.t_srs = gpslib.get_rev_conversion(t_srs=s_srs)

        pts = np.array(transform(np.vstack((self.x_coord, self.y_coord)).transpose()))
        self.long, self.lat = pts[:, 0], pts[:, 1]

    @property
    def datetime(self):
        """Get pythonic version of the acquisition time of each trace."""
        return np.array([datetime.datetime(1970, 1, 1) +
                         datetime.timedelta(days=int(dd)) +
                         datetime.timedelta(days=dd % 1)
                         for dd in self.decday], dtype=np.datetime64)
