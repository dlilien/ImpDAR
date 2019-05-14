#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
"""
import numpy as np
from .LastTrace import LastTrace
from .LeaderTrailer import LeaderTrailer
from .PickParameters import PickParameters


class Picks():
    """Information about picks

    This object holds all picks for a given radargram. The main containers are matrices holding information about the indices in the data matrix a given pick lies.

    Attributes
    ----------
    samp1: nsamp x tnum array
        Min/max above the center of each pick. A new row is added for each new pick.
    samp2: nsamp x tnum array
        Min/max at the center of each pick. A new row is added for each new pick.
    samp3: nsamp x tnum array
        Min/max below the center of each pick. A new row is added for each new pick.
    time: nsamp x tnum array
        In StoDeep used to contain TWTT samp2. Since this is redundant I'm deptrecating it, and it will be zeros for impdar processed data.
    power: nsamp x tnum array
        Power across a pick.
    picknums: list of length nsamp
        The number of each pick.
    lasttrace: impdar.lib.LastTrace.LastTrace
        Information about the end of the trace
    lt: impdar.lib.LeaderTrailer.LeaderTrailer
        StoDeep legacy for compatibility, unused
    pickparams: impdar.lib.PickParameters.PickParameters
        This structure contains important information used in picking, such as frequency for picks.
    """
    attrs = ['samp1', 'samp2', 'samp3', 'time', 'power', 'picknums']
    flatten = [False, False, False, False, False, True]
    spec_attrs = ['lasttrace', 'lt', 'pickparams']

    def __init__(self, radardata, pick_struct=None):
        if pick_struct is not None:
            # Loading from a file
            for attr, flat in zip(self.attrs, self.flatten):
                setattr(self, attr, pick_struct[attr][0][0])
                if flat:
                    setattr(self, attr, getattr(self, attr).flatten())

                # Convert matlab zeros to Nones
                if getattr(self, attr).shape == (1, 1) and getattr(self, attr)[0][0] == 0:
                    setattr(self, attr, None)
            self.lasttrace = LastTrace(pick_struct['lasttrace'])
            self.lt = LeaderTrailer(radardata, pick_struct['lt'])
            self.pickparams = PickParameters(radardata, pick_struct['pickparams'])
            self.picknums = self.picknums.tolist()
        else:
            # Blank initialization
            self.samp1 = None
            self.samp2 = None
            self.samp3 = None
            self.time = None
            self.power = None
            self.picknums = None
            self.lasttrace = LastTrace()
            self.lt = LeaderTrailer(radardata)
            self.pickparams = PickParameters(radardata)

        # These are not in StoInterpret, but I'm using them to keep life more object oriented
        self.radardata = radardata
        # This will contain the handles for all the lines plotted so we can do some selection
        self.lines = []

    def add_pick(self, picknum=0):
        """Add a new pick

        This method handles any complexity in adding a new pick. If no picks exist, it creates the matrices. Otherwise, it just adds rows. If the last pick hasn't been used, this just recycles that pick.

        Parameters
        ----------
        picknum: int, optional
            Number to call the new pick. Default zero

        Returns
        -------
        index: int
            The index of the row number containing the new pick

        Raises
        ------
        ValueError if the picknum already exists--we do not deal with repeats

        """
        if self.samp1 is None:
            # We have no matrices yet
            self.samp1 = np.zeros((1, self.radardata.tnum))
            self.samp2 = np.zeros((1, self.radardata.tnum))
            self.samp3 = np.zeros((1, self.radardata.tnum))
            self.time = np.zeros((1, self.radardata.tnum))
            self.power = np.zeros((1, self.radardata.tnum))
            self.samp1[0, :] = np.nan
            self.samp2[0, :] = np.nan
            self.samp3[0, :] = np.nan
            self.time[0, :] = np.nan
            self.power[0, :] = np.nan
            self.picknums = [picknum]
            self.lasttrace.add_pick(-9999, 0)
        elif np.all(np.isnan(self.samp1[-1, :])):
            # If the last pick is blank, we just overwrite it. Zero the pick.
            self.samp1[-1, :] = np.nan
            self.samp2[-1, :] = np.nan
            self.samp3[-1, :] = np.nan
            self.time[-1, :] = np.nan
            self.power[-1, :] = np.nan
            self.picknums[-1] = picknum
        else:
            # If loading from matlab, need to cast picknums as a list
            if type(self.picknums) == np.ndarray:
                self.picknums = self.picknums.flatten().tolist()
            if picknum in self.picknums:
                raise ValueError('We already have that pick')

            # We are just adding a row to the existing matrices of samples etc.
            self.samp1 = np.vstack((self.samp1, np.zeros((1, self.radardata.tnum))))
            self.samp2 = np.vstack((self.samp2, np.zeros((1, self.radardata.tnum))))
            self.samp3 = np.vstack((self.samp3, np.zeros((1, self.radardata.tnum))))
            self.time = np.vstack((self.time, np.zeros((1, self.radardata.tnum))))
            self.power = np.vstack((self.power, np.zeros((1, self.radardata.tnum))))
            self.samp1[-1, :] = np.nan
            self.samp2[-1, :] = np.nan
            self.samp3[-1, :] = np.nan
            self.time[-1, :] = np.nan
            self.power[-1, :] = np.nan
            self.lasttrace.add_pick(-9999, 0)

            self.picknums.append(picknum)
        # We return the row number of the sample, which gives access to all its info
        return self.samp1.shape[0]

    def update_pick(self, picknum, pick_info):
        """Update a pick with new information

        Rather than go in and manually update a pick every time it is changed, we take in all information about an individual pick simultaneously and this method updates the pick's information in the Pick object.

        Parameters
        ----------
        picknum: int
            The pick number to update. Must be in picknums
        pick_info: 5xtnum np.ndarray
            Array where rows are upper pick bound, pick center, lower pick bound, time (deprecated, usually left as zeros or nans), and power across the pick.

        Raises
        ------
        ValueError if picknum is not in picknums or if the shape of the pick_info is bad.

        """
        try:
            ind = self.picknums.index(picknum)
        except ValueError:
            raise ValueError('picknum provided is not a pick; you must you use a picknum not an index')

        if pick_info.shape != (5, self.radardata.tnum):
            raise ValueError('pick_info must be a 5xtnum array')

        self.samp1[ind, :] = pick_info[0, :]
        self.samp2[ind, :] = pick_info[1, :]
        self.samp3[ind, :] = pick_info[2, :]
        self.time[ind, :] = pick_info[3, :]
        self.power[ind, :] = pick_info[4, :]

    def to_struct(self):
        """Convert to a format writable to a .mat file

        Returns
        -------
        mat: dict
            Dictionary of attributes with special types converted to be matlab compatible
        """
        mat = {}
        for attr in self.attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                mat[attr] = 0
        for attr in self.spec_attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr).to_struct()
            else:
                mat[attr] = 0
        return mat
