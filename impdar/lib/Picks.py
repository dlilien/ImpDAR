#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the MIT license.

"""

"""
import numpy as np
from .LastTrace import LastTrace
from .LeaderTrailer import LeaderTrailer
from .PickParameters import PickParameters


class Picks():
    """Information about picks"""
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
