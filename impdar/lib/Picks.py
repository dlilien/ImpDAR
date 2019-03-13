#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the MIT license.

"""

"""
from .LastTrace import LastTrace
from .LeaderTrailer import LeaderTrailer
from .PickParameters import PickParameters


class Picks():
    """Information about picks"""
    attrs = ['samp1', 'samp2', 'samp3', 'time', 'power', 'picknums']
    spec_attrs = ['lasttrace', 'lt', 'pickparams']

    def __init__(self, radardata, pick_struct=None):
        if pick_struct is not None:
            # Loading from a file
            for attr in self.attrs:
                setattr(self, attr, pick_struct[attr][0][0])
                # Convert matlab zeros to Nones
                if getattr(self, attr) == np.zeros((1, 1)):
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
            self.picknums = [picknum]
            self.lasttrace.add_pick(-9999, 0)

        elif np.all(np.isnan(self.samp1[-1, :])):
            # If the last pick is blank, we just overwrite it. Zero the pick.
            self.mod_line(self.samp1.shape[0], 0, 0)
        else:
            # We are just adding a row to the existing matrices of samples etc.
            self.samp1 = np.vstack((self.samp1, np.zeros((1, self.radardata.tnum))))
            self.samp2 = np.vstack((self.samp2, np.zeros((1, self.radardata.tnum))))
            self.samp3 = np.vstack((self.samp3, np.zeros((1, self.radardata.tnum))))
            self.time = np.vstack((self.time, np.zeros((1, self.radardata.tnum))))
            self.power = np.vstack((self.power, np.zeros((1, self.radardata.tnum))))
            self.lasttrace.add_pick(-9999, 0)
            self.picknums.append(picknum)
        # We return the row number of the sample, which gives access to all its info
        return self.samp1.shape[0]

    def update_pick(self, picknum, pick_info):
        self.samp1[picknum, :] = pick_info[0, :]
        self.samp2[picknum, :] = pick_info[1, :]
        self.samp3[picknum, :] = pick_info[2, :]
        self.time[picknum, :] = pick_info[3, :]
        self.power[picknum, :] = pick_info[4, :]

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
