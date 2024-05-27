#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
"""Legacy class for appending bits of neighboring profiles."""

from .Crop import Crop


class LeaderTrailer():
    """The lt structure from StoInterpret.

    Parameters
    ----------
    radardata: impdar.lib.RadarData.RadarData
        The radardata object for which this contains info.
    lt_struct: dict
        Input dictionary of information that we turn into the struct.
    """

    attrs = ['llength', 'tlength', 'ltmatrix']

    def __init__(self, radardata, lt_struct=None):
        if lt_struct is not None:
            for attr in self.attrs:
                setattr(self, attr, lt_struct[0][0][attr])
            self.crop = Crop(radardata)
        else:
            self.llength = 0
            self.tlength = 0
            self.ltmatrix = 0
            self.crop = Crop(radardata)

    def to_struct(self):
        """Export to dictionary for Matlab.

        Returns
        -------
        dict:
            The data in a format for use with scipy.io.savemat
        """
        mat = {}
        for attr in self.attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                mat[attr] = 0
        mat['crop'] = self.crop.to_struct()
        return mat
