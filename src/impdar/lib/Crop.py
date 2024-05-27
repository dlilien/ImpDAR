#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
"""Object to track modifications to data dimensions."""

import numpy as np


class Crop():
    """Crop information that tells us how the data have been modified.

    Parameters
    ----------
    radardata: impdar.lib.RadarData.RadarData
        Radardata object for this to be affiliated with

    Attributes
    ----------
    tnum: int
        The number of traces
    maxsnum: int
        The maximum snum
    mintt: float
        The minimum twtt
    maxstt: float
        The maximum twtt
    """

    attrs = ['tnum', 'maxsnum', 'mintt', 'maxtt']

    def __init__(self, radardata):
        self.tnum = radardata.tnum
        self.maxsnum = radardata.snum
        self.mintt = np.min(radardata.travel_time)
        self.maxtt = np.max(radardata.travel_time)

    def to_struct(self):
        """Export to format for matlab.

        Returns
        -------
        mat: dict
            Attributes in form digestible for scipy.io.savemat
        """
        mat = {}
        for attr in self.attrs:
            # I don't guard against None here since None should not be allowed
            mat[attr] = getattr(self, attr)
        return mat
