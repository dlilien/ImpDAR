#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Crop information
"""

import numpy as np


class Crop():
    """Crop information. I have no idea what this is for but it is retained for backwards compatibility"""
    attrs = ['tnum', 'maxsnum', 'mintt', 'maxtt']

    def __init__(self, radardata):
        self.tnum = radardata.tnum
        self.maxsnum = radardata.snum
        self.mintt = np.min(radardata.travel_time)
        self.maxtt = np.max(radardata.travel_time)

    def to_struct(self):
        mat = {}
        for attr in self.attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                mat[attr] = 0
        return mat
