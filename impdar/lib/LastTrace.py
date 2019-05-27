#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.


import numpy as np


class LastTrace():
    """The sample and trace of the last trace when picking
    
    Attributes
    ----------
    snum: int
        The last sample picked
    tnum: int
        The last trace picked
    """
    attrs = ['snum', 'tnum']

    def __init__(self, lasttrace_struct=None):
        if lasttrace_struct is not None:
            for attr in self.attrs:
                val = lasttrace_struct[0][0][attr][0][0].flatten()
                if len(val) == 1 and val[0] == -9999:
                    val = None
                setattr(self, attr, val)
        else:
            self.snum = None
            self.tnum = None

    def add_pick(self, snum, tnum):
        if self.snum is None:
            self.snum = [snum]
            self.tnum = [tnum]
        else:
            # we may load these as arrays and then want to append
            if type(self.snum) == np.ndarray:
                self.snum = self.snum.flatten().tolist()
            if type(self.tnum) == np.ndarray:
                self.tnum = self.tnum.flatten().tolist()
            
            # If we have started anew, we only need these lines
            self.snum.append(int(snum))
            self.tnum.append(int(tnum))

    def mod_line(self, ind, snum, tnum):
        if (self.snum is None) or (self.tnum is None):
            raise AttributeError('need snum and tnum defined')
        if len(self.snum) <= ind or len(self.snum) <= ind:
            raise ValueError('Index is too large for snum/tnum')
        self.snum[ind] = snum
        self.tnum[ind] = tnum

    def to_struct(self):
        mat = {}
        for attr in self.attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                mat[attr] = -9999
        return mat
