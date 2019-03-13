#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
LastTrace
"""


class LastTrace():
    """The sample and trace of the last trace for picking"""
    attrs = ['snum', 'tnum']

    def __init__(self, lasttrace_struct=None):
        if lasttrace_struct is not None:
            for attr in self.attrs:
                setattr(self, attr, lasttrace_struct[0][0][attr][0][0].flatten())
        else:
            self.snum = None
            self.tnum = None

    def add_pick(self, snum, tnum):
        if self.snum is None:
            self.snum = [snum]
            self.tnum = [tnum]
        else:
            self.snum.append(snum)
            self.tnum.append(tnum)

    def mod_line(self, ind, snum, tnum):
        self.snum[ind] = snum
        self.tnum[ind] = tnum

    def to_struct(self):
        mat = {}
        for attr in self.attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                mat[attr] = 0
        return mat
