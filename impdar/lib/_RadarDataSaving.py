#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
A class that RadarData inherits so that it has convenient saving methods
"""

from scipy.io import savemat
from .RadarFlags import RadarFlags


class RadarDataSaving:

    def save(self, fn):
        """Save the radar data

        Parameters
        ----------
        fn: str
            Filename. Should have a .mat extension
        """
        mat = {}
        for attr in self.attrs_guaranteed:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                # this guards against error in matlab format
                mat[attr] = 0
        for attr in self.attrs_optional:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
        if hasattr(self, 'picks') and self.picks is not None:
            mat['picks'] = self.picks.to_struct()
        if self.flags is not None:
            mat['flags'] = self.flags.to_matlab()
        else:
            # We want the structure available to prevent read errors from corrupt files
            mat['flags'] = RadarFlags().to_matlab()
        savemat(fn, mat)
