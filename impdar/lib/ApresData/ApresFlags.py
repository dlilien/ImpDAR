#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Flags to keep track of processing steps
"""

import numpy as np

class ApresFlags():
    """Flags that indicate the processing that has been used on the data.

    These are used for figuring out whether different processing steps have been performed. They also contain some information about the input arguments for some (but not all) of the processing steps.

    Attributes
    ----------
    batch: bool
        Legacy indication of whether we are batch processing. Always False.
    agc: bool
        Automatic gain control has been applied.
    reverse: bool
        Data have been reversed.
    restack: bool
        Data have been restacked.
    rgain: bool
        Data have a linear range gain applied.
    bpass: 3x1 :class:`numpy.ndarray`
        Elements: (1) 1 if bandpassed; (2) Low; and (3) High (MHz) bounds
    hfilt: 2x1 :class:`numpy.ndarray`
        Elements: (1) 1 if horizontally filtered; (2) Filter type
    interp: 2x1 :class:`numpy.ndarray`
        Elements: (1) 1 if constant distance spacing applied (2) The constant spacing (m)
    mig: 2x1 :class: String
        None if no migration done, mtype if migration done.
     """

    def __init__(self):
        self.file_read_code = None
        self.range = 0
        self.stack = 1
        self.attrs = ['file_read_code','phase2range','stack']
        self.attr_dims = [None,None,None]

    def to_matlab(self):
        """Convert all associated attributes into a dictionary formatted for use with :func:`scipy.io.savemat`
        """
        outmat = {att: getattr(self, att) for att in self.attrs}
        return outmat

    def from_matlab(self, matlab_struct):
        """Associate all values from an incoming .mat file (i.e. a dictionary from :func:`scipy.io.loadmat`) with appropriate attributes
        """
        for attr, attr_dim in zip(self.attrs, self.attr_dims):
            setattr(self, attr, matlab_struct[attr][0][0][0])
            # Use this because matlab inputs may have zeros for flags that
            # were lazily appended to be arrays, but we preallocate
            if attr_dim is not None and getattr(self, attr).shape[0] == 1:
                setattr(self, attr, np.zeros((attr_dim, )))
