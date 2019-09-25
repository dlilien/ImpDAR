#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""

import numpy as np
from scipy.io import savemat
from .ApresFlags import ApresFlags
from .ApresHeader import ApresHeader

def save_apres(self, fn):
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

    if self.flags is not None:
        mat['flags'] = self.flags.to_matlab()
    else:
        # We want the structure available to prevent read errors from corrupt files
        mat['flags'] = ApresFlags().to_matlab()

    if self.header is not None:
        mat['header'] = self.header.to_matlab()
    else:
        # We want the structure available to prevent read errors from corrupt files
        mat['header'] = ApresHeader().to_matlab()

    # Make sure not to expand the size of the data due to type conversion
    if hasattr(self, 'data_dtype') and self.data_dtype is not None and self.data_dtype != mat['data'].dtype:
        # Be carefuly of obliterating NaNs
        # We will use singles instead of ints for this guess
        if (self.data_dtype in [int, np.int8, np.int16]) and np.any(np.isnan(mat['data'])):
            print('Warning: new file is float16 rather than ', self.data_dtype, ' since we now have NaNs')
            mat['data'] = mat['data'].astype(np.float16)
        elif (self.data_dtype in [np.int32]) and np.any(np.isnan(mat['data'])):
            print('Warning: new file is float32 rather than ', self.data_dtype, ' since we now have NaNs')
            mat['data'] = mat['data'].astype(np.float32)
        elif (self.data_dtype in [np.int64]) and np.any(np.isnan(mat['data'])):
            print('Warning: new file is float64 rather than ', self.data_dtype, ' since we now have NaNs')
            mat['data'] = mat['data'].astype(np.float64)
        else:
            mat['data'] = mat['data'].astype(self.data_dtype)
    savemat(fn, mat)

