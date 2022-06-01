#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Benjamin Hills <bhills@uw.edu>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import os

import numpy as np
import h5py

from scipy.io import savemat
from .ApresFlags import ApresFlags
from .ApresHeader import ApresHeader


def save(self, fn):
    """Save the radar data

    Parameters
    ----------
    self: class
        ApresData object
    fn: str
        Filename. Extension can be h5 or legacy mat.
    """

    ext = os.path.splitext(fn)[1]
    if ext in ['.h5', '.hdf5']:
        return save_h5(self, fn)
    elif ext == '.mat':
        return save_mat(self, fn)
    else:
        raise ValueError('File extension choices are .h5 and .mat (legacy)')


def save_mat(self, fn):
    """Save the radar data as an ImpDAR .mat file

    Parameters
    ----------
    self: class
        ApresData object
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

    if 'header' in vars(self):
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


def save_h5(self, fn):
    """Save the radar data as an h5 file

    Parameters
    ----------
    self: class
        ApresData object
    fn: str
        Filename. Should have a .h5 extension
    """
    with h5py.File(fn, 'w') as f:
        save_as_h5_group(self, f, 'dat')


def save_as_h5_group(self, h5_file_descriptor, groupname='dat'):
    """Save to a group in h5 file (useful to group multiple datasets to one file

    Parameters
    ----------
    self: class
        ApresData object
    h5_file_descriptor: open hdf5 file
        The file object to write to. Can also be a group (if data should be a subgroup).
    groupname: str, optional
        The name this (sub)group should have.
    """
    grp = h5_file_descriptor.create_group(groupname)
    for attr in self.attrs_guaranteed:
        val = getattr(self, attr)
        if val is not None:
            if hasattr(val, 'shape') and np.any([s != 1 for s in val.shape]):
                if val.dtype == 'O':
                    if hasattr(self, 'data_dtype') and self.data_dtype is not None:
                        dtype = self.data_dtype
                    else:
                        dtype = 'f'
                else:
                    dtype = val.dtype
                grp.create_dataset(attr, data=val, dtype=dtype)
            else:
                grp.attrs.create(attr, val)
        else:
            grp.attrs[attr] = h5py.Empty('f')
    for attr in self.attrs_optional:
        if hasattr(self, attr) and getattr(self, attr) is not None:
            val = getattr(self, attr)
            if hasattr(val, 'shape') and np.any([s != 1 for s in val.shape]):
                if val.dtype == 'O':
                    if hasattr(self, 'data_dtype') and self.data_dtype is not None:
                        dtype = self.data_dtype
                    else:
                        dtype = 'f'
                else:
                    # override the dtype for data
                    if attr == 'data':
                        dtype = self.data_dtype
                    dtype = val.dtype
                grp.create_dataset(attr, data=val, dtype=dtype)
            else:
                grp.attrs.create(attr, val)

    if self.flags is not None:
        self.flags.write_h5(grp)
    else:
        ApresFlags().write_h5(grp)

    if self.header is not None:
        self.header.write_h5(grp)
    else:
        ApresHeader().write_h5(grp)
