#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

Load St Olaf .mat files and convert to the .mat ImpDAR format

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

March 31 2020

"""

import numpy as np
from scipy.io import loadmat
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags

def load_stomat(fn_sto, dname=None, *args, **kwargs):
    """Load a St Olaf .mat file into ImpDAR"""

    # Import the matlab file and setup a RadarData class to fill
    sto_mat = loadmat(fn_sto)
    sto_data = RadarData(None)
    sto_data.fn = fn_sto

    # Get variables that are already in the St Olaf .mat file
    sto_data.dt = sto_mat['dt'][0][0]
    sto_data.chan = sto_mat['chan'][0][0]
    sto_data.snum = sto_mat['snum'][0][0]
    sto_data.tnum = sto_mat['tnum'][0][0]
    sto_data.trace_num = sto_mat['trace_num'][0]
    sto_data.trig_level = sto_mat['trig_level'][0]
    sto_data.pressure = sto_mat['pressure'][0]
    sto_data.travel_time = sto_data.dt * 1.0e6 * np.arange(sto_data.snum)
    sto_data.lat = sto_mat['lat'][0]
    sto_data.long = sto_mat['long'][0]
    sto_data.elev = sto_mat['elev'][0]
    sto_data.decday = sto_mat['decday'][0]
    sto_data.trace_int = sto_mat['trace_int'][0]
    sto_data.dist = sto_mat['dist'][0]
    # Try the x/y coordinates but not neccesarry
    try:
        sto_data.x_coord = sto_mat['x_coord'][0]
        sto_data.y_coord = sto_mat['y_coord'][0]
    except:
        Warning('No projected coordinate system (x,y).')
    # Sometimes the trigger is saved as a number and sometimes as an array
    # ImpDAR wants it as an array
    trig = sto_mat['trig'][0]
    if len(trig) == sto_data.tnum:
        sto_data.trig = trig
    elif len(trig) == 1:
        sto_data.trig = trig[0]*np.ones(sto_data.tnum)

    # Get the image, could be saved under a variety of names
    if dname is not None:
        sto_data.data = sto_mat[dname]
    else:
        dnames = ['filtdata','interp_data','nmo_dat','migdata']
        for name in dnames:
            if name in sto_mat.keys():
                print('Found data image saved under:',name)
                sto_data.data = sto_mat[name]
            else:
                continue

    # TODO: Import the flags
    """
    sto_data.flags = RadarFlags()
    sto_flags = sto_mat['flags']
    'bpass'
    'hfilt'
    'rgain'
    'agc'
    'restack'
    'reverse'
    'crop'
    'nmo'
    'interp'
    'mig'
    'elev'
    """

    # Check attributed in the RadarData class and return
    sto_data.check_attrs()
    return sto_data
