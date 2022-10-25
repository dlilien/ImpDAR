#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Load ApRES data as a profile

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

July 11 2022
"""

import glob
import numpy as np

import numpy as np
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags

from ..ApresData.load_apres import load_apres
from ..ImpdarError import ImpdarError


def load_apres_profile(fns_apres, *args, **kwargs):
    """
    Load Apres acquisitiuon into a profile.
    """

    apres_obj = load_apres(fns_apres)
    apres_obj.apres_range(2)

    # Setup a RadarData class to fill
    apres_data = RadarData(None)
    apres_data.fn = fns_apres[0]

    # Try to pull in as much as possible directly from the apres object
    for attr in vars(apres_obj):
        if attr == 'data':
            continue
        if attr in apres_data.attrs_guaranteed:
            setattr(apres_data, attr, getattr(apres_obj,attr))
        if attr in apres_data.attrs_optional:
            setattr(apres_data, attr, getattr(apres_obj,attr))

    # Reshape the data array into a profile
    apres_data.data = np.reshape(apres_obj.data,(apres_obj.bnum*apres_obj.cnum, apres_obj.snum))
    apres_data.data = np.transpose(apres_data.data)
    apres_data.data.dtype = complex
    apres_data.snum = apres_data.data.shape[0]
    apres_data.tnum = apres_data.data.shape[1]
    apres_data.trace_num = np.arange(apres_data.tnum)

    # travel time in the apres object is set up for fmcw data
    # convert using the range attr
    apres_data.travel_time = apres_obj.Rcoarse/(apres_obj.header.ci/2.)
    apres_data.travel_time *= 1e6

    # reformat arrays for values that are only taken once per chirp
    apres_data.decday = apres_obj.chirp_time.flatten()
    apres_data.lat = np.transpose(np.tile(apres_obj.lat,
                                          (apres_obj.cnum, 1))).flatten()
    apres_data.long = np.transpose(np.tile(apres_obj.long,
                                           (apres_obj.cnum, 1))).flatten()
    if apres_obj.elev is None:
        apres_data.elev = np.zeros_like(apres_data.lat)
    elif np.shape(apres_obj.elev) == np.shape(apres_obj.lat):
        apres_data.elev = np.transpose(np.tile(apres_obj.elev,
                                           (apres_obj.cnum, 1))).flatten()

    # Set the along-track coordinates
    try:
        apres_data.get_projected_coords()
    except:
        apres_data.dist = np.zeros(apres_data.tnum)
    apres_data.trace_int = np.gradient(apres_data.dist)

    # fill and empty pressure array
    apres_data.pressure = np.zeros(apres_data.tnum)
    apres_data.trig = np.nan*np.zeros(apres_data.tnum)
    apres_data.trig_level = np.nan*np.zeros(apres_data.tnum)
    apres_data.chan = 0

    # Initialize the flags
    apres_data.flags = RadarFlags()

    # Check attributed in the RadarData class and return
    apres_data.check_attrs()
    return apres_data
