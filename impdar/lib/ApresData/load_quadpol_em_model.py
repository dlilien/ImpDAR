#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Load quad-polarized model output

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Feb 4 2021
"""

import numpy as np
from datetime import datetime
from . import QuadPolData

# -----------------------------------------------------------------------------------------------------

def load_quadpol_fujita(model):
    """
    Load processed modeled apres profiles from all four polarizations: hh, hv, vh, vv
    into one data object

    Parameters
    ----------
    model : class
        output from the effective medium model

    Returns
    -------
    quadpol_data: class
        quad-polarized apres data object
    """

    # load into the QuadPolData object
    quadpol_data = QuadPolData(None)
    quadpol_data.shh = model.shh
    quadpol_data.shv = model.shv
    quadpol_data.svh = model.svh
    quadpol_data.svv = model.svv
    quadpol_data.range = model.range

    now = datetime.now()
    timezero = datetime(1, 1, 1, 0, 0, 0)
    offset = now-timezero
    quadpol_data.decday = offset.days + offset.seconds/(3600.*24.) + 377. # Matlab compatable

    quadpol_data.snum = len(model.shh)
    v = model.c/np.sqrt(model.epsr)
    quadpol_data.travel_time = quadpol_data.range/v
    quadpol_data.dt = np.mean(np.gradient(quadpol_data.travel_time))

    quadpol_data.data_dtype = quadpol_data.shh.dtype

    return quadpol_data
