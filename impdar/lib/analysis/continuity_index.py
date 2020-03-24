#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Continuity index for layers in radar data.
Karlsson et al. (2012)

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 26 2019

"""

import numpy as np

# ----------------------------------------------------------------------------

def continuity_index(dat,b_ind,s_ind=None,cutoff_ratio=None):
    """
    Karlsson Continuity Method

    Based on Karlsson et al. (2012)
    This method gives a value for the continuity of radar layers

    Parameters
    ----------
    b_ind:  int
        bed pick index
    s_ind:  int; optional
        surface pick index
    cutoff_ratio:   float; optional
        assigns the number of samples that are removed from top and bottom of the trace

    Output
    ---------
    conttinuity_index: array
    """

    P = 10*np.log10(dat.data**2.)

    bpick = dat.picks.samp1[b_ind]
    if s_ind is None:
        spick = np.zeros_like(bpick)
    else:
        spick = dat.picks.samp1[s_ind]

    # empty continuity index array
    cont = np.empty((dat.tnum,)).astype(float)
    # calculate the continuity index for each trace
    for tr in range(dat.tnum):
        # Nan if the picks are nan
        if np.isnan(bpick[tr]) or np.isnan(spick[tr]):
            cont[tr] = np.nan
        else:
            # get data from between the surface and bed
            b = int(bpick[tr])
            s = int(spick[tr])
            p_ext=P[s:b,tr]
            # cutoff based on the assigned ratio
            if cutoff_ratio is not None:
                cut=int(len(p_ext)*cutoff_ratio)
                p_ext=p_ext[cut:-cut]

            # Nan if sampling criteria are not met
            if len(p_ext) < 10 or len(p_ext) > dat.snum or np.any(~np.isfinite(p_ext)):
                cont[tr] = np.nan

            # calculate the continuity index based on Karlsson et al. (2012) eq. 1
            else:
                cont[tr]=np.mean(abs(np.gradient(p_ext)))

    dat.continuity_index = cont
