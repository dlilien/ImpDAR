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

def continuityKarlsson(P,s_ind,b_ind,lat,lon,cutoff_ratio,win=20,uice=168.,eps=3.2):
    """
    Karlsson Continuity Method

    Based on Karlsson et al. (2012)
    This method gives a value for the continuity of radar layers

    Parameters
    ----------
    P:              uncorrected power
    s_ind:          surface pick index
    b_ind:          bed pick index
    cutoff_ratio:   assigns the number of samples that are removed from top and bottom of the trace

    Output
    ---------
    cont:
    cont_filt:

    """

    # empty continuity index array
    cont = np.empty_like(b_ind).astype(float)
    cont[:] = np.nan
    # calculate the continuity index for each trace
    for tr in range(len(P[0])):
        spick=int(s_ind[tr])
        bpick=int(b_ind[tr])
        if bpick-spick<10 or bpick>len(P[:,0]) or np.isnan(bpick-spick):
            continue
        else:
            # get data from between the surface and bed
            p_ext=P[spick:bpick,tr]
            # cutoff based on the assigned ratio
            cut=int(len(p_ext)/cutoff_ratio)
            p_ext=p_ext[cut:-cut]
            if np.any(~np.isfinite(p_ext)):
                continue
            # calculate the continuity index based on Karlsson et al. (2012) eq. 1
            cont[tr]=np.mean(abs(np.gradient(p_ext)))
    # smoother--simple moving boxcar
    cont_filt = np.convolve(cont, np.ones((win,))/win, mode='valid')

    return cont,cont_filt

