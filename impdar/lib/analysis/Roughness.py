#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""

Roughness calculation for picked layers.

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 26 2019

"""

import numpy as np
from scipy.signal import detrend,medfilt
from scipy.special import i0

def kirchhoff_roughness(dat,picknum,freq,filt_n=101,eps=3.15):
    """
    Roughness by Kirchhoff Theory
    Christianson et al. (2016), equation C2

    Paramaters
    ----------
    freq:   float
        antenna frequency
    filt_n: int; optional
        number of traces included in the median filter
    eps:    float; optional
        relative permittivity of ice
    """

    if 'interp' not in vars(dat.flags):
        raise KeyError('Do interpolation before roughness calculation.')

    # calculate the speed and wavelength
    eps0 = 8.8541878128e-12        # vacuum permittivity
    mu0 = 1.25663706212e-6         # vacuum permeability
    u = 1./np.sqrt(eps*eps0*mu0)   # speed of light in ice
    lam = u/freq                   # wavelength m

    # get a pick depth
    if 'z' in vars(dat.picks):
        Z = dat.picks.z
    else:
        print('Warning: setting pick depth for constant velocity in ice.')
        Z = dat.picks.time*u/2/1e6

    # Find window size based on the width of the first Fresnel zone
    D1 = np.sqrt(2.*lam*(np.nanmean(Z)/np.sqrt(eps))) # Width of Fresnel zone
    dx = dat.trace_int[0]                                       # m spacing between traces
    N = int(round(D1/(2.*dx)))                                   # number of traces in the Fresnel window

    # -----------------------------------------------------------------------------

    # Define the bed geometry
    bed_raw = dat.elev - Z[picknum]
    bed_filt = medfilt(bed_raw,filt_n)

    # RMS bed roughness; Christianson et al. (2016) equation C2
    ED1 = np.nan*np.empty((len(bed_filt),))
    for n in range(N,len(bed_filt)-N+1):
        b = bed_filt[n-N:n+N].copy()
        b = b[np.where(~np.isnan(b))]
        if len(b) <= 1:
            ED1[n] = np.nan
        else:
            b_ = detrend(b)
            b_sum = 0.
            for i in range(len(b)):
                b_sum += (b_[i])**2.
            ED1[n] = np.sqrt((1/(len(b)-1.))*b_sum)


    # Find the power reduction by Kirchoff theory
    # Christianson et al. (2016), equation C1

    g = 4.*np.pi*ED1/lam
    b = (i0((g**2.)/2.))**2.
    pn = np.exp(-(g**2.))*b

    return ED1,pn
