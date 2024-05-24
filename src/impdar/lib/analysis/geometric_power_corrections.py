#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""

Power Corrections from the glaciological literature.

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 26 2019

"""

import numpy as np

# ----------------------------------------------------------------------------

def power_correction(dat,eps=[],d_eps=[],u=1.69e8,h_aircraft=0.):
    """
    Geometric spreading correction for radar power.
    Optionally includes refractive focusing

    Parameters
    ---------
    eps:    array; optional
        permittivity (relative)
    d_eps:  array; optional
        depths for permittivity boundaries
    u:      float; optional
        speed of light in ice
    h_aircraft: float; optional
        height of aircraft, airborne surveys need a correction for refractive focusing from air to ice


    Output
    ---------
    corrected_power
    """

    # get a pick depth
    if 'z' in vars(dat.picks):
        Z = dat.picks.z
    else:
        print('Warning: setting pick depth for constant velocity in ice.')
        Z = dat.picks.time*u/2./1e6

    # spreading correction for a spherical wave
    spherical_loss = (2.*Z)**2.

    q = np.ones_like(Z)
    if len(d_eps) > 0:
        if d_eps[0] != 0:
            raise KeyError('The first depth needs to be 0.')
        # correct for focusing from air to firn
        if h_aircraft > 0.:
            qadd = refractive_focusing(h_aircraft,2.*(Z+h_aircraft),1.,eps[0])
            q*=qadd
        # correct for focusing within the firn
        for i in range(len(eps)-1):
            qadd = refractive_focusing(d_eps[i],2.*Z,eps[i],eps[i+1])
            q*=qadd

    # power correction including spreading and refractive gains
    dat.picks.corrected_power = dat.picks.power * spherical_loss/q

# ----------------------------------------------------------------------------

def refractive_focusing(z1,z2,eps1,eps2):
    """
    Refractive focusing at an interface
    Bogorodsky et al., 1985; equation 3.8

    Parameters
    ---------
    z1:     float
        scalar      Thickness above interface (m)
    z2:     scalar      Thickness below interface (m)
    eps1:   scalar      Permittivity above interface (relative)
    eps2:   scalar      Permittivity below interface (relative)

    Output
    ---------
    q:      scalar              refractive focusing coefficient
    """
    q = ((z1+z2)/(z1+z2*np.sqrt(eps1/eps2)))**2.
    if hasattr(q,'__len__'):
        q[z2 <= z1] = 1.
    else:
        if z2 <= z1:
            q = 1.
    return q
