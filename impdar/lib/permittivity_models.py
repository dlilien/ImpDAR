#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

Permittivity model options for snow and firn densities
Implemented in ImpDAR

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

July 25 2019

"""

# ----------------------------------------------------------------------------

def snow_permittivity(rho,fs,m=0.,fw=9.07e9):
    """
    Calculate the dielectric permittivity of snow
    Kendra et al. (1998), IEEE

    Parameters
    ---------
    rho:        snow density       (g/cm3)
    fs:         radar frequency     (Hz)
    m:          snow wetness        (%)
    fw:         relaxation frequency of water at 0C

    Output
    ---------
    eps_ds:     snow permittivity
    """

    # permittivity of dry snow, Kendra eq. 13
    eps_s = 1. + 1.7*rho + 0.7*rho**2.
    # added permittivity from wetness, Kendra eq. 14
    eps_s += 0.02*m**1.015+(.073*m**1.31)/(1+(fs/fw))

    return eps_s

def firn_permittivity(rhof,rhoi=917.,epsi_real=3.12,epsi_imag=-9.5):
    """
    Calculate the dielectric permittivity of firn with the DECOMP mixing model
    Wilhelms (2005), GRL

    Parameters
    ---------
    rhof:           firn density       (kg/m3)
    rhoi:           ice density       (kg/m3)
    epsi_real:      real permittivity of ice (relative)
    epsi_imag:      imaginary permittivity of ice (relative)

    Output
    ---------
    eps_f:          firn permittivity   (relative)
    """

    # Wilhelms (2005), end of section 2
    lhs = 1. + (rhof/rhoi)*((epsi_real-1j*epsi_imag)**(1/3.)-1)
    eps_f = lhs**3.

    return eps_f
