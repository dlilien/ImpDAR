#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Empirical attenuation calculations from Hills et al. (2020)

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 26 2019

"""

import numpy as np
from scipy import stats

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
### Single-Reflector Methods
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

def attenuation_method2(dat,picknum,sigPc=0.,sigZ=0.,Cint=.95,u=1.69e8,*args, **kwargs):
    """
    ### Method 2 from the attenuation framework (Hills et al., 2020) ###
    Method 1 is the same but only for the bed reflector, as is most common.

    Based on Jacobel et al. (2009)
    This method fits a line to the measured power for an individual reflector (in log space)
    The resulting attenuation rate represents a depth-averaged value
    over the depth range which the layer spans.

    Parameters
    ----------
    picknum: int
        pick number to do the attenuation calculation on
    sigPc:  float; optional
        standard deviation in measured power (error) used to constrain the regression
        defaults to 0 (i.e. simple regression)
    sigZ:   float; optional
        standard deviation in measured depth (error) used to constrain the regression
        defaults to 0 (i.e. simple regression)
    Cint:   float; optional
        confidence interval with which to describe the resulting attenuation error
        default 95%
    u:      float; optional
        light velocity in ice

    Output
    ----------
    N:    float
        One-way attenuation rate (dB/km)
    Nerr: float
        Error in one-way attenuation rate (dB/km)
    """

    # get a pick depth
    if 'z' in vars(dat.picks):
        Z = dat.picks.z
    else:
        print('Warning: setting pick depth for constant velocity in ice.')
        Z = dat.picks.time*u/2/1e6

    # Get pick from index and remove nan values
    Pc = 10.*np.log10(dat.picks.corrected_power[picknum])
    Z = Z[picknum]
    idx = ~np.isnan(Pc) & ~np.isnan(Z)
    Pc = Pc[idx]
    Z = Z[idx]

    # Convert to km
    if np.any(Z > 10.):
        Z/=1000.
    if sigZ > .1:
        sigZ/=1000.

    Szz = np.sum((Z-np.mean(Z))**2.)
    Spp = np.sum((Pc-np.mean(Pc))**2.)
    Szp = np.sum((Z-np.mean(Z))*(Pc-np.mean(Pc)))

    if sigZ == 0 and sigPc == 0:
        # Simple regression
        N = -(Szp)/Szz
        alpha = np.mean(Pc) + N*np.mean(Z)
        # Error based on vertical distance from line only
        Pc_err = np.sum((Pc - ((-N)*Z + alpha))**2.)
        sigN = np.sqrt(Pc_err/Szz/(len(Z)-2))
        tscore = stats.t.ppf(1.-(1.-Cint)/2., len(Z)-2)
        Nerr = tscore*sigN
    else:
        # Deming regression after Casella and Berger (2002) section 12.2
        lam = (sigZ**2.)/(sigPc**2.)
        # Regression slope, eq. 12.2.16
        N = -(-Szz+lam*Spp+np.sqrt((Szz-lam*Spp)**2.+4.*lam*Szp**2.))/(2.*lam*Szp)
        alpha = np.mean(Pc) + N*np.mean(Z)
        # Standard deviation in slope 12.2.22
        sigN = np.sqrt(((1.+lam*N**2.)**2.*(Szz*Spp-Szp**2.))/((Szz-lam*Spp)**2.+4.*lam*Szp**2.))
        tscore = stats.t.ppf(1.-(1.-Cint)/2., len(Z)-2)
        # Error using Gleser's Modification with 95% confidence interval
        Nerr = tscore*sigN/(np.sqrt(len(Z)-2))

    # Final Output as a one-way rate in dB/km
    N *= 1/2.
    Nerr *= 1/2.

    return N,Nerr

# -----------------------------------------------------------------------------------------------------

def attenuation_method3(dat,picknum,Ns=np.arange(30.),Nh_target=1.,Cw=0.1,win_init=100,win_step=100,u=1.69e8):
    """
    ### Method 3 from the attenuation framework (Hills et al., 2020) ###

    Based on Schroeder et al. (2016a, 2016b)
    This method decorrolates the attenuation-corrected power from the ice thickness
    Assumes constant reflectivity

    Parameters
    ----------
    picknum: int
        pick number to do the attenuation calculation on
    Ns:     array; optional
        Attenuation rates to test (one-way in dB/km)
    Nh_target:  float; optional
        Radiometric resolution target
    Cw:         float; optional
        Minimum correlation coefficient threshold
    win_init:   int; optional
        Initial number of traces for window size
    win_step:   int; optional
        Number of traces to increase the window size at each step
    u:          float; optional
        light velocity in ice

    Output
    ----------
    N_result:  array
        One-way attenuation rate (dB/km)
    win_result:   array
        resulting window size (number of traces)
    """

    # get a pick depth
    if 'z' in vars(dat.picks):
        Z = dat.picks.z
    else:
        print('Warning: setting pick depth for constant velocity in ice.')
        Z = dat.picks.time*u/2/1e6

    # Get pick from index and remove nan values
    Pc = 10*np.log10(dat.picks.corrected_power[picknum])
    Z = Z[picknum]
    idx = ~np.isnan(Pc) & ~np.isnan(Z)
    Pc = Pc[idx]
    Z = Z[idx]

    # Convert to km
    if np.any(Z > 10.):
        Z/=1000.

    # Create empty arrays to fill for the resulting attenuation rate and window size
    N_result = np.zeros((dat.tnum,))
    win_result = np.zeros((dat.tnum,))
    C = np.zeros_like(Ns)
    # Loop through all the traces
    for tr in range(win_init//2,dat.tnum-win_init//2):
        # zero out the correlation coefficient array
        C[:] = 0.
        # Initial window size
        win = win_init
        # Radiometric Resolution (needs to converge onto Nh_target before the attenuation rate is accepted)
        Nh = Nh_target + 1.
        # while radiometric resolution is outside target range and window is fully within the profile
        while Nh > Nh_target and win//2<=tr and win//2<=(len(Z)-tr):
            # thickness and power in the window
            z = Z[tr-win//2:tr+win//2]
            pc = Pc[tr-win//2:tr+win//2]
            # loop through all the possible attenuation rates
            sum2 = np.sqrt(sum((z-np.mean(z))**2.))
            # TODO: I think I could substantially speed things up in here
            for j,Nj in enumerate(Ns):
                # attenuation-corrected power, Schroeder et al. (2016) eq. 4
                pa = pc + 2.*z*Nj
                # calculate the correlation coefficient, Schroeder et al. (2016) eq. 5
                sum1 = sum((z-np.mean(z))*(pa-np.mean(pa)))
                sum3 = np.sqrt(sum((pa-np.mean(pa))**2.))
                C[j] = abs(sum1/(sum2*sum3))
            # Whichever value has the lowest correlation coefficient is chosen
            Cm = np.min(C)
            Nm = Ns[C==Cm]
            C0 = C[Ns==0]
            # If the minimum correlation coefficient is below threshold, Cw,
            # and the zero correlation coefficient is above
            # then update the radiometric resolution
            if Cm < Cw and C0 > Cw:
                Nh = np.max(Ns[C<Cw])-np.min(Ns[C<Cw])
            win += win_step

        N_result[tr] = Nm
        win_result[tr] = win

    return N_result,win_result

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
### Multiple-Reflector Methods
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

def attenuation_method5(dat,picknums,win=1,sigPc=0,sigZ=0,Cint=.95,u=1.69e8,*args,**kwargs):
    """
    ### Method 5 from the attenuation framework (Hills et al., 2020) ###

    Based on MacGregor et al. (2014) and Matsuoka et al. (2010)
    This method fits a line to the measured power for internal reflectors (in log space)

    Parameters
    ----------
    picknums: array
        pick numbers to do the attenuation calculation on
    win:    int; optional
        number of traces to use in the regression (default: 1)
    sigPc:  float; optional
        standard deviation in measured power (error) used to constrain the regression
        defaults to 0 (i.e. simple regression)
    sigZ:   float; optional
        standard deviation in measured depth (error) used to constrain the regression
        defaults to 0 (i.e. simple regression)
    Cint:   float; optional
        confidence interval with which to describe the resulting attenuation error
        default 95%
    u:      float; optional
        light velocity in ice

    Output
    ----------
    N_result:  array
        horizontal attenuation profile (as a one-way rate in dB/km)
    Nerr_result: array
        Error in one-way attenuation rate (dB/km)
    """

    # get a pick depth
    if 'z' in vars(dat.picks):
        Z = dat.picks.z
    else:
        print('Warning: setting pick depth for constant velocity in ice.')
        Z = dat.picks.time*u/2/1e6

    # Convert to km
    if np.any(Z > 10.):
        Z/=1000.
    if sigZ > .1:
        sigZ/=1000.

    # create empty arrays for output
    N_result = np.nan*np.empty((dat.tnum,))
    Nerr_result = np.nan*np.empty((dat.tnum,))
    a_result = np.nan*np.empty((dat.tnum,))
    # calculate the attenuation rate for each desired trace (or window)
    for tr in np.arange(win//2,dat.tnum-win//2):
        # grab the data within the window
        # Get pick from index and remove nan values
        pc = np.squeeze(10.*np.log10(dat.picks.corrected_power[picknums,tr-win//2:tr+win//2+1]))
        z = np.squeeze(Z[picknums,tr-win//2:tr+win//2+1])
        idx = ~np.isnan(pc) & ~np.isnan(z)
        pc = pc[idx]
        z = z[idx]

        # If there are not enough picked layers for this window, output nan
        if len(pc)<5:
            N_result[tr] = np.nan
            Nerr_result[tr] = np.nan
            a_result[tr] = np.nan
        else:
            Szz = np.sum((z-np.mean(z))**2.)
            Spp = np.sum((pc-np.mean(pc))**2.)
            Szp = np.sum((z-np.mean(z))*(pc-np.mean(pc)))

            if sigZ == 0 and sigPc == 0:
                # Simple regression
                N = -(Szp)/Szz
                alpha = np.mean(pc) + N*np.mean(z)
                # Error based on vertical distance from line only
                pc_err = np.sum((pc - ((-N)*z + alpha))**2.)
                sigN = np.sqrt(pc_err/Szz/(len(z)-2))
                tscore = stats.t.ppf(1.-(1.-Cint)/2., len(z)-2)
                Nerr = tscore*sigN
            else:
                # Deming regression after Casella and Berger (2002) section 12.2
                lam = (sigZ**2.)/(sigPc**2.)
                # Regression slope, eq. 12.2.16
                N = -(-Szz+lam*Spp+np.sqrt((Szz-lam*Spp)**2.+4.*lam*Szp**2.))/(2.*lam*Szp)
                alpha = np.mean(pc) + N*np.mean(z)
                # Standard deviation in slope 12.2.22
                sigN = np.sqrt(((1.+lam*N**2.)**2.*(Szz*Spp-Szp**2.))/((Szz-lam*Spp)**2.+4.*lam*Szp**2.))
                tscore = stats.t.ppf(1.-(1.-Cint)/2., len(z)-2)
                # Error using Gleser's Modification with 95% confidence interval
                Nerr = tscore*sigN/(np.sqrt(len(z)-2))

            # Final Output as a one-way rate in dB/km
            N_result[tr] = N*.5 #one-way attenuation rate
            Nerr_result[tr] = Nerr*.5

    return N_result,Nerr_result

# -----------------------------------------------------------------------------------------------------

def attenuation_method6a(dat,picknums,att_ds,win=500.,sigPc=0,sigZ=0,Cint=.95,u=1.69e8,*args,**kwargs):
    """
    ### Method 6 from the attenuation framework (Hills et al., 2020) ###

    This method groups all the picks from all traces together.
    A window of fixed size then moves from top to bottom of the profile,
    fitting the regression to all data within the window for each step.

    Parameters
    ----------
    picknums:   array
        picks to include in the calculation
    att_ds: array
        depths at which to center attenuation calculations
    win:    float
        window over which to calculate the attenuation rate
    sigPc:  float; optional
        standard deviation in measured power (error) used to constrain the regression
        defaults to 0 (i.e. simple regression)
    sigZ:   float; optional
        standard deviation in measured depth (error) used to constrain the regression
        defaults to 0 (i.e. simple regression)
    Cint:   float; optional
        confidence interval with which to describe the resulting attenuation error
        default 95%
    u:      float; optional
        light velocity in ice

    Output
    ---------
    N_result:      1-d array for attenuation rate in dB/km
    Nerr_result:      1-d array for attenuation rate error in dB/km
    """

    # get a pick depth
    if 'z' in vars(dat.picks):
        Z = dat.picks.z
    else:
        print('Warning: setting pick depth for constant velocity in ice.')
        Z = dat.picks.time*u/2/1e6

    # Get pick from index and remove nan values
    Pc = 10.*np.log10(dat.picks.corrected_power[picknums].flatten())
    Z = Z[picknums].flatten()
    idx = ~np.isnan(Pc) & ~np.isnan(Z)
    Pc = Pc[idx]
    Z = Z[idx]

    # Convert to km
    if np.any(Z > 10.):
        Z/=1000.
    if np.any(att_ds>10.):
        att_ds/=1000.
    if win>10.:
        win/=1000.

    # Create empty arrays to fill for the output attenuation rate and window size
    N_result = np.zeros_like(att_ds).astype(float)
    Nerr_result = np.zeros_like(att_ds).astype(float)

    # loop through all the depths
    for i,att_d in enumerate(att_ds):
        z = Z[np.logical_and(Z>(att_d-win/2),Z<(att_d+win/2))]
        pc = Pc[np.logical_and(Z>(att_d-win/2),Z<(att_d+win/2))]
        if len(z)<5:
            N_result[i] = np.nan
            Nerr_result[i] = np.nan
            continue

        # Sum of squares
        Szz = np.sum((z-np.mean(z))**2.)
        Spp = np.sum((pc-np.mean(pc))**2.)
        Szp = np.sum((z-np.mean(z))*(pc-np.mean(pc)))

        if sigZ == 0 and sigPc == 0:
            # Simple regression
            N = -(Szp)/Szz
            alpha = np.mean(pc) + N*np.mean(z)
            # Error based on vertical distance from line only
            pc_err = np.sum((pc - ((-N)*z + alpha))**2.)
            sigN = np.sqrt(pc_err/Szz/(len(z)-2))
            tscore = stats.t.ppf(1.-(1.-Cint)/2., len(z)-2)
            Nerr = tscore*sigN
        else:
            # Deming regression after Casella and Berger (2002) section 12.2
            lam = (sigZ**2.)/(sigPc**2.)
            # Regression slope, eq. 12.2.16
            N = -(-Szz+lam*Spp+np.sqrt((Szz-lam*Spp)**2.+4.*lam*Szp**2.))/(2.*lam*Szp)
            alpha = np.mean(pc) + N*np.mean(z)
            # Standard deviation in slope 12.2.22
            sigN = np.sqrt(((1.+lam*N**2.)**2.*(Szz*Spp-Szp**2.))/((Szz-lam*Spp)**2.+4.*lam*Szp**2.))
            tscore = stats.t.ppf(1.-(1.-Cint)/2., len(z)-2)
            # Error using Gleser's Modification with 95% confidence interval
            Nerr = tscore*sigN/(np.sqrt(len(z)-2))

        # Fill in the result array
        N_result[i] = .5*N
        Nerr_result[i] = .5*Nerr

    return N_result, Nerr_result

# -----------------------------------------------------------------------------------------------------

def attenuation_method6b(dat,picknums,att_ds,Ns=np.arange(30.),Nh_target=1.,Cw=0.1,win_init=100.,win_step=100.,u=1.69e8,*args,**kwargs):
    """
    ### Method 6b from the attenuation framework (Hills et al., 2020) ###

    Based on Schroeder et al. (2016a, 2016b)
    This method minimizes the correlation coeffiecient between attenuation
    rate and ice thickness

    Here, we are using the Schroeder method (i.e. from method 3) but in the vertical.
    All picks from all traces are grouped together and the window moves from
    top to bottom of the ice column, adjusting size to optimize the fit.

    Parameters
    ----------
    picknums: int
        pick number to do the attenuation calculation on
    att_ds
    Ns:     array; optional
        Attenuation rates to test (one-way in dB/km)
    Nh_target:  float; optional
        Radiometric resolution target
    Cw:         float; optional
        Minimum correlation coefficient threshold
    win_init:   float; optional
        Initial number of traces for window size
    win_step:   float; optional
        Number of traces to increase the window size at each step
    u:      float; optional
        light velocity in ice

    Output
    ----------
    N_result:  array
        One-way attenuation rate (dB/km)
    win_result:   array
        resulting window size (m)
    """

    # get a pick depth
    if 'z' in vars(dat.picks):
        Z = dat.picks.z
    else:
        print('Warning: setting pick depth for constant velocity in ice.')
        Z = dat.picks.time*u/2/1e6

    # Get pick from index and remove nan values
    Pc = 10.*np.log10(dat.picks.corrected_power[picknums].flatten())
    Z = Z[picknums].flatten()
    idx = ~np.isnan(Pc) & ~np.isnan(Z)
    Pc = Pc[idx]
    Z = Z[idx]

    # Convert to km
    if np.any(Z > 10.):
        Z/=1000.
    if np.any(att_ds>10.):
        att_ds/=1000.
    if win_init>10.:
        win_init/=1000.
        win_step/=1000.

    # Create empty arrays to fill for the output attenuation rate and window size
    N_result = np.zeros_like(att_ds)
    win_result = np.zeros_like(att_ds)
    C = np.zeros_like(Ns)
    # loop through all the depths
    for i,att_d in enumerate(att_ds):
        # current depth for attenuation calc
        att_d = att_ds[i]
        # Correlation Coefficient (starts empty)
        C[:] = np.nan
        # Initial window size
        win = win_init
        # Radiometric Resolution (needs to converge to Nh_target)
        Nh = Nh_target + 1.
        while Nh > Nh_target and att_d-win/2>=np.nanmin(Z) and att_d+win/2<=np.nanmax(Z):
            # thickness and power in the window
            z = Z[np.argwhere(abs(Z-att_d)<win/2)]
            pc = Pc[np.argwhere(abs(Z-att_d)<win/2)]
            # loop through all the possible attenuation rates
            sum2 = np.sqrt(np.nansum((z-np.nanmean(z))**2.))
            # TODO: I think I could substantially speed things up in here
            for n in range(len(Ns)):
                # attenuation-corrected power, Schroeder et al. (2016) eq. 4
                pa = pc + 2.*z*Ns[n]
                # calculate the correlation coefficient, Schroeder et al. (2016) eq. 5
                sum1 = np.nansum((z-np.nanmean(z))*(pa-np.nanmean(pa)))
                sum3 = np.sqrt(np.nansum((pa-np.nanmean(pa))**2.))
                if np.any(np.isnan([sum1,sum2,sum3])):
                    C[n] = np.nan
                else:
                    C[n] = abs(sum1/(sum2*sum3))
            # Whichever value has the lowest correlation coefficient is chosen
            Cm = np.nanmin(C)
            Nm = Ns[C==Cm]
            C0 = C[Ns==0]
            if Cm < Cw and C0 > Cw:
                Nh = (np.max(Ns[C<Cw])-np.min(Ns[C<Cw]))/2.
            # get ready for the next iteration
            win += win_step
        # output
        N_result[i] = Nm
        win_result[i] = win*1000.
    return N_result, win_result

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
### Secondary Reflection Methods
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

def attenuation_method7(dat,primary_picknum,secondary_picknum,Rib=-.22,Rfa=-17,u=1.69e8,*args,**kwargs):
    """
    ### Method 7 from the attenuation framework (Hills et al., 2020) ###

    # Based on MacGregor (2011) and Christianson et al. (2016)
    # This method fits a line to the to the measured power from the basal reflector (in log space)

    Parameters
    ----------
    primary_picknum:    int
        index for primary reflection
    secondary_picknum:  int
        index for secondary reflection
    Rib:    float; optional
        reflectivity of the ice-bed interface (in dB);
        default -0.22 for ice-seawater interface from Christianson et al. (2016; Appendix A)
    Rfa:    flaot; optional
        reflectivity of the firn-air interface (in dB);
        default -17 from Christianson et al. (2016; Appendix A)
    u:      float; optional
        light velocity in ice

    Output
    ---------
    mean(N):    flaot
        mean of the calculated attenuation rates
    std(N):     flaot
        standard deviation of the calculated attenuation rates
    """

    # get a pick depth
    if 'z' in vars(dat.picks):
        Z = dat.picks.z
    else:
        print('Warning: setting pick depth for constant velocity in ice.')
        Z = dat.picks.time*u/2/1e6

    # Convert to km
    if np.any(Z > 10.):
        Z/=1000.

    # Get pick from index and remove nan values
    P1 = dat.picks.corrected_power[primary_picknum]
    P2 = dat.picks.corrected_power[secondary_picknum]
    Z1 = Z[primary_picknum]
    Z2 = Z[secondary_picknum]
    idx = ~np.isnan(P1) & ~np.isnan(P2) & ~np.isnan(Z1) & ~np.isnan(Z2)
    P1 = P1[idx]
    P2 = P2[idx]
    Z1 = Z1[idx]
    Z2 = Z2[idx]

    # Check that the secondary depth is double the primary
    if not abs(np.nanmean(Z1)*2. - np.nanmean(Z2)) < .1*np.nanmean(Z1):
        raise ValueError('The secondary reflection is not twice as deep as the primary.')

    # convert all terms out of log space in order to use eq. A4
    Rfa = 10**(Rfa/10.)
    Rib = 10**(Rib/10.)
    # Calculate the attenuation length scale with Christianson et al. (2016) eq. A4
    La = -2.*Z1/np.log((4./(Rib*Rfa))*(P2/P1))
    # Then attenuation rate is (following Jacobel et al. (2009))
    N = 10.*np.log10(np.e)/La
    return np.nanmean(N),np.nanstd(N)
