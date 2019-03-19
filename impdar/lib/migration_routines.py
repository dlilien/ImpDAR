#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Migration routines for ImpDAR

These are mostly based on older scripts in SeisUnix:
https://github.com/JohnWStockwellJr/SeisUnix/wiki
Options are:
    Kirchhoff (diffraction summation)
    Stolt (frequency wavenumber, constant velocity)
    Gazdag (phase shift, either constant or depth-varying velocity)


Author:
Benjamin Hills
benjaminhhills@gmail.com
University of Washington
Earth and Space Sciences

Mar 12 2019

"""

import numpy as np
import time

# -----------------------------------------------------------------------------

def migrationKirchhoff(dat,vel=1.69e8,vel_fn=None,nearfield=False):
    """

    Kirchhoff Migration (Berkhout 1980; Schneider 1978; Berryhill 1979)

    This migration method uses an integral solution to the scalar wave equation Yilmaz (2001) eqn 4.5.
    The algorithm cycles through all every sample in each trace, creating a hypothetical diffraciton
    hyperbola for that location,
        t(x)^2 = t(0)^2 + (2x/v)^2
    To migrate, we integrate the power along that hyperbola and assign the solution to the apex point.
    There are two terms in the integral solution, Yilmaz (2001) eqn 4.5, a far-field term and a
    near-field term. Most algorithms ignor the near-field term because it is small.

    Parameters
    ---------
    dat: data as a dictionary in the ImpDAR format
    vel: wave velocity, default is for ice
    nearfield: boolean to indicate whether or not to use the nearfield term in summation

    Output
    ---------
    dat.data: migrated data

    """

    print('Kirchhoff Migration (diffraction summation) of %.0fx%.0f matrix'%(dat.tnum,dat.snum))
    # check that the arrays are compatible
    if np.size(dat.data,1) != dat.tnum or np.size(dat.data,0) != dat.snum:
        raise ValueError('The input array must be of size (tnum,snum)')
    # start the timer
    start = time.time()
    # Calculate the time derivative of the input data
    gradD = np.gradient(dat.data,dat.travel_time/1e6,axis=0)
    # Create an empty array to fill with migrated data
    migdata = np.zeros_like(dat.data)
    # Loop through all traces
    for xi in range(dat.tnum):
        print('Migrating trace number:',xi)
        # get the trace distance
        x = dat.dist[xi]
        # Loop through all samples
        for ti in range(dat.snum):
            # get the sample time
            t = dat.travel_time[ti]/1e6
            # convert to depth
            z = vel*t/2.
            # get the radial distances between input point and output point
            rs = np.sqrt((dat.dist-x)**2.+z**2.)
            # find the cosine of the angle of the tangent line, correct for obliquity factor
            with np.errstate(invalid='ignore'):
                costheta = z/rs
            # get the exact indices from the array (closest to rs)
            Didx = [np.argmin(abs(dat.travel_time/1e6-2.*r/vel)) for r in rs]
            # integrate the farfield term
            gradDhyp = np.array([gradD[Didx[i],i] for i in range(len(Didx))])
            gradDhyp[2.*rs/vel>max(dat.travel_time/1e6)] = 0.    # zero points that are outside of the domain
            integral = np.nansum(gradDhyp*costheta/vel)  # TODO: Yilmaz eqn 4.5 has an extra r in this weight factor???
            # integrate the nearfield term
            if nearfield == True:
                Dhyp = np.array([dat.data[Didx[i],i] for i in range(len(Didx))])
                Dhyp[2.*rs/vel>max(dat.travel_time/1e6)] = 0.    # zero points that are outside of the domain
                integral += np.nansum(Dhyp*costheta/rs**2.)
            # sum the integrals and output
            migdata[ti,xi] = 1./(2.*np.pi)*integral
    dat.data = migdata.copy()
    # print the total time
    print('Kirchhoff Migration of %.0fx%.0f matrix complete in %.2f seconds'
          %(dat.tnum,dat.snum,time.time()-start))
    return dat

# -----------------------------------------------------------------------------

def migrationStolt(dat,vel=1.68e8,vel_fn=None,nearfield=False):
    """

    Stolt Migration (Stolt, 1978, Geophysics)

    This is by far the fastest migration method. It is a simple transformation from
    frequency-wavenumber (FKx) to wavenumber-wavenumber (KzKx) space.

    Parameters
    ---------
    dat: data as a dictionary in the ImpDAR format
    vel: wave velocity, default is for ice

    Output
    ---------
    dat.data: migrated data

    """

    print('Stolt Migration (f-k migration) of %.0fx%.0f matrix'%(dat.tnum,dat.snum))
    # check that the arrays are compatible
    if np.size(dat.data,1) != dat.tnum or np.size(dat.data,0) != dat.snum:
        raise ValueError('The input array must be of size (tnum,snum)')
    # save the start time
    start = time.time()
    # pad the array with zeros up to the next power of 2 for discrete fft
    nt = 2**(np.ceil(np.log(dat.snum)/np.log(2))).astype(int)
    nx = 2**(np.ceil(np.log(dat.tnum)/np.log(2))).astype(int)
    # 2D Forward Fourier Transform to get data in frequency-wavenumber space, FK = D(kx,z=0,ws)
    FK = np.fft.fft2(dat.data,(nt,nx))
    # get the temporal frequencies
    ws = 2.*np.pi*np.fft.fftfreq(nt, d=dat.dt)
    # get the horizontal wavenumbers
    kx = 2.*np.pi*np.fft.fftfreq(nx, d=np.mean(dat.trace_int))
    # interpolate from frequency (ws) into wavenumber (kz)
    from scipy.interpolate import interp2d
    interp_real = interp2d(kx,ws,FK.real)
    interp_imag = interp2d(kx,ws,FK.imag)
    # interpolation will move from frequency-wavenumber to wavenumber-wavenumber, KK = D(kx,kz,t=0)
    KK = np.zeros_like(FK)
    print('Interpolating from temporal frequency (ws) to vertical wavenumber (kz):')
    # for all temporal frequencies
    for zj in range(nt):
        kzj = ws[zj]*2./vel
        if zj%100 == 0:
            print('Interpolating',int(ws[zj]/1e6/2/np.pi),'MHz')
        # for all horizontal wavenumbers
        for xi in range(len(kx)):
            kxi = kx[xi]
            # migration conversion to wavenumber (Yilmaz equation C.53)
            wsj = vel/2.*np.sqrt(kzj**2.+kxi**2.)
            # get the interpolated FFT values, real and imaginary, S(kx,kz,t=0)
            KK[zj,xi] = interp_real(kxi,wsj) + 1j*interp_imag(kxi,wsj)
    # all vertical wavenumbers
    kz = ws*2./vel
    # grid wavenumbers for scaling calculation
    kX,kZ = np.meshgrid(kx,kz)
    # scaling for obliquity factor (Yilmaz equation C.56)
    with np.errstate(invalid='ignore'):
        scaling = kZ/np.sqrt(kX**2.+kZ**2.)
    KK *= scaling
    # the DC frequency should be 0.
    KK[0,0] = 0.+0j
    # 2D Inverse Fourier Transform to get back to distance spce, D(x,z,t=0)
    dat.data = np.real(np.fft.ifft2(KK))
    # Cut array to input matrix dimensions
    dat.data = dat.data[:dat.snum,:dat.tnum]
    # print the total time
    print('Stolt Migration of %.0fx%.0f matrix complete in %.2f seconds'
          %(dat.tnum,dat.snum,time.time()-start))
    return dat

# -----------------------------------------------------------------------------

def migrationPhaseShift(dat,vel=1.69e8,vel_fn=None,nearfield=False):
    """

    Phase-Shift Migration
    case with velocity, v=constant (Gazdag 1978, Geophysics)
    case with layered velocities, v(z) (Gazdag 1978, Geophysics)
    case with vertical and lateral velocity variations, v(x,z) (Ristow and Ruhl 1994, Geophysics)

    Phase-shifting migration for constant, layered, and laterally varying velocity structures.
    This method works down from the surface, using the imaging principle
    (i.e. summing over all frequencies to get the solution at t=0)
    for the migrated section at each step.

    Parameters
    ---------
    dat: data as a dictionary in the ImpDAR format
    vel: v(x,z)
        Up to 2-D array with three columns for velocities (m/s), and z/x (m).
        Array structure is velocities in first column, z location in second, x location in third.
        If uniform velocity (i.e. vel=constant) input constant
        If layered velocity (i.e. vel=v(z)) input array with shape (#vel-points, 2) (i.e. no x-values)
    vel_fn: filename for layered velocity input

    Output
    ---------
    dat.data: migrated data

    **
    The foundation of this script was taken from two places:
    Matlab code written by Andreas Tzanis, Dept. of Geophysics, University of Athens (2005)
    Seis Unix script sumigffd.c, Credits: CWP Baoniu Han, July 21th, 1997
    **

    """

    print('Phase-Shift Migration of %.0fx%.0f matrix'%(dat.tnum,dat.snum))
    # check that the arrays are compatible
    if np.size(dat.data,1) != dat.tnum or np.size(dat.data,0) != dat.snum:
        raise ValueError('The input array must be of size (tnum,snum)')
    # save the start time
    start = time.time()
    # pad the array with zeros up to the next power of 2 for discrete fft
    nt = 2**(np.ceil(np.log(dat.snum)/np.log(2))).astype(int)
    # get frequencies and wavenumbers
    kx = 2.*np.pi*np.fft.fftfreq(dat.tnum,d=np.mean(dat.trace_int))
    kx = np.fft.fftshift(kx)
    ws = 2.*np.pi*np.fft.fftfreq(nt,d=dat.dt)
    ws = np.fft.fftshift(ws)
    # 2D Forward Fourier Transform to get data in frequency-wavenumber space, FK = D(kx,z=0,ws)
    FK = np.fft.fft2(dat.data,(nt,dat.tnum))
    FK = np.fft.fftshift(FK)
    # Velocity structure from input
    if vel_fn is not None:
        try:
            vel = np.genfromtxt(vel_fn)
            print('Velocities loaded from %s.'%vel_fn)
        except:
            raise TypeError('File %s was given for input velocity array, but cannot be loaded. Please reformat.'%vel_fn)
    vmig = getVelocityProfile(dat,vel)
    # Migration by phase shift, frequency-wavenumber (FKx) to time-wavenumber (TKx)
    TK = phaseShift(dat, vmig, vel, kx, ws, FK)
    # Transform from time-wavenumber (TKx) to time-space (TX) domain to get migrated section
    dat.data = np.fft.ifft(np.fft.ifftshift(TK,1)).real
    # print the total time
    print('Phase-Shift Migration of %.0fx%.0f matrix complete in %.2f seconds'
          %(dat.tnum,dat.snum,time.time()-start))
    return dat

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

### Supporting functions

def phaseShift(dat, vmig, vels_in, kx, ws, FK):
    """

    Phase-Shift migration to get from frequency-wavenumber (FKx) space to time-wavenumber (TKx) space.
    This is for either constant or layered velocity v(z).

    Parameters
    ---------
    dat: data as a dictionary in the ImpDAR format
    vmig: migration velocity (m/s)
        can be constant or 1-D or 2-D array
    vels_in: v(x,z)
        Up to 2-D array with three columns for velocities (m/s), and z/x (m).
        Array structure is velocities in first column, z location in second, x location in third.
        If uniform velocity (i.e. vel=constant) input constant
        If layered velocity (i.e. vel=v(z)) input array with shape (#vel-points, 2) (i.e. no x-values)
    kx: horizontal wavenumbers
    ws: temporal frequencies
    FK: 2-D array of the data image in frequency-wavenumber space (FKx)

    Output
    ---------
    TK: 2-D array of the migrated data image in time-wavenumber space (TKx).

    **
    The foundation of this script was taken from Matlab code written by Andreas Tzanis,
    Dept. of Geophysics, University of Athens (2005)
    **

    """

    # initialize the time-wavenumber array to be filled with complex values
    TK = np.zeros((dat.snum,len(kx)))+0j

    # Uniform velocity case, vmig=constant
    if not hasattr(vmig,"__len__"):
        print('Constant velocity %s m/usec'%(vmig/1e6))
        # iterate through all frequencies
        for iw in range(len(ws)):
            w = ws[iw]
            if w == 0.0:
                w = 1e-10/dat.dt
            if iw%100 == 0:
                print('Frequency',int(w/1e6/(2.*np.pi)),'MHz')
            # remove frequencies outside of the domain
            vkx2 = (vmig*kx/2.)**2.
            ik = np.argwhere(vkx2 < w**2.)
            FFK = FK[iw,ik]
            # get the phase for shift
            phase = (-w*dat.dt*np.sqrt(1.0 - vkx2[ik]/w**2.)).real
            cc = np.conj(np.cos(phase)+1j*np.sin(phase))
            # Accumulate output image (time-wavenumber space) summed over all frequencies
            for itau in range(dat.snum):
                 FFK *= cc
                 TK[itau,ik] += FFK

    else:
        # Layered velocity case, vmig=v(z)
        if len(np.shape(vmig)) == 1:
            if len(vmig) != dat.snum:
                raise ValueError('Interpolated velocity profile is not the length of the number of samples in a trace.')
            print('1-D velocity structure, Gazdag Migration')
            print(np.shape(vels_in)[0]-1,'layers with velocities',
                    ' '.join('%.2e'%v for v in vels_in[:,0]),'(m/s), at depths',
                    ' '.join('%.1f'%t for t in vels_in[:,1]),'(m).')
            # iterate through all output travel times
            for itau in range(dat.snum):
                tau = dat.travel_time[itau]/1e6
                if itau%100 == 0:
                    print('Time %.2e of %.2e' %(tau,dat.travel_time[-1]/1e6))
                # iterate through all frequencies
                for iw in range(len(ws)):
                    w = ws[iw]
                    if w == 0.0:
                        w = 1e-10/dat.dt
                    # cosine squared
                    coss = 1.0+0j - (0.5 * vmig[itau]*kx/w)**2.
                    # calculate phase for shift
                    phase = (-w*dat.dt*np.sqrt(coss)).real
                    cc = np.conj(np.cos(phase)+1j*np.sin(phase))
                    FK[iw] *= cc
                    # zero if outside domain
                    idx = coss <= (tau/dat.travel_time[-1]/1e6)**2.
                    FK[iw,idx] = complex(0.0,0.0)
                    # sum over all frequencies
                    TK[itau] += FK[iw]

        # Layered velocity case, vmig=v(x,z)
        elif np.shape(vmig)[1] == 3:
            print('2-D velocity structure, Fourier Finite-Difference Migration')
            # iterate through all output travel times
            for itau in range(dat.snum):
                tau = dat.travel_time[itau]/1e6
                if itau%100 == 0:
                    print('Time %.2e of %.2e' %(tau,dat.travel_time[-1]/1e6))
                # iterate through all frequencies
                for iw in range(len(ws)):
                    w = ws[iw]
                    if w == 0.0:
                        w = 1e-10/dat.dt

                    ###
                    # add FD operator
                    # add thin lens term
                    ###

                    # cosine squared
                    coss = 1.0+0j - (0.5 * vmig[itau]*kx/w)**2.
                    # calculate phase for shift
                    phase = (-w*dat.dt*np.sqrt(coss)).real
                    cc = np.conj(np.cos(phase)+1j*np.sin(phase))
                    FK[iw] *= cc
                    # zero if outside domain
                    idx = coss <= (tau/dat.travel_time[-1]/1e6)**2.
                    FK[iw,idx] = complex(0.0,0.0)
                    # sum over all frequencies
                    TK[itau] += FK[iw]

    # Cut to original array size
    TK = TK[:,:dat.tnum]
    # Normalize for inverse FFT
    TK /= dat.snum
    return TK

# -----------------------------------------------------------------------------

def finiteDiff(dat, vel, vel_fn, kx, ws, FK):
    """

    Fourier Finite-Difference migration to get from frequency-wavenumber (FKx) space to time-wavenumber (TKx) space.
    This is for either constant or variable velocity v(x,z).

    Parameters
    ---------
    dat: data as a dictionary in the ImpDAR format
    vmig: migration velocity (m/s)
    vel: 3-column array of wave velocities (m/s) and x/z locations (m)
        If only one layer (i.e. constant velocity) input one layer with zeroes for x/z.
    vel_fn: filename for layered velocity input
    kx: horizontal wavenumbers
    ws: temporal frequencies
    FK: 2-D array of the data image in frequency-wavenumber space (FKx)

    Output
    ---------
    TK: 2-D array of the migrated data image in time-wavenumber space (TKx).

    **
    The foundation of this script was taken from Seis Unix script sumigffd.c
    Credits: CWP Baoniu Han, July 21th, 1997
    **

    """

    return

# -----------------------------------------------------------------------------

def getVelocityProfile(dat,vels_in):
    """

    Map the layered velocity structure into the shape of the data.

    Parameters
    ---------
    dat: data as a dictionary in the ImpDAR format
    vels_in: v(x,z)
        Up to 2-D array with three columns for velocities (m/s), and z/x (m).
        Array structure is velocities in first column, z location in second, x location in third.
        If uniform velocity (i.e. vel=constant) input constant
        If layered velocity (i.e. vel=v(z)) input array with shape (#vel-points, 2) (i.e. no x-values)

    Output
    ---------
    vmig: 2-D array of migration velocities (m/s), shape is (#traces, #samples).
        If constant input velocity, output is constant.
        If only z-component in input velocity array, output is v(z)

    """

    if not hasattr(vels_in,"__len__"):
        return vels_in

    nlay,dimension = np.shape(vels_in)
    vel_v = vels_in[:,0]
    vel_z = vels_in[:,1]

    ### Layered Velocity
    if nlay > 1 and dimension == 2:
        # Depth of layer boundaries from input velocity locations
        zs = np.max(vel_v)*dat.travel_time/1e6/2.    # depth array for maximum possible penetration
        if vel_z[0] > np.nanmin(zs):
            vel_z[0] = np.nanmin(zs)
        if vel_z[-1] < np.nanmax(zs):
            vel_z[-1] = np.nanmax(zs)
        # Compute times from input velocity/location array (vels_in)
        vel_t = 2.*vel_z/vel_v
        # Interpolate to get t(z) for maximum penetration depth array
        from scipy.interpolate import interp1d
        tinterp = interp1d(vel_z,vel_t)
        tofz = tinterp(zs)
        # Compute z(t) from monotonically increasing t
        zinterp = interp1d(tofz,zs)
        zoft = zinterp(dat.travel_time/1e6)
        # TODO: does this need more rigorous testing? Will the interpolated range always be big enough?
        if dat.travel_time[-1]/1e6 > tofz[-1]:
            raise ValueError('Two-way travel time array extends outside of interpolation range')
        # Compute vmig(t) from z(t)
        vmig = 2.*np.gradient(zoft,dat.dt)

    ### Lateral Velocity Variations
    elif nlay > 1 and dimension == 3:
        vel_x = vels_in[:,2]
        print(vel_x)
        vmig = vel_v

    return vmig
