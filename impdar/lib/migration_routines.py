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

def migrationKirchhoff(dat,vel=1.69e8,nearfield=False,**kwargs):
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
    dat.migdata: migrated data

    """

    print('Kirchhoff Migration (diffraction summation) of %.0fx%.0f matrix'%(dat.tnum,dat.snum))
    # check that the arrays are compatible
    if np.size(dat.data,1) != dat.tnum or np.size(dat.data,0) != dat.snum:
        raise ValueError('The input array must be of size (tnum,snum)')
    # start the timer
    start = time.time()
    # Calculate the time derivative of the input data
    gradD = np.gradient(dat.data,dat.travel_time,axis=0)
    # Create an empty array to fill with migrated data
    dat.migdata = np.zeros_like(dat.data)
    # Loop through all traces
    for xi in range(dat.tnum):
        print('Migrating trace number:',xi)
        # get the trace distance
        x = dat.dist[xi]
        # Loop through all samples
        for ti in range(dat.snum):
            # get the sample time
            t = dat.travel_time[ti]
            # convert to depth
            z = vel*t/2.
            # get the radial distances between input point and output point
            rs = np.sqrt((dat.dist-x)**2.+z**2.)
            # find the cosine of the angle of the tangent line, correct for obliquity factor
            costheta = z/rs
            # get the exact indices from the array (closest to rs)
            Didx = [np.argmin(abs(dat.travel_time-2.*r/vel)) for r in rs]
            # integrate the farfield term
            gradDhyp = np.array([gradD[Didx[i],i] for i in range(len(Didx))])
            gradDhyp[2.*rs/vel>max(dat.travel_time)] = 0.    # zero points that are outside of the domain
            integral = np.nansum(gradDhyp*costheta/vel)  # TODO: Yilmaz eqn 4.5 has an extra r in this weight factor???
            # integrate the nearfield term
            if nearfield == True:
                Dhyp = np.array([dat.data[Didx[i],i] for i in range(len(Didx))])
                Dhyp[2.*rs/vel>max(dat.travel_time)] = 0.    # zero points that are outside of the domain
                integral += np.nansum(Dhyp*costheta/rs**2.)
            # sum the integrals and output
            dat.migdata[ti,xi] = 1/(2.*np.pi)*integral
    # print the total time
    print('Kirchhoff Migration of %.0fx%.0f matrix complete in %.2f seconds'
          %(len(dat.dist),len(dat.travel_time),time.time()-start))
    return dat

# -----------------------------------------------------------------------------

def migrationStolt(dat,vel=1.68e8):
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
    dat.migdata: migrated data

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
    # 2D Forward Fourier Transform to get data in frequency-wavenumber space, FK = D(kx,z=0,omega)
    FK = np.fft.fft2(dat.data,(nt,nx))
    # get the temporal frequencies
    omega = np.fft.fftfreq(nt, d=dat.dt)
    # get the horizontal wavenumbers
    kx = np.fft.fftfreq(nx, d=np.mean(dat.trace_int))
    # interpolate from frequency (omega) into wavenumber (kz)
    from scipy.interpolate import interp2d
    interp_real = interp2d(kx,omega,FK.real)
    interp_imag = interp2d(kx,omega,FK.imag)
    # interpolation will move from frequency-wavenumber to wavenumber-wavenumber, KK = D(kx,kz,t=0)
    KK = np.zeros_like(FK)
    print('Interpolating from temporal frequency ($\omega$) to vertical wavenumber (kz):')
    # for all temporal frequencies
    for zj in range(len(omega)):
        kzj = omega[zj]*2./vel
        if zj%100 == 0:
            print('Interpolating',int(omega[zj]/1e6),'MHz')
        # for all horizontal wavenumbers
        for xi in range(len(kx)):
            kxi = kx[xi]
            # migration conversion to wavenumber (Yilmaz equation C.53)
            omegaj = vel/2.*np.sqrt(kzj**2.+kxi**2.)
            # get the interpolated FFT values, real and imaginary, S(kx,kz,t=0)
            KK[zj,xi] = interp_real(kxi,omegaj) + 1j*interp_imag(kxi,omegaj)
    # all vertical wavenumbers
    kz = omega*2./vel
    # grid wavenumbers for scaling calculation
    kX,kZ = np.meshgrid(kx,kz)
    # scaling for obliquity factor (Yilmaz equation C.56)
    scaling = kZ/np.sqrt(kX**2.+kZ**2.)
    KK *= scaling
    # the DC frequency should be 0.
    KK[0,0] = 0.+0j
    # 2D Inverse Fourier Transform to get back to distance spce, D(x,z,t=0)
    dat.migdata = np.real(np.fft.ifft2(KK))
    # Cut array to input matrix dimensions
    dat.migdata = dat.migdata[:dat.snum,:dat.tnum]
    # print the total time
    print('Stolt Migration of %.0fx%.0f matrix complete in %.2f seconds'
          %(dat.tnum,dat.snum,time.time()-start))
    return dat

# -----------------------------------------------------------------------------

def migrationGazdag(dat,vels_in=np.array([[1.69e8,0]])):
    """

    Gazdag Migration (Gazdag 1978, Geophysics)

    Phase-shifting migration for constant or layered velocity structures.
    This method works down from the surface, using the imaging principle
    (i.e. summing over all frequencies to get the solution at t=0)
    for the migrated section at each step.

    Parameters
    ---------
    dat: data as a dictionary in the ImpDAR format
    vels_in: 2-D array of layered wave velocities (m/s) and layer thickness (m)
        Structure is velocities in first column, layer thicknesses in second
        If only one layer (i.e. constant velocity) input one layer with zero thickness.

    Output
    ---------
    dat.migdata: migrated data


    **
    The foundation of this script was taken from Matlab code written by Andreas Tzanis,
    Dept. of Geophysics, University of Athens (2005)
    **

    """

    print('Gazdag Migration (phase-shift migration) of %.0fx%.0f matrix'%(dat.tnum,dat.snum))
    # check that the arrays are compatible
    if np.size(dat.data,1) != dat.tnum or np.size(dat.data,0) != dat.snum:
        raise ValueError('The input array must be of size (tnum,snum)')
    # save the start time
    start = time.time()
    # Velocity structure from input
    nlay, vpairs = np.shape(vels_in)
    vmig = getVelocityProfile(vels_in,nlay,dat.snum,dat.tnum,dat.dt,dat.travel_time)
    # Fourier transform to frequency-wavenumber (FKx) domain
    FK = np.fft.fft2(dat.data)
    FK = np.fft.fftshift(FK)
    # Migration by phase shift, frequency-wavenumber (FKx) to time-wavenumber (TKx)
    if nlay == 1:
        TK = phaseShiftConstantVel(dat, vmig, FK)
    elif nlay > 1:
        TK = phaseShiftLayeredVel(dat, vmig, FK)
    # Transform from time-wavenumber (TKx) to time-space (TX) domain to get migrated section
    dat.migdata = np.fft.ifft(np.fft.ifftshift(TK,1)).real
    # print the total time
    print('Gazdag Migration of %.0fx%.0f matrix complete in %.2f seconds'
          %(dat.tnum,dat.snum,time.time()-start))
    return dat.migdata




def phaseShiftConstantVel(dat, vmigc, FK):
    """

    Phase-Shift migration to get from frequency-wavenumber (FKx) space to time-wavenumber (TKx) space.
    This is in the case of a constant velocity medium. For variable velocity see phaseShiftLayeredVel().

    Parameters
    ---------
    ns: number of samples
    ntr: number of traces
    dt: time step between samples (s)
    dx: horizontal step between traces (m)
    vmigc: constant migration velocity (m/s)
    FK: 2-D array of the data image in frequency-wavenumber space (FKx)

    Output
    ---------
    TK: 2-D array of the migrated data image in time-wavenumber space (TKx).

    **
    The foundation of this script was taken from Matlab code written by Andreas Tzanis,
    Dept. of Geophysics, University of Athens (2005)
    **

    """

    # initialize the time-wavenumber array to be filled
    TK = np.zeros((dat.snum,dat.tnum))+0j
    # get frequencies, TODO: why the 2pi?
    dx = np.mean(dat.trace_int)
    kx = 2.*np.pi*np.fft.fftfreq(dat.tnum,d=dx)
    kx = np.fft.fftshift(kx)
    ws = 2.*np.pi*np.fft.fftfreq(dat.snum,d=dat.dt)
    ws = np.fft.fftshift(ws)
    for iw in range(dat.snum):
        w = ws[iw]
        if w == 0.0:
            w = 1e-10/dat.dt
        if iw%100 == 0:
            print('Frequency',int(w/1e6/2.*np.pi),'MHz')
        vkx2 = (vmigc*vmigc * kx**2.)/4.
        ik = np.argwhere(vkx2 < w**2.)
        FFK = FK[iw,ik]
        phase = (-w*dat.dt*np.sqrt(1.0 - vkx2[ik]/w**2.)).real
        cc    = np.conj(np.cos(phase)+1j*np.sin(phase))
        # Accumulate image summed over all frequencies
        for itau in range(dat.snum):
             FFK *= cc
             TK[itau,ik] += FFK
    # Normalize for inverse FFT
    TK /= dat.snum
    return TK



def phaseShiftLayeredVel(dat, vmigv, FK):
    """

    Phase-Shift migration to get from frequency-wavenumber (FKx) space to time-wavenumber (TKx) space.
    This is in the case of a layered velocity medium (i.e. horizontally homogenous, but vertically
    variable). For constant velocity see phaseShiftConstantVel().

    Parameters
    ---------
    ns: number of samples
    ntr: number of traces
    dt: time step between samples (s)
    dx: horizontal step between traces (m)
    vmigv: 1-D array of migration velocities of length=ns (m/s)
    FK: 2-D array of the data image in frequency-wavenumber space (FKx)

    Output
    ---------
    TK: 2-D array of the migrated data image in time-wavenumber space (TKx).

    **
    The foundation of this script was taken from Matlab code written by Andreas Tzanis,
    Dept. of Geophysics, University of Athens (2005)
    **

    """
    # Set up iteration constants
    dx = np.mean(dat.trace_int)
    ns = dat.snum
    dt = dat.dt
    ntr = dat.tnum

    nw, w0, dw = ns, -np.pi/dt, 2.*np.pi/(ns*dt)
    nx, kx0 = ntr, -np.pi/dx
    dkx = 2.*np.pi/(nx*dx)
    ntau, dtau = ns, dt
    ft = 0
    ftau, tmax = ft, ft+(ntau-1)*dtau
    # initialize image array
    TK = np.zeros((ns,ntr))+0j

    kx = np.arange(kx0,-kx0,dkx)
    for ikx in range(nx):
        kx = kx0 + (ikx)*dkx
        for itau in range(ntau):
            tau = ftau + (itau)*dtau
            for iw in range(nw):
                w = w0 + (iw)*dw
                if w == 0.0:
                    w = 1e-10/dt
                coss = 1.0 - (0.5 * vmigv[itau]*kx/w)**2.
                if coss > (tau/tmax)**2.:
                    phase = (-w*dt*np.sqrt(coss)).real
                    cc = np.conj(np.cos(phase)+1j*np.sin(phase))
                    FK[iw,ikx] *= cc
                else:
                    FK[iw,ikx] = complex(0.0,0.0)
                TK[itau,ikx] += FK[iw,ikx]
    # Normalize for inverse FFT
    TK /= nw
    return TK




def getVelocityProfile(vels_in,nlay,ns,ntr,dt,ts,firstz=0.,firstt=0.):
    """

    Map the layered velocity structure into the shape of the data.

    Parameters
    ---------
    vels_in: 2-D array of layered wave velocities (m/s) and layer thickness (m)
        Structure is velocities in first column, layer thicknesses in second
        If only one layer (i.e. constant velocity) input one layer with zero thickness.
    nlay: number of layers
    ns: number of samples
    ntr: number of traces
    dt: time step between samples (s)
    ts: 1-D array of sample times (s)
    firstz: uppermost depth (m)
    firstt: start time of trace (s)

    Output
    ---------
    vmig: 1-D array of migration velocities, same shape as ts. (m/s)
        If only one layer in input velocity array, output is constant.

    **
    The foundation of this script was taken from Matlab code written by Andreas Tzanis,
    Dept. of Geophysics, University of Athens (2005)
    **

    """

    # velocity layers
    layer_velocity  = vels_in[:,0]
    layer_thickness = vels_in[:,1]
    # constant velocity
    if nlay == 1:
        vmig = layer_velocity[0]
    elif nlay > 1:
        # Depth of layer boundaries from input thicknesses
        zs = np.max(layer_velocity)*ts/2.    # depth array for maximum possible penetration
        zboundaries = np.append(np.insert(np.cumsum(layer_thickness[:nlay-1]),0,0),max(zs))
        # Compute layer times from input velocity/thickness array (vels_in)
        # TODO: this could be shortened
        tofz_layers = np.empty((nlay+1))
        tofz_layers[0] = firstt
        for iz in np.arange(0,nlay-1):
            tofz_layers[iz+1] = tofz_layers[iz] + 2.0*layer_thickness[iz]/layer_velocity[iz]
        tofz_layers[nlay] = tofz_layers[nlay-1] + 2.0*(max(zs)-sum(layer_thickness[:-1]))/layer_velocity[nlay-1]
        # Interpolate to get t(z) for maximum penetration depth array
        from scipy.interpolate import interp1d
        tinterp = interp1d(zboundaries,tofz_layers)
        tofz = tinterp(zs)
        # Compute z(t) from t(z)
        zinterp = interp1d(tofz,zs)
        zoft = zinterp(np.arange(firstt,firstt+ns*dt,dt))
        vfz = layer_velocity[0]             # initial velocity
        vlz = layer_velocity[nlay-1]        # final velocity at depth z
        dz = np.mean(np.gradient(zs))
        lz  = firstz+(ns-1)*dz
        idx = np.argwhere(ts < tofz[0])             # out of range values
        zoft[idx] = 0.5*ts[idx]*vfz
        idx = np.argwhere(ts >= tofz[-1])
        zoft[idx] = lz + 0.5*(ts[idx] - tofz[-1])*vlz
        # Compute vmig(t) from z(t)
        vmig = np.empty((ns))
        for it in range(ns-1):
            vmig[it] = 2.0*(zoft[it+1] - zoft[it])/dt
        vmig[ns-1] = vmig[ns-2]
    return vmig


