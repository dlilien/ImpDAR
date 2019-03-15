#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 11:01:05 2019

@author: benhills
"""

import numpy as np
import time

# -----------------------------------------------------------------------------

def migrationKirchhoff(xs,ts,D,vel=1.68e8,nearfield=False):
    print('Kirchhoff Migration (diffraction summation) of %.0fx%.0f matrix'%(len(xs),len(ts)))
    # start the timer
    start = time.time()
    # Calculate the time derivative of the input data
    gradD = np.gradient(D,ts,axis=0)
    # Create an empty array to fill with migrated data
    migD = np.zeros_like(D)
    # Loop through all traces
    for xi in range(len(D[0])):
        print('Migrating trace number:',xi)
        # get the trace distance
        x = xs[xi]
        # Loop through all samples
        for ti in range(len(D)):
            # get the sample time        
            t = ts[ti]
            # convert to depth
            z = vel*t/2.
            # get the radial distances between input point and output point
            rs = np.sqrt((xs-x)**2.+z**2.)
            # find the cosine of the angle of the tangent line 
            # correct for 'obliquity factor'
            costheta = z/rs
            # get the exact indices from the array (closest to rs)
            Didx = [np.argmin(abs(ts-2.*r/vel)) for r in rs]
            # integrate the farfield term
            gradDhyp = np.array([gradD[Didx[i],i] for i in range(len(Didx))])
            gradDhyp[2.*rs/vel>max(ts)] = 0.    # zero points that are outside of the domain
            integral = np.nansum(gradDhyp*costheta/vel)  # TODO: Yilmaz eqn 4.5 has an extra r in this weight factor???
            # integrate the nearfield term
            if nearfield == True:
                Dhyp = np.array([D[Didx[i],i] for i in range(len(Didx))])
                Dhyp[2.*rs/vel>max(ts)] = 0.    # zero points that are outside of the domain
                integral += np.nansum(Dhyp*costheta/rs**2.)
            # sum the integrals and output
            migD[ti,xi] = 1/(2.*np.pi)*integral
    # print the total time
    print('Kirchhoff Migration of %.0fx%.0f matrix complete in %.2f seconds'
          %(len(xs),len(ts),time.time()-start))
    return migD

# -----------------------------------------------------------------------------

def migrationStolt(xs,ts,D,vel=1.68e8):
    print('Stolt Migration (f-k migration) of %.0fx%.0f matrix'%(len(xs),len(ts)))
    # save the start time
    start = time.time()
    # pad the array with zeros up to the next power of 2 for discrete fft
    nt = 2**(np.ceil(np.log(len(ts))/np.log(2))).astype(int)
    nx = 2**(np.ceil(np.log(len(xs))/np.log(2))).astype(int)
    # 2D Forward Fourier Transform to get S(kx,z=0,omega)
    S = np.fft.fft2(D,(nt,nx))
    # get the temporal frequencies
    dt = np.mean(np.gradient(ts))
    omega = np.fft.fftfreq(nt, d=dt)
    # get the horizontal wavenumbers
    xstep = xs[1]-xs[0]
    kx = np.fft.fftfreq(nx, d=xstep)
    # interpolate from frequency (omega) into wavenumber (kz)
    from scipy.interpolate import interp2d
    interp_real = interp2d(kx,omega,S.real)
    interp_imag = interp2d(kx,omega,S.imag)
    Sinterp = np.zeros_like(S)
    print('Interpolating from temporal frequency ($\omega$) to vertical wavenumber (kz):')
    # for all temporal frequencies
    for zj in range(len(omega)):
        kzj = omega[zj]*2./vel
        if zj%100 == 0:
            print(round(kzj,2),'1/m')
        # for all horizontal wavenumbers
        for xi in range(len(kx)):
            kxi = kx[xi]
            # migration conversion to wavenumber (Yilmaz equation C.53)
            omegaj = vel/2.*np.sqrt(kzj**2.+kxi**2.)
            # get the interpolated FFT values, real and imaginary, S(kx,kz,t=0)
            Sinterp[zj,xi] = interp_real(kxi,omegaj) + 1j*interp_imag(kxi,omegaj)
    # all vertical wavenumbers
    kz = omega*2./vel
    # grid wavenumbers for scaling calculation
    kX,kZ = np.meshgrid(kx,kz)
    # scaling for obliquity factor (Yilmaz equation C.56)    
    scaling = kZ/np.sqrt(kX**2.+kZ**2.)
    Sinterp *= scaling
    # the DC frequency should be 0.
    Sinterp[0,0] = 0. + 1j*0.
    # 2D Inverse Fourier Transform to get D(x,z,t=0)
    migD = np.real(np.fft.ifft2(Sinterp))
    # Cut back to input matrix dimensions
    migD = migD[:len(ts),:len(xs)]
    # print the total time
    print('Stolt Migration of %.0fx%.0f matrix complete in %.2f seconds'
          %(len(xs),len(ts),time.time()-start))    
    return migD
    
# -----------------------------------------------------------------------------



def main(xs,ts,D,vels_in):
    """ 
     GAZDAGMIG : Gazdag phase-shifting migration for constant or layered 
                 velocity structures.  
           ==>   For the sake of speed, the construction of the image is 
                 usually performed with the FORTRAN '90 MEX-file 
                 "gazdag.f90". For MS Windows OS, the MEX-file is provided 
                 ready to use (gazdag.dll). For Linux OS or any other flavour 
                 of Unix, the MEX-file should be built by the user. At any 
                 rate,(very much) slower native M-code to perform the same 
                 tasks is attached in the subfunctions "gazdagcv" and 
                 "gazdaglv". This will take over if the MEX-file is not 
                 available.   
     
        Usage  : dmig  = gazdagmig(d, dt, dx, vofh) 
     
       Inputs  :  
            d  : The zero- or common-offset GPR section  
           dt  : The sampling rate  
           dx  : The trace spacing (spatial sampling rate)  
      vofh(n,2): The 1-D velocity model of "nlay" velocity - 
                 thickness pairs, for example: 
                 [ 0.1    1 ;       ... 1st layer 
                   0.08   2 ;       ... 2nd layer 
                   ... 
                   0.18   0 ]       ... n'th layer==basal halfspace. 
                  A uniform halfspace is given as a single layer stucture with 
                  zero thickness. Velovity values are given in m/ns and 
                  thichnesses in m.  
     
       Outputs  :  
          dmig  : The migrated section 
     
       Requires : yxtoxy.m 
     
          Uses  : FORTRAN'90 MEX function "gazdag.<mex> with source code 
                  "gazdag.f90", and the attached subfunctions "gazdagcv" and 
                  "gazdaglv"  
     
      Author    : Andreas Tzanis, 
                  Dept. of Geophysics, 
                  University of Athens 
                  atzanis@geol.uoa.gr 
     
      Copyright (C) 2005, Andreas Tzanis. All rights reserved. 
     
        This program is free software; you can redistribute it and/or modify 
        it under the terms of the GNU General Public License as published by 
        the Free Software Foundation; either version 2 of the License, or 
        (at your option) any later version. 
     
        This program is distributed in the hope that it will be useful, 
        but WITHOUT ANY WARRANTY; without even the implied warranty of 
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
        GNU General Public License for more details. 
     
        You should have received a copy of the GNU General Public License 
        along with this program; if not, write to the Free Software 
        Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. 
    """
    
    # Data structure
    [ns,ntr]=np.shape(D)
    # Velocity structure from input
    [nlay, vpairs]=np.shape(vels_in)
    dt = np.mean(np.gradient(ts))
    dx = np.mean(np.gradient(xs))
    vmig = getVelocityProfile(vels_in,ns,ntr,dt,ts)
    # Fourier transform to frequency-wavenumber (FKx) domain 
    fk  = np.fft.fft2(D)
    fk  = np.fft.fftshift(fk)
    # Migration by phase shift
    if nlay == 1:
        img = phaseShiftConstantVel(ns, ntr, dt, dx, vmig, fk)
    elif nlay >1:
        img = phaseShiftLayeredVel(ns, ntr, dt, dx, vmig, fk)
    # Transform from time-wavenumber (TKx) to time-space (TX) domain to get migrated section 
    migD = np.fft.ifft(np.fft.ifftshift(img,1)).real
    return migD





# -----------------------------------------------------------------------------
    
def getVelocityProfile(vs_in,nlay,ns,ntr,dt,ts,firstz=0.,firstt=0.):
    # velocity layers
    layer_velocity  = vs_in[:,0]
    layer_thickness = vs_in[:,1] 
    # Compute migration velocity vmig(t) from v(z)  
    if nlay == 1: 
        vmig = layer_velocity[0]
    elif nlay > 1:
        zmax = np.max(vs_in[:,0])*ts[-1]/2.    # max possible penetration 
        dz = (zmax-firstz)/(ns-1)
        zs = np.arange(firstz,zmax+dz,dz)
        if ns != len(zs):
            raise Exception('Error creatind depth array from velocities.')
        # Compute t(z) from V(z) 
        tofz_layers = np.empty((nlay+1))
        tofz_layers[0] = firstt
        for iz in np.arange(0,nlay-1):
            tofz_layers[iz+1] = tofz_layers[iz] + 2.0*layer_thickness[iz]/layer_velocity[iz]
        tofz_layers[nlay] = tofz_layers[nlay-1] + 2.0*(zmax-sum(layer_thickness[:-1]))/layer_velocity[nlay-1]
        zboundaries = np.append(np.insert(np.cumsum(layer_thickness[:nlay-1]),0,0),zmax)
        # interpolation
        from scipy.interpolate import interp1d
        tinterp = interp1d(zboundaries,tofz_layers)
        tofz = tinterp(zs)
        zinterp = interp1d(tofz,zs)
        zoft = zinterp(np.arange(firstt,firstt+ns*dt,dt))
        # Compute z(t) from t(z)
        vfz = layer_velocity[0]             # initial velocity 
        vlz = layer_velocity[nlay-1]        # final velocity at depth z 
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

# -----------------------------------------------------------------------------
 
def phaseShiftConstantVel(ns, ntr, dt, dx, vmigc, fk):
    """
    This is the MATLAB code to compute the FKx -> TKx image of Gazdag 
    Phase-shifting migration for zero-offset GPR data in uniform halfspaces 
    (migration velocity  vmigc is constant). 
    *** The active implementation is vectorized (as much as the 
        integration process will allow and is reasonably fast. It will do a 
        [512 x 512] job in about 36 seconds (as opposed to the 5 seconds      needed for the MEX-file which works with typical -FORTRAN- 
        sequential programing).  
    *** The inactive (commented out) implementation is traditional sequential 
        programing (the same as in FORTRAN). It is VERY SLOW! 
    *** Note that the continuation operator "cc" is conjugated, opposite to 
        what books usually say), to account for the (engineering) definition 
        of the Fourier kernel in MATLAB   
    
    Author: Andreas Tzanis (C)  
            November 2005 
    """
    # Set up physical constants  
    nw, w0, dw = ns, -np.pi/dt, 2.*np.pi/(ns*dt)
    nx, kx0 = ntr, -np.pi/dx
    dkx = 2.*np.pi/(nx*dx)
    ntau, dtau = ns, dt
    # initialize image array
    img = np.zeros((ns,ntr))+0j
     
    kx = np.arange(kx0,-kx0,dkx)
    for iw in range(nw):
        w = w0 + (iw)*dw 
        if w == 0.0:
            w = 1e-10/dt
        vkx2 = (vmigc*vmigc * kx**2.)/4.
        ik = np.argwhere(vkx2 < w**2.)
        ffk = fk[iw,ik]
        phase = (-w*dtau*np.sqrt(1.0 - vkx2[ik]/w**2.)).real
        cc    = np.conj(np.cos(phase)+1j*np.sin(phase))
        # Accumulate image summed over all frequencies 
        for itau in range(ntau):
             ffk *= cc
             img[itau,ik] += ffk
    # Normalize for inverse FFT 
    img = img/nw; 
    return img

# -----------------------------------------------------------------------------

def phaseShiftLayeredVel(ns, ntr, dt, dx, vmigv, fk):
    """
    
    This is the MATLAB code to compute the FKx -> TKx image of Gazdag 
    Phase-shifting migration for zero-offset GPR data in LAYERED halfspaces 
    (migration velocity  vmigc is a function of time). 
    *** This routine is programmed sequentially, as in FORTRAN (same code 
        with gazdag.f90 and is DAMNED SLOW! 
    *** Note that the continuation operator "cc" is conjugated, (opposite to 
        what books usually say), to account for the (engineering) definition 
        of the Fourier kernel in MATLAB    
    
    Author: Andreas Tzanis (C)  
            November 2005 
    
    """
    # Set up physical constants  
    nw, w0, dw = ns, -np.pi/dt, 2.*np.pi/(ns*dt)
    nx, kx0 = ntr, -np.pi/dx
    dkx = 2.*np.pi/(nx*dx)
    ntau, dtau = ns, dt
    ft = 0
    ftau, tmax = ft, ft+(ntau-1)*dtau
    # initialize image array
    img = np.zeros((ns,ntr))+0j
    
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
                    fk[iw,ikx] *= cc
                else:
                    fk[iw,ikx] = complex(0.0,0.0); 
                img[itau,ikx] += fk[iw,ikx]
    # Normalize for inverse FFT 
    img /= nw; 
    return img