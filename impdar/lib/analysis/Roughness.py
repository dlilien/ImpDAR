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
import matplotlib.pyplot as plt
import pickle
from scipy.signal import detrend,medfilt
from scipy.special import i0
from scipy.io import loadmat
from RadarFunctions import *

import numpy as np
import matplotlib.pyplot as plt
import glob,os
from scipy.signal import detrend,medfilt
from scipy.io import loadmat


# -----------------------------------------------------------------------------

# Load Files

FROM_FILES = True
if FROM_FILES:
    dataDir1 = '/Users/benhills/Drive/Research/Data/GroundRadarData/south_pole_lake/2018_2019/storadar/migtk/'
    dataDir2 = '/Users/benhills/Drive/Research/Data/GroundRadarData/south_pole_lake/ITASE/ITASE-2_07-4_to_S._Pole/'
    pickDir = './Pickfiles/'

    with open('./Pickfiles/MatchedLayers.pickle', 'rb') as f:
        layers_dict = pickle.load(f)
    with open('./Pickfiles/LineCoords.pickle', 'rb') as f:
        line_coords = pickle.load(f)
    bed_pick = layers_dict[len(layers_dict)-1]
    # Loop through all the files
    for i in range(19):
        if i == 0:
            fname = 's5_pole_7_051'
        else:
            fname = 's074_pole_'+str(i)+'_ch2'

        # Get the pick file
        with open(pickDir+fname+'picks.pickle', 'rb') as f:
            fdict = pickle.load(f)
        x = fdict['x']
        y = fdict['y']
        bd = fdict[len(fdict)-3][0]
        bpower = fdict[len(fdict)-3][1]
        bed_pick[fname] = np.array([[-99],bd,bpower,x,y])

# -----------------------------------------------------------------------------


# Constants
u_ice = 1.69e8     #m/s
freq = 3e6         #Hz
lam = u_ice/freq   #m
eps = 3.2          # relative permittivity

# Find window size
H = 2850.          # Approximate Ice thickness
D1 = np.sqrt(2.*lam*(H/np.sqrt(eps))) # Width of Fresnel zone
dx = 3.                         # m spacing between traces
N = int(round(D1/(2.*dx)))

# -----------------------------------------------------------------------------

ROUGH = True
if ROUGH:
    # Calculate 1-d roughness one line at a time
    for line in bed_pick.keys() - {'mean_depth'}:

        # get the surface elevation
        try:
            elev = bed_pick[line][3]
            # coordinates
            x = line_coords[line][0]
            y = line_coords[line][1]
        except:
            # get the data file for elevation
            try:
                migdata,elev,time,dist,vdist = loadStoMigData(dataDir1+line+'_migtk.mat',uice=168.,CReSIS=False)
                # coordinates
                x = line_coords[line][0]
                y = line_coords[line][1]
            except:
                migdata,elev,time,dist,vdist = loadStoMigData(dataDir2+line+'_mig.mat',uice=168.,CReSIS=False)
                # TODO: this is a hack to fix different geoids
                if line[:4] == 's074':
                    elev -= 18.5
                # coordinates
                if FROM_FILES:
                    x = bed_pick[line][3]
                    y = bed_pick[line][4]
                else:
                    x = bed_pick[line][0]
                    y = bed_pick[line][1]

        # Define the bed geometry
        h = bed_pick[line][1]
        bed_raw = elev - h
        bed = medfilt(bed_raw,401)

        # RMS bed roughness
        ED1 = np.empty((len(bed),))
        ED1[:] =  np.nan
        for n in range(N,len(bed)-N+1):
            z = bed[n-N:n+N]
            z = z[np.where(~np.isnan(z))]
            if len(z) <= 1:
                ED1[n] = np.nan
            else:
                z_ = detrend(z)
                total = 0.
                for i in range(len(z)):
                    total += (z_[i])**2.
                ED1[n] = np.sqrt((1/(len(z)-1.))*total)

        # Find the power reduction by Kirchoff theory
        g = 4*np.pi*ED1/lam
        b = (i0((g**2.)/2.))**2.
        pn = np.exp(-(g**2.))*b

        if line == 's5_pole_7_051':
            hold = 1. - pn[2723:3233]
            pn[2723:3233] = 1. - .1*hold

        power = bed_pick[line][2]
        Pc_filt = medfilt(power - 10.*np.log10(pn),401)
        bed_pick[line] = np.array([x,y,bed_raw,bed,h,power,Pc_filt,elev,ED1,pn])

    with open('./Pickfiles/BedPick_Roughness.pickle', 'wb') as f:
        pickle.dump(bed_pick,f)

# -----------------------------------------------------------------------------

PLOT = True
if PLOT:
    # get the roughness calcs from above
    with open('./Pickfiles/BedPick_Roughness.pickle', 'rb') as f:
        bed_pick = pickle.load(f)
    # loop through and plot
    X = np.array([])
    Y = np.array([])
    ED1 = np.array([])
    PN = np.array([])
    BP = np.array([])
    Pfilt = np.array([])
    for line in bed_pick.keys() - {'mean_depth'}:

        X = np.append(X,bed_pick[line][0])
        Y = np.append(Y,bed_pick[line][1])
        BP = np.append(BP,bed_pick[line][5])
        ED1 = np.append(ED1,bed_pick[line][8])
        PN = np.append(PN,bed_pick[line][9])
        Pfilt = np.append(Pfilt,bed_pick[line][6])
    """
    ax1 = plt.subplot(221)
    plt.scatter(X,Y,c=ED1,vmin=0,vmax=5)
    plt.colorbar()
    plt.title('Kirchhoff Roughness')
    ax2 = plt.subplot(222)
    plt.scatter(X,Y,c=PN,cmap='Blues',vmin=.6,vmax=1)
    plt.colorbar()
    plt.title('Power Reduction')
    ax3 = plt.subplot(223)
    plt.scatter(X,Y,c=BP - 10.*np.log10(PN),vmin=20,vmax=50,cmap='hot')
    plt.colorbar()
    plt.title('Corrected Power')
    """

    ax4 = plt.subplot(121)
    plt.scatter(X,Y,c=Pfilt,vmin=10,vmax=40,cmap='hot')
    ax2 = plt.subplot(122)
    plt.scatter(X,Y,c=BP,vmin=10,vmax=40,cmap='hot')
    plt.colorbar()
    plt.title('Filtered Power')


# for both the upper and lower grids
for d in ['Upper/','Lower/']:
    os.chdir("../../Data/"+d)

    # create output arrays
    x_out = dict.fromkeys(['along','across'],[np.array([]),np.array([])])
    y_out = dict.fromkeys(['along','across'],[np.array([]),np.array([])])
    z_out = dict.fromkeys(['along','across'],[np.array([]),np.array([])])

    # for all subdirectories at each grid
    for sub_d in os.listdir('.'):
        if os.path.isdir(sub_d) and sub_d[:2] == '20':
            os.chdir(sub_d)
        else:
            continue

        # for all files in the subdirectory
        for f in glob.glob("*_mig.mat"):
            try:
                # load the mat file with the migrated data in it and the pick file
                mfile = loadmat(f)
                file2 = f[:-4]+'_pick.mat'
                pfile = loadmat(file2)
                print "file loaded: " + f
            except:
                print 'file does not exist.'

            # data variables
            elev = mfile['elev'][0]
            dist = mfile['dist'][0]
            picktime = pfile['picks'][0][0][3][0]
            power = pfile['picks'][0][0][4][0]
            lat = pfile['geocoords'][0][0][0][0]
            lon = pfile['geocoords'][0][0][1][0]
            x_coord = pfile['geocoords'][0][0][3][0]
            y_coord = pfile['geocoords'][0][0][4][0]

            # distance in the file is not a real value, it is a counter
            # calculate distance using haversine formula
            lat*=np.pi/180.
            lon*=np.pi/180.
            R = 6371000.            #Radius of the earth
