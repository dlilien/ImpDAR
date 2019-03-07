"""
@author Joshua Driscol
@date 03/04/2019

This script removes the pretrigger from raw .mat files (assuming that there is one).
The raw .mat files are generated using `impdar load pe _____` for a PulseEkko radar, for example.
For every file in the folder, read in the .mat file, access its attributes as dictionary keys.
Then update the data array, snum, and travel_time attributes.

"""

from scipy.io import loadmat, savemat
import numpy as np
import os

this_dir = './'

#remove pretrigger
for subdir, dirs, files in os.walk(this_dir):
    for file in files:
        if file.endswith('.mat'):
            mat_file = loadmat(file)
            trig = int(mat_file['trig'])
            data = mat_file['data']
            snum = mat_file['snum']
            travel_time = mat_file['travel_time']
            
            mat_file['data'] = data[trig:,:]
            mat_file['snum'] = snum - trig + 1
            mat_file['travel_time'] = travel_time[1:-trig+1]
            savemat(file, mat_file)
