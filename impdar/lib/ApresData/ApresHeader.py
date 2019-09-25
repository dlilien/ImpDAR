#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Header for ApRES data

This code is based on a series of Matlab scripts from Craig Stewart,
Keith Nicholls, and others.
The ApRES (Automated phase-sensitive Radio Echo Sounder) is a self-contained
instrument from BAS.

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

Sept 23 2019

"""

import numpy as np
import re

# --------------------------------------------------------------------------------------------

class ApresHeader():
    """
    Class for parameters from the header file.
    """
    def __init__(self):
        """Initialize data paramaters"""
        self.fsysclk = 1e9
        self.fs = 4e4
        self.fn = None
        self.header_string = None
        self.file_format = None
        self.noDwellHigh = None
        self.noDwellLow = None
        self.f0 = None
        self.f_stop = None
        self.ramp_up_step = None
        self.ramp_down_step = None
        self.tstep_up = None
        self.tstep_down = None
        self.snum = None
        self.nsteps_DDS = None
        self.chirp_length = None
        self.chirp_grad = None
        self.nchirp_samples = None
        self.ramp_dir = None

    # --------------------------------------------------------------------------------------------

    def read_header(self,fn_apres,max_header_len=2000):
        """
        Read the header string, to be partitioned later

        Parameters
        ---------
        fn_apres: string
            file name to update with
        max_header_len: int
            maximum length of header to read (can be too long)

        Output
        ---------
        """
        self.fn = fn_apres
        fid = open(fn_apres,'rb')
        self.header_string = str(fid.read(max_header_len))
        fid.close()

    # --------------------------------------------------------------------------------------------

    def get_file_format(self):
        """
        Determine fmcw file format from burst header using keyword presence
        There are a few different formats through the years.

        ### Original Matlab script Notes ###
        Craig Stewart
        2013-10-20
        Updated by Keith Nicholls, 2014-10-22: RMB2
        """

        if 'SW_Issue=' in self.header_string: # Data from RMB2 after Oct 2014
            self.file_format = 5
        elif 'SubBursts in burst:' in self.header_string: # Data from after Oct 2013
            self.file_format = 4
        elif '*** Burst Header ***' in self.header_string: # Data from Jan 2013
            self.file_format = 3
        elif 'RADAR TIME' in self.header_string: # Data from Prototype FMCW radar (nov 2012)
            self.file_format = 2
        else:
            raise TypeError('Unknown file format - check file')

    def update_parameters(self,fn_apres=None):
        """
        Update the parameters with the apres file header


        ### Original Matlab Notes ###
        Extract from the hex codes the actual paramaters used by RMB2
        The contents of config.ini are copied into a data header.
        Note this script assumes that the format of the hex codes have quotes
        e.g. Reg02="0D1F41C8"

        Checks for a sampling frequency of 40 or 80 KHz.  Apart from Lai Bun's
        variant (WDC + Greenland) it could be hard coded to 40 KHz.

        However, there is no check made by the system that the N_ADC_SAMPLES
        matches the requested chirp length

        NOT COMPLETE - needs a means of checking for profile mode, where multiple sweeps
        per period are transmitted- see last line
        """

        if self.header_string is None:
            if fn_apres is None:
                raise TypeError('Must input file name if the header has not been read yet.')
            else:
                self.read_header(fn_apres)
        if self.file_format is None:
            self.get_file_format()

        loc1 = [m.start() for m in re.finditer('Reg0', self.header_string)]
        loc2 = [m.start() for m in re.finditer('="', self.header_string)]

        for k in range(len(loc1)):
            case = self.header_string[loc1[k]:loc2[k]]

            if case == 'Reg01':
                # Control Function Register 2 (CFR2) Address 0x01 Four bytes
                # Bit 19 (Digital ramp enable)= 1 = Enables digital ramp generator functionality.
                # Bit 18 (Digital ramp no-dwell high) 1 = enables no-dwell high functionality.
                # Bit 17 (Digital ramp no-dwell low) 1 = enables no-dwell low functionality.
                # With no-dwell high, a positive transition of the DRCTL pin initiates a positive slope ramp, which
                # continues uninterrupted (regardless of any activity on the DRCTL pin) until the upper limit is reached.
                # Setting both no-dwell bits invokes a continuous ramping mode of operation;
                loc3 = self.header_string[loc2[k]+2:].find('"')
                val = self.header_string[loc2[k]+2:loc2[k]+loc3+2]
                val = bin(int(val, 16))
                val = val[::-1]
                self.noDwellHigh = int(val[18])
                self.noDwellLow = int(val[17])

            #elif case == 'Reg08':
            #    # Phase offset word Register (POW) Address 0x08. 2 Bytes dTheta = 360*POW/2^16.
            #    val = char(reg{1,2}(k));
            #    H.phaseOffsetDeg = hex2dec(val(1:4))*360/2^16;

            elif case == 'Reg0B':
                # Digital Ramp Limit Register Address 0x0B
                # Digital ramp upper limit 32-bit digital ramp upper limit value.
                # Digital ramp lower limit 32-bit digital ramp lower limit value.
                loc3 = self.header_string[loc2[k]+2:].find('"')
                val = self.header_string[loc2[k]+2:loc2[k]+loc3+2]
                self.f0 = int(val[8:], 16)*self.fsysclk/(2**32)
                self.f_stop = int(val[:8], 16)*self.fsysclk/(2**32)

            elif case == 'Reg0C':
                # Digital Ramp Step Size Register Address 0x0C
                # Digital ramp decrement step size 32-bit digital ramp decrement step size value.
                # Digital ramp increment step size 32-bit digital ramp increment step size value.
                loc3 = self.header_string[loc2[k]+2:].find('"')
                val = self.header_string[loc2[k]+2:loc2[k]+loc3+2]
                self.ramp_up_step = int(val[8:], 16)*self.fsysclk/(2**32)
                self.ramp_down_step = int(val[:8], 16)*self.fsysclk/(2**32)

            elif case == 'Reg0D':
                # Digital Ramp Rate Register Address 0x0D
                # Digital ramp negative slope rate 16-bit digital ramp negative slope value that defines the time interval between decrement values.
                # Digital ramp positive slope rate 16-bit digital ramp positive slope value that defines the time interval between increment values.
                loc3 = self.header_string[loc2[k]+2:].find('"')
                val = self.header_string[loc2[k]+2:loc2[k]+loc3+2]
                self.tstep_up = int(val[4:], 16)*4/self.fsysclk
                self.tstep_down = int(val[:4], 16)*4/self.fsysclk

        strings = ['SamplingFreqMode=','N_ADC_SAMPLES=']
        output = np.empty((len(strings))).astype(str)
        for i,string in enumerate(strings):
            if string in self.header_string:
                search_start = self.header_string.find(string)
                search_end = self.header_string[search_start:].find('\\')
                output[i] = self.header_string[search_start+len(string):search_end+search_start]

        self.fs = output[0]
        if self.fs == 1:        # if self.fs > 70e3:
            self.fs = 8e4       #     self.fs = 80e3
        else:                   # else
            self.fs = 4e4       #     self.fs = 40e3

        self.snum = int(output[1])

        self.nsteps_DDS = round(abs((self.f_stop - self.f0)/self.ramp_up_step)) # abs as ramp could be down
        self.chirp_length = int(self.nsteps_DDS * self.tstep_up)
        self.nchirp_samples = round(self.chirp_length * self.fs)

        # If number of ADC samples collected is less than required to collect
        # entire chirp, set chirp length to length of series actually collected
        if self.nchirp_samples > self.snum:
            self.chirp_length = self.snum / self.fs

        self.chirp_grad = 2.*np.pi*(self.ramp_up_step/self.tstep_up) # chirp gradient (rad/s/s)
        if self.f_stop > 400e6:
            self.ramp_dir = 'down'
        else:
            self.ramp_dir = 'up'

        if self.noDwellHigh and self.noDwellLow:
            self.ramp_dir = 'upDown'
            self.nchirpsPerPeriod = np.nan # self.nchirpSamples/(self.chirpLength)

# --------------------------------------------------------------------------------------------

    def to_matlab(self):
        """Convert all associated attributes into a dictionary formatted for use with :func:`scipy.io.savemat`
        """
        outmat = {att: getattr(self, att) for att in vars(self)}
        return outmat

    def from_matlab(self, matlab_struct):
        """Associate all values from an incoming .mat file (i.e. a dictionary from :func:`scipy.io.loadmat`) with appropriate attributes
        """
        for attr, attr_dim in zip(self.attrs, self.attr_dims):
            setattr(self, attr, matlab_struct[attr][0][0][0])
            # Use this because matlab inputs may have zeros for flags that
            # were lazily appended to be arrays, but we preallocate
            if attr_dim is not None and getattr(self, attr).shape[0] == 1:
                setattr(self, attr, np.zeros((attr_dim, )))

        for attr in self.bool_attrs:
            setattr(self, attr, True if matlab_struct[attr][0][0][0] == 1 else 0)
