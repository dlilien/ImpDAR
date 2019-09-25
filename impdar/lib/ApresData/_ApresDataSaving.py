#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Load ApRES data

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
import datetime
from . import ApresData
from . import ApresFlags

# --------------------------------------------------------------------------------------------

class apres_parameters:
    """
    Class for parameters from the header file.
    """
    def __init__(self):
        """Initialize data paramaters"""
        self.fsysclk = 1e9
        self.fs = 4e4

    def update_parameters(self,fn_apres,max_header_len=2000):
        """
        Update the parameters with the apres file header

        Parameters
        ---------
        fn_apres: string
            file name to update with
        max_header_len: int
            maximum length of header to read (can be too long)

        Output
        ---------

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
        fid = open(fn_apres,'rb')
        header = str(fid.read(max_header_len))
        fid.close()
        loc1 = [m.start() for m in re.finditer('Reg0', header)]
        loc2 = [m.start() for m in re.finditer('="', header)]

        for k in range(len(loc1)):
            case = header[loc1[k]:loc2[k]]

            if case == 'Reg01':
                # Control Function Register 2 (CFR2) Address 0x01 Four bytes
                # Bit 19 (Digital ramp enable)= 1 = Enables digital ramp generator functionality.
                # Bit 18 (Digital ramp no-dwell high) 1 = enables no-dwell high functionality.
                # Bit 17 (Digital ramp no-dwell low) 1 = enables no-dwell low functionality.
                # With no-dwell high, a positive transition of the DRCTL pin initiates a positive slope ramp, which
                # continues uninterrupted (regardless of any activity on the DRCTL pin) until the upper limit is reached.
                # Setting both no-dwell bits invokes a continuous ramping mode of operation;
                loc3 = header[loc2[k]+2:].find('"')
                val = header[loc2[k]+2:loc2[k]+loc3+2]
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
                loc3 = header[loc2[k]+2:].find('"')
                val = header[loc2[k]+2:loc2[k]+loc3+2]
                self.startFreq = int(val[8:], 16)*self.fsysclk/(2**32)
                self.stopFreq = int(val[:8], 16)*self.fsysclk/(2**32)

            elif case == 'Reg0C':
                # Digital Ramp Step Size Register Address 0x0C
                # Digital ramp decrement step size 32-bit digital ramp decrement step size value.
                # Digital ramp increment step size 32-bit digital ramp increment step size value.
                loc3 = header[loc2[k]+2:].find('"')
                val = header[loc2[k]+2:loc2[k]+loc3+2]
                self.rampUpStep = int(val[8:], 16)*self.fsysclk/(2**32)
                self.rampDownStep = int(val[:8], 16)*self.fsysclk/(2**32)

            elif case == 'Reg0D':
                # Digital Ramp Rate Register Address 0x0D
                # Digital ramp negative slope rate 16-bit digital ramp negative slope value that defines the time interval between decrement values.
                # Digital ramp positive slope rate 16-bit digital ramp positive slope value that defines the time interval between increment values.
                loc3 = header[loc2[k]+2:].find('"')
                val = header[loc2[k]+2:loc2[k]+loc3+2]
                self.tstepUp = int(val[4:], 16)*4/self.fsysclk
                self.tstepDown = int(val[:4], 16)*4/self.fsysclk

        strings = ['SamplingFreqMode=','N_ADC_SAMPLES=']
        output = np.empty((len(strings))).astype(str)
        for i,string in enumerate(strings):
            if string in header:
                search_start = header.find(string)
                search_end = header[search_start:].find('\\')
                output[i] = header[search_start+len(string):search_end+search_start]

        self.fs = output[0]
        if self.fs == 1:        # if self.fs > 70e3:
            self.fs = 8e4       #     self.fs = 80e3
        else:                   # else
            self.fs = 4e4       #     self.fs = 40e3

        self.snum = int(output[1])

        self.nstepsDDS = round(abs((self.stopFreq - self.startFreq)/self.rampUpStep)) # abs as ramp could be down
        self.chirpLength = int(self.nstepsDDS * self.tstepUp)
        self.nchirpSamples = round(self.chirpLength * self.fs)

        # If number of ADC samples collected is less than required to collect
        # entire chirp, set chirp length to length of series actually collected
        if self.nchirpSamples > self.snum:
            self.chirpLength = self.snum / self.fs

        self.K = 2.*np.pi*(self.rampUpStep/self.tstepUp) # chirp gradient (rad/s/s)
        if self.stopFreq > 400e6:
            self.rampDir = 'down'
        else:
            self.rampDir = 'up'

        if self.noDwellHigh and self.noDwellLow:
            self.rampDir = 'upDown'
            self.nchirpsPerPeriod = np.nan # self.nchirpSamples/(self.chirpLength)

# --------------------------------------------------------------------------------------------

def file_format(fn_apres,max_header_len=2000):
    """
    Determine fmcw file format from burst header using keyword presence
    There are a few different formats through the years.

    Parameters
    ---------
    fn_apres: string
        file name to check
    max_header_len: int
        length to read looking for header (can be too long)

    Output
    ---------
    fmt: int
        file format (1-5)

    ### Original Matlab script Notes ###
    Craig Stewart
    2013-10-20
    Updated by Keith Nicholls, 2014-10-22: RMB2
    """

    with open(fn_apres,'rb') as fid:
        header = str(fid.read(max_header_len))
    if 'SW_Issue=' in header: # Data from RMB2 after Oct 2014
        fmt = 5
    elif 'SubBursts in burst:' in header: # Data from after Oct 2013
        fmt = 4
    elif '*** Burst Header ***' in header: # Data from Jan 2013
        fmt = 3
    elif 'RADAR TIME' in header: # Data from Prototype FMCW radar (nov 2012)
        fmt = 2
    else:
        TypeError('Unknown file format - check file')
    return fmt


# -----------------------------------------------------------------------------------------------------

def load_burst_rmb(apres_data,fn_apres,burst=1,fs=40000,max_header_len=2000,burst_pointer=0):
    """
    Load bursts from the apres acquisition.
    Normally, this should be called from the load_apres function.

    Parameters
    ---------
    apres_data: class
        data object
    fn_apres: string
        file name to load
    burst: int
        number of bursts to load
    fs: int
        sampling frequency
    max_header_len: int
        maximum length to read for header (can be too long)
    burst_pointer: int
        where to start reading the file for bursts

    Output
    ---------
    apres_data: class
        data object

    ### Original Matlab Script Notes ###
    Read FMCW data file from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)

    Corrected so that Sampling Frequency has correct use (ie, not used in
    this case)
    """

    apres_data.file_read_code = 0

    try:
        fid = open(fn_apres,'rb')
    except:
        # Unknown file
        apres_data.file_read_code = -1;
        TypeError('Cannot read file', fn_apres)

    # Get the total length of the file
    fid.seek(0,2)
    file_len = fid.tell()
    burst_count = 1

    # Read bursts in a loop
    while burst_count <= burst and burst_pointer <= file_len - max_header_len:
        # Go to burst pointer and read the header for the burst
        fid.seek(burst_pointer)
        header = str(fid.read(max_header_len))

        try:
            # Read header values
            strings = ['N_ADC_SAMPLES=','NSubBursts=','Average=','nAttenuators=','Attenuator1=',
                      'AFGain=','TxAnt=','RxAnt=']
            output = np.empty((len(strings))).astype(str)
            for i,string in enumerate(strings):
                if string in header:
                    search_start = header.find(string)
                    search_end = header[search_start:].find('\\')
                    output[i] = header[search_start+len(string):search_end+search_start]
            # Write header values to data object
            apres_data.snum = int(output[0])
            apres_data.n_subbursts = int(output[1])
            apres_data.average = int(output[2])
            apres_data.n_attenuators = int(output[3])
            apres_data.attenuator1 = np.array(output[4].split(',')).astype(int)[apres_data.n_attenuators-1]
            apres_data.attenuator2 = np.array(output[5].split(',')).astype(int)[apres_data.n_attenuators-1]
            apres_data.tx_ant = np.array(output[6].split(',')).astype(int)
            apres_data.rx_ant = np.array(output[7].split(',')).astype(int)

            apres_data.tx_ant = apres_data.tx_ant[apres_data.tx_ant==1]
            apres_data.rx_ant = apres_data.rx_ant[apres_data.rx_ant==1]

            if apres_data.average != 0:
                apres_data.chirps_in_burst = 1
            else:
                apres_data.chirps_in_burst = apres_data.n_subbursts*len(apres_data.tx_ant)*\
                                            len(apres_data.rx_ant)*apres_data.n_attenuators

            search_string = '*** End Header ***'
            search_ind = header.find(search_string)
            burst_pointer += search_ind + len(search_string)

        except:
            apres_data.file_read_code = -2
            apres_data.burst = burst_count

        words_per_burst = apres_data.chirps_in_burst*apres_data.snum

        if burst_count < burst and burst_pointer <= file_len - max_header_len:
            if apres_data.average != 0:
                burst_pointer += apres_data.chirps_in_burst*apres_data.snum*4
            else:
                burst_pointer += apres_data.chirps_in_burst*apres_data.snum*2

        burst_count += 1

    # --- Section --- #

    strings = ['Time stamp=','Temp1=','Temp2=','BatteryVoltage=']
    output = np.empty((len(strings))).astype(str)

    for i,string in enumerate(strings):
        if string in header:
            search_start = header.find(string)
            search_end = header[search_start:].find('\\')
            output[i] = header[search_start+len(string):search_end+search_start]

    if 'Time stamp' not in header:
        apres_data.file_read_code = -4
    else:
        apres_data.time_stamp = datetime.datetime.strptime(output[0],'%Y-%m-%d %H:%M:%S')

    apres_data.temperature1 = float(output[1])
    apres_data.temperature2 = float(output[2])
    apres_data.battery_voltage = float(output[3])

    # --- Section --- #

    if burst_count == burst+1:
        if apres_data.average == 2:
            apres_data.data = fid.read(words_per_burst)
        elif apres_data.average == 1:
            fid.seek(burst_pointer+1)
            apres_data.data = fid.read(words_per_burst)
        else:
            apres_data.data = np.fromfile(fid,dtype='uint16',count=words_per_burst)
            #apres_data.data = fid.read(16)
        if fid.tell()-(burst_pointer-1) < words_per_burst:
            apres_data.file_read_code = 2

        apres_data.data[apres_data.data<0] = apres_data.data[apres_data.data<0] + 2**16.
        apres_data.data *= int(2.5/2**16.)

        if apres_data.average == 2:
            apres_data.data /= (apres_data.n_subbursts*apres_data.n_attenuators)

        apres_data.start_ind = np.transpose(np.arange(0,apres_data.snum*apres_data.chirps_in_burst,apres_data.snum))
        apres_data.end_ind = apres_data.start_ind + apres_data.snum-1
        apres_data.burst = burst

    else:
        apres_data.burst = burst_count -1
        apres_data.file_read_code = -4

    fid.close()

    # Clean temperature record (wrong data type?)
    if apres_data.temperature1 > 300:
        apres_data.temperature1 -= 512
    if apres_data.temperature2 > 300:
        apres_data.temperature2 -= 512

    apres_data.data_code = 0

    return apres_data

# -----------------------------------------------------------------------------------------------------

def load_apres(fn_apres,burst=1,fs=40000, *args, **kwargs):
    """
    # Craig Stewart
    # 2013 April 24
    # 2013 September 30 - corrected error in vif scaling
    # 2014/5/20 time stamp moved here from fmcw_derive_parameters (so that this
    # is not overwritted later)
    # 2014/5/21 changed how radar chirp is defined (now using chirp gradient as
    # fundamental parameter)
    # 2014/5/22 fixed bug in chirptime
    # 2014/8/21 moved make odd length to external (called from fmcw_range)
    # 2014/10/22 KWN - edited to allow for new headers in RMB2 files
    """


    ## Load data and reshape array
    if fn_apres[-4:] == '.mat':
        apres_data = ApresData(fn_apres)
    else:
        apres_data = ApresData(None)
        apres_data.file_format = file_format(fn_apres)
        apres_data.data = load_burst_rmb(apres_data, fn_apres, burst, fs, apres_data.file_format)

    H = apres_parameters()

    # Extract just good chirp data from voltage record and rearrange into
    # matrix with one chirp per row
    # note: you can't just use reshape as we are also cropping the 20K samples
    # of sync tone etc which occur after each 40K of chirp.
    AttSet = apres_data.attenuator1 + 1j*apres_data.attenuator2 # unique code for attenuator setting

    ## Add metadata to structure

    # Sampling parameters
    if ischar(apres_data.file_format):
        apres_data.er = 3.18;
        apres_data.dt = 1/apres_data.fs;
        apres_data.ci = 3e8/np.sqrt(apres_data.er);
        apres_data.lambdac = apres_data.ci/apres_data.fc;
        # Load each chirp into a row

        apres_data.vif = np.zeros(apres_data.ChirpsInBurst,apres_data.SamplesPerChirp); # preallocate array
        chirpInterval = 1.6384/(24*3600); # days
        apres_data.Endind = apres_data.Startind + apres_data.SamplesPerChirp - 1;
        for chirp in range(apres_data.ChirpsInBurst):
            apres_data.vif[chirp,:] = apres_data.data[apres_data.Startind[chirp]:apres_data.Endind[chirp]]
            apres_data.chirpNum[chirp,0] = chirp                                                # chirp number in burst
            apres_data.chirpAtt[chirp,0] = AttSet[1+mod(chirp-1,numel(AttSet))]                 # attenuator setting for chirp
            apres_data.chirpTime[chirp,0] = apres_data.TimeStamp + chirpInterval*(chirp-1)      # time of chirp
    else:
        apres_data.SamplesPerChirp = apres_data.snum
        apres_data.fs = 4e4         # sampling frequency
        apres_data.f0 = 2e8         # start frequency
        #apres_data.fc = 3e8        # start frequency
        apres_data.K = 2*np.pi*2e8  # chirp gradient in rad/s/s (200MHz/s)
        #apres_data.f0 = apres_data.f0 + (apres_data.K/(4*pi))/apres_data.fs; # start frequency
        apres_data.processing = {}

        H.update_parameters(fn_apres)

        if apres_data.file_format == 5 or apres_data.file_format == 4:
            apres_data.K = H.K
            apres_data.f0 = H.startFreq
            apres_data.fs = H.fs
            apres_data.f1 = H.startFreq + H.chirpLength * H.K/2/np.pi
            apres_data.SamplesPerChirp = round(H.chirpLength * H.fs);
            apres_data.T = H.chirpLength
            apres_data.B = H.chirpLength * H.K/2/np.pi
            apres_data.fc = H.startFreq + apres_data.B/2
            apres_data.dt = 1./H.fs
            apres_data.er = 3.18
            apres_data.ci = 3e8/np.sqrt(apres_data.er);
            apres_data.lambdac = apres_data.ci/apres_data.fc;
            apres_data.Nsamples = H.nchirpSamples;
            # Load each chirp into a row
            apres_data.Endind = apres_data.Startind + apres_data.SamplesPerChirp - 1;

            apres_data.vif = np.zeros((apres_data.ChirpsInBurst,apres_data.SamplesPerChirp)) # preallocate array
            chirpInterval = 1.6384/(24*3600); # days
            for chirp in range(apres_data.ChirpsInBurst):
                apres_data.vif[chirp,:] = apres_data.data[apres_data.Startind[chirp]:apres_data.Endind[chirp]]
                apres_data.chirpNum[chirp,0] = chirp                                            # chirp number in burst
                apres_data.chirpAtt[chirp,0] = AttSet(1+mod(chirp-1,numel(AttSet)))             # attenuator setting for chirp
                apres_data.chirpTime[chirp,0] = apres_data.TimeStamp + chirpInterval*(chirp-1)  # time of chirp
        else:
            apres_data.er = 3.18;
            # Load each chirp into a row

            apres_data.Endind = apres_data.Startind + apres_data.SamplesPerChirp - 1
            apres_data.vif = np.zeros((apres_data.ChirpsInBurst,apres_data.SamplesPerChirp)) # preallocate array
            chirpInterval = 1.6384/(24*3600) # days

            for chirp in range(apres_data.ChirpsInBurst):
                apres_data.vif[chirp,:] = apres_data.data[apres_data.Startind[chirp]:apres_data.Endind[chirp]]
                apres_data.chirpNum[chirp,0] = chirp                                            # chirp number in burst
                apres_data.chirpAtt[chirp,0] = AttSet[1+mod(chirp-1,numel(AttSet))]             # attenuator setting for chirp
                apres_data.chirpTime[chirp,0] = apres_data.TimeStamp + chirpInterval*(chirp-1)  # time of chirp

            apres_data.ChirpsInBurst = np.shape(apres_data.vif,0)
            apres_data.SamplesPerChirp = np.shape(apres_data.vif,1)
            apres_data.dt = 1./apres_data.fs                            # sample interval (s)
            apres_data.T = (np.shape(apres_data.vif,1)-1)/apres_data.fs     # period between first and last sample
            #apres_data.T = size(apres_data.vif,1)/apres_data.fs; # period of sampling (cls test 26 aug 2014)
            # - this makes the amplitude of the fft centred at the right range, but phase wrong

            apres_data.f1 = apres_data.f0 + apres_data.T*apres_data.K/(2*np.pi) # stop frequency
            #apres_data.f1 = apres_data.f0 + apres_data.dt*(apres_data.SamplesPerChirp-1)*apres_data.K/(2*pi); # stop frequency

            #apres_data.B = apres_data.f1-apres_data.f0; # bandwidth (hz)
            #apres_data.B = apres_data.T*(apres_data.K/(2*pi)); # bandwidth (hz)
            apres_data.B = (np.shape(apres_data.vif,1)/apres_data.fs)*(apres_data.K/(2*np.pi)) # bandwidth (hz)

            apres_data.fc = np.mean([apres_data.f0, apres_data.f1]); # Centre frequency
            #apres_data.fc = apres_data.f0 + apres_data.B/2; # Centre frequency
            apres_data.ci = 3e8/np.sqrt(apres_data.er); # velocity in material
            apres_data.lambdac = apres_data.ci/apres_data.fc; # Centre wavelength

    # Create time and frequency stamp for samples
    apres_data.t = apres_data.dt*(range(np.shape(apres_data.vif,2)-1)); # sampling times (rel to first)
    apres_data.f = apres_data.f0 + apres_data.t*apres_data.K/(2*np.pi);

    # Calibrate
    #ca13 = [1 6]; # 2013
    #ca14 = [1 2]; # 2014
    #ca = [1 4];
    #apres_data = fmcw_cal(apres_data,ca13);

    apres_data.flags = ApresFlags()
    #apres_data.tnum = dzt_data.data.shape[1]
    #apres_data.trace_num = np.arange(dzt_data.data.shape[1]) + 1
    #apres_data.trig_level = np.zeros((dzt_data.tnum, ))
    #apres_data.pressure = np.zeros((dzt_data.tnum, ))
    #apres_data.dt = dzt_data.range / dzt_data.snum * 1.0e-9
    #apres_data.travel_time = np.atleast_2d(np.arange(0,
    #                                               dzt_data.range / 1.0e3,
    #                                               dzt_data.dt * 1.0e6)).transpose()
    #apres_data.travel_time += dzt_data.dt * 1.0e6
