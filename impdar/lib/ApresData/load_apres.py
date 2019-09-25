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
import datetime
from . import ApresData
from . import ApresFlags

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
        # TODO: fix this in the __init__ file
        apres_data = ApresData(fn_apres)
    else:
        apres_data = ApresData(None)
        apres_data.header.update_parameters()
        apres_data.data = load_burst_rmb(apres_data, fn_apres, burst, fs, apres_data.header.file_format)

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

        if apres_data.file_format == 5 or apres_data.file_format == 4:
            apres_data.K = apres_data.header.K
            apres_data.f0 = apres_data.header.startFreq
            apres_data.fs = apres_data.header.fs
            apres_data.f1 = apres_data.header.startFreq + apres_data.header.chirpLength * apres_data.header.K/2/np.pi
            apres_data.SamplesPerChirp = round(apres_data.header.chirpLength * apres_data.header.fs);
            apres_data.T = apres_data.header.chirpLength
            apres_data.B = apres_data.header.chirpLength * apres_data.header.K/2/np.pi
            apres_data.fc = apres_data.header.startFreq + apres_data.B/2
            apres_data.dt = 1./apres_data.header.fs
            apres_data.er = 3.18
            apres_data.ci = 3e8/np.sqrt(apres_data.er);
            apres_data.lambdac = apres_data.ci/apres_data.fc;
            apres_data.Nsamples = apres_data.header.nchirpSamples;
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

# -----------------------------------------------------------------------------------------------------

def load_burst_rmb(self,burst=1,fs=40000,max_header_len=2000,burst_pointer=0):
    """
    Load bursts from the apres acquisition.
    Normally, this should be called from the load_apres function.

    Parameters
    ---------
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

    ### Original Matlab Script Notes ###
    Read FMCW data file from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)

    Corrected so that Sampling Frequency has correct use (ie, not used in
    this case)
    """

    if self.header.fn is None:
        raise TypeError('Read in the header before loading data.')

    self.header.file_read_code = 0

    try:
        fid = open(self.header.fn,'rb')
    except:
        # Unknown file
        self.header.file_read_code = -1;
        raise TypeError('Cannot open file', self.header.fn)

    # Get the total length of the file
    fid.seek(0,2)
    file_len = fid.tell()
    burst_count = 1

    # --- Read bursts in a loop --- #
    while burst_count <= burst and burst_pointer <= file_len - max_header_len:
        # Go to burst pointer and read the header for the burst
        fid.seek(burst_pointer)

        try:
            # Read header values
            strings = ['N_ADC_SAMPLES=','NSubBursts=','Average=','nAttenuators=','Attenuator1=',
                      'AFGain=','TxAnt=','RxAnt=']
            output = np.empty((len(strings))).astype(str)
            last_read_ind = 0
            for i,string in enumerate(strings):
                if string in self.header.header_string:
                    search_start = self.header.header_string.find(string)
                    search_end = self.header.header_string[search_start:].find('\\')
                    output[i] = self.header.header_string[search_start+len(string):search_end+search_start]
                    last_read_ind = max([last_read_ind,search_end])
            # Write header values to data object
            self.snum = int(output[0])
            self.n_subbursts = int(output[1])
            self.average = int(output[2])
            self.n_attenuators = int(output[3])
            self.attenuator1 = np.array(output[4].split(',')).astype(int)[self.n_attenuators-1]
            self.attenuator2 = np.array(output[5].split(',')).astype(int)[self.n_attenuators-1]
            self.tx_ant = np.array(output[6].split(',')).astype(int)
            self.rx_ant = np.array(output[7].split(',')).astype(int)

            self.tx_ant = self.tx_ant[self.tx_ant==1]
            self.rx_ant = self.rx_ant[self.rx_ant==1]

            if self.average != 0:
                self.chirps_in_burst = 1
            else:
                self.chirps_in_burst = self.n_subbursts*len(self.tx_ant)*\
                                            len(self.rx_ant)*self.n_attenuators

            # End of burst
            search_string = '*** End Header ***'
            search_ind = self.header.header_string.find(search_string)
            burst_pointer += search_ind + len(search_string)

        except:
            # If the burst read is unsuccessful exit with an updated read code
            self.header.file_read_code = -2
            self.burst = burst_count

        words_per_burst = self.chirps_in_burst*self.snum

        # Move the burst pointer
        if burst_count < burst and burst_pointer <= file_len - max_header_len:
            if self.average != 0:
                burst_pointer += self.chirps_in_burst*self.snum*4
            else:
                burst_pointer += self.chirps_in_burst*self.snum*2

        burst_count += 1

    # --- Get remaining information from burst header --- #

    # Look for a few different strings and save output
    strings = ['Time stamp=','Temp1=','Temp2=','BatteryVoltage=']
    output = np.empty((len(strings))).astype(str)
    for i,string in enumerate(strings):
        if string in self.header.header_string:
            search_start = self.header.header_string.find(string)
            search_end = self.header.header_string[search_start:].find('\\')
            output[i] = self.header.header_string[search_start+len(string):search_end+search_start]

    if 'Time stamp' not in self.header.header_string:
        self.header.file_read_code = -4
    else:
        self.time_stamp = datetime.datetime.strptime(output[0],'%Y-%m-%d %H:%M:%S')
        timezero = datetime.datetime(1, 1, 1, 0, 0, 0)
        day_offset = self.time_stamp - timezero
        self.decday = day_offset.days + 377. # Matlab compatable

    self.temperature1 = float(output[1])
    self.temperature2 = float(output[2])
    self.battery_voltage = float(output[3])

    # --- Read in the actual data --- #

    # Go to the end of the header
    end_byte = b'*** End Header ***'
    data_ind = fid.read(5000).rfind(end_byte) + len(end_byte)
    fid.seek(data_ind)

    # Only if all the bursts were read
    if burst_count != burst+1:
        self.burst = burst_count - 1
        self.header.file_read_code = -4
    else:
        # TODO: Check the other readers for average == 1 or average == 2
        if self.average == 2:
            self.data = np.fromfile(fid,dtype='uint32',count=words_per_burst)
        elif self.average == 1:
            fid.seek(burst_pointer+1)
            self.data = np.fromfile(fid,dtype='float4',count=words_per_burst)
        else:
            self.data = np.fromfile(fid,dtype='uint16',count=words_per_burst)

        if fid.tell()-(burst_pointer-1) < words_per_burst:
            self.header.file_read_code = 2

        self.data[self.data<0] = self.data[self.data<0] + 2**16.
        self.data = self.data.astype(float) * 2.5/2**16.

        if self.average == 2:
            self.data /= (self.n_subbursts*self.n_attenuators)

        self.start_ind = np.transpose(np.arange(0,self.snum*self.chirps_in_burst,self.snum))
        self.end_ind = self.start_ind + self.snum
        self.burst = burst

    fid.close()

    # Clean temperature record (wrong data type?)
    if self.temperature1 > 300:
        self.temperature1 -= 512
    if self.temperature2 > 300:
        self.temperature2 -= 512

