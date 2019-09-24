#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Load ApRES data

This code is based on a series of Matlab scripts
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
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags
from scipy.io import loadmat

# --------------------------------------------------------------------------------------------

def file_format(fn,max_header_len=2000):
    """
    Determine fmcw file format from burst header using keyword presence

    Craig Stewart
    2013-10-20
    Updated by Keith Nicholls, 2014-10-22: RMB2
    """

    with open(fn,'rb') as fid:
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

# --------------------------------------------------------------------------------------------

class apres_parameters:
    def __init__(self):
        """Initialize data paramaters"""
        self.fsysclk = 1e9
        self.fs = 4e4

    def update_parameters(fn_apres):
        fid = open(fn,'rb')
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
                val = header[loc2[k]+2:loc2[k]+loc3]
                print(int(val, 16))
                #val = dec2bin(hex2dec(val)); val = fliplr(val);
                #H.noDwellHigh = str2num(val(18+1));
                #H.noDwellLow = str2num(val(17+1));

            #elif case == 'Reg08':
            #    # Phase offset word Register (POW) Address 0x08. 2 Bytes dTheta = 360*POW/2^16.
            #    val = char(reg{1,2}(k));
            #    H.phaseOffsetDeg = hex2dec(val(1:4))*360/2^16;

            elif case == 'Reg0B=':
                # Digital Ramp Limit Register Address 0x0B
                # Digital ramp upper limit 32-bit digital ramp upper limit value.
                # Digital ramp lower limit 32-bit digital ramp lower limit value.
                loc3 = header[loc2[k]+2:].find('"')
                val = header[loc2[k]+2:loc2[k]+loc3]
                self.startFreq = hex2dec(val(9:end))*self.fsysclk/(2**32);
                self.stopFreq = hex2dec(val(1:8))*self.fsysclk/(2**32);

            elif case == 'Reg0C=':
                # Digital Ramp Step Size Register Address 0x0C
                # Digital ramp decrement step size 32-bit digital ramp decrement step size value.
                # Digital ramp increment step size 32-bit digital ramp increment step size value.
                loc3 = header[loc2[k]+2:].find('"')
                val = header[loc2[k]+2:loc2[k]+loc3]
                self.rampUpStep = hex2dec(val(9:end))*self.fsysclk/2^32;
                self.rampDownStep = hex2dec(val(1:8))*self.fsysclk/2^32;

            elif case == 'Reg0D=':
                # Digital Ramp Rate Register Address 0x0D
                # Digital ramp negative slope rate 16-bit digital ramp negative slope value that defines the time interval between decrement values.
                # Digital ramp positive slope rate 16-bit digital ramp positive slope value that defines the time interval between increment values.
                loc3 = header[loc2[k]+2:].find('"')
                val = header[loc2[k]+2:loc2[k]+loc3]
                self.tstepUp = hex2dec(val(5:end))*4/self.fsysclk;
                self.tstepDown = hex2dec(val(1:4))*4/self.fsysclk;

        loc = strfind(A,'SamplingFreqMode=');
        searchCR = strfind(A(loc(1):end),[char(10)]);
        self.fs = sscanf(A(loc(1)+length(['SamplingFreqMode=']):searchCR(1)+loc(1)),'%d\n');
        if self.fs == 1         # if self.fs > 70e3:
            self.fs = 8e4       #     self.fs = 80e3
        else                    # else
            self.fs = 4e4       #     self.fs = 40e3

        loc = strfind(A,'N_ADC_SAMPLES=')
        searchCR = strfind(A(loc(1):end),[char(10)])
        self.Nsamples = sscanf(A(loc(1)+length(['N_ADC_SAMPLES=']):searchCR(1)+loc(1)),'%d\n')

        self.nstepsDDS = round(abs((H.stopFreq - H.startFreq)/self.rampUpStep)) # abs as ramp could be down
        self.chirpLength = self.nstepsDDS * self.tstepUp
        self.nchirpSamples = round(self.chirpLength * self.fs)

        # If number of ADC samples collected is less than required to collect
        # entire chirp, set chirp length to length of series actually collected
        if self.nchirpSamples > self.Nsamples
            self.chirpLength = self.Nsamples / self.fs

        self.K = 2.*np.pi*(self.rampUpStep/self.tstepUp) # chirp gradient (rad/s/s)
        if(self.stopFreq > 400e6)
            self.rampDir = 'down'
        else
            self.rampDir = 'up'

        if self.noDwellHigh and self.noDwellLow:
            self.rampDir = 'upDown'
            self.nchirpsPerPeriod = np.nan # self.nchirpSamples/(self.chirpLength)

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

    apres_data = RadarData(None)
    H = apres_paramaters()

    ## Load data and reshape array
    if fn_apres[-4:] == '.mat':
        apres_data.file_format = 'mat'
        apres_mat = loadmat(fn_apres) # note when you load a mat file you get whatever burst was stored in this - not the one you selected
        apres_data.data = apres_mat['data']
    else:
        apres_data.file_format = file_format(fn_apres)
        apres_data.data = load_burst_rmb(fn_apres, burst, fs, apres_data.file_format)

    """

    # Extract just good chirp data from voltage record and rearrange into
    # matrix with one chirp per row
    # note: you can't just use reshape as we are also cropping the 20K samples
    # of sync tone etc which occur after each 40K of chirp.
    AttSet = apres_data.Attenuator_1 + 1i*apres_data.Attenuator_2; # unique code for attenuator setting


    ## Add metadata to structure


    # Sampling parameters
    apres_data.filename = filename;
    if ~ischar(FileFormat)
        apres_data.SamplesPerChirp = apres_data.Nsamples;
        apres_data.fs = 4e4; # sampling frequency
        apres_data.f0 = 2e8; # start frequency
        #apres_data.fc = 3e8; # start frequency
        apres_data.K = 2*pi*2e8; # chirp gradient in rad/s/s (200MHz/s)
        #apres_data.f0 = apres_data.f0 + (apres_data.K/(4*pi))/apres_data.fs; # start frequency
        apres_data.processing = {};

        H.update_parameters(apres_data.filename, file_fmt);

        if FileFormat == 5 || FileFormat == 4
            apres_data.K = H.K;
            apres_data.f0 = H.startFreq;
            apres_data.fs = H.fs;
            apres_data.f1 = H.startFreq + H.chirpLength * H.K/2/pi;
            apres_data.SamplesPerChirp = round(H.chirpLength * H.fs);
            apres_data.T = H.chirpLength;
            apres_data.B = H.chirpLength * H.K/2/pi;
            apres_data.fc = H.startFreq + apres_data.B/2;
            apres_data.dt = 1/H.fs;
            apres_data.er = 3.18;
            apres_data.ci = 3e8/sqrt(apres_data.er);
            apres_data.lambdac = apres_data.ci/apres_data.fc;
            apres_data.Nsamples = H.nchirpSamples;
            # Load each chirp into a row
            apres_data.Endind = apres_data.Startind + apres_data.SamplesPerChirp - 1;

            apres_data.vif = zeros(apres_data.ChirpsInBurst,apres_data.SamplesPerChirp); # preallocate array
            chirpInterval = 1.6384/(24*3600); # days
            for chirp = 1:apres_data.ChirpsInBurst
                apres_data.vif(chirp,:) = apres_data.v(apres_data.Startind(chirp):apres_data.Endind(chirp));
                apres_data.chirpNum(chirp,1) = chirp; # chirp number in burst
                apres_data.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); # attenuator setting for chirp
                apres_data.chirpTime(chirp,1) = apres_data.TimeStamp + chirpInterval*(chirp-1); # time of chirp
            end
        else
            apres_data.er = 3.18;
            # Load each chirp into a row

            apres_data.Endind = apres_data.Startind + apres_data.SamplesPerChirp - 1;
            apres_data.vif = zeros(apres_data.ChirpsInBurst,apres_data.SamplesPerChirp); # preallocate array
            chirpInterval = 1.6384/(24*3600); # days
            for chirp = 1:apres_data.ChirpsInBurst
                apres_data.vif(chirp,:) = apres_data.v(apres_data.Startind(chirp):apres_data.Endind(chirp));
                apres_data.chirpNum(chirp,1) = chirp; # chirp number in burst
                apres_data.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); # attenuator setting for chirp
                apres_data.chirpTime(chirp,1) = apres_data.TimeStamp + chirpInterval*(chirp-1); # time of chirp
            end
            apres_data.ChirpsInBurst = size(apres_data.vif,1);
            apres_data.SamplesPerChirp = size(apres_data.vif,2);
            apres_data.dt = 1/apres_data.fs; # sample interval (s)
            apres_data.T = (size(apres_data.vif,2)-1)/apres_data.fs; # period between first and last sample
            #apres_data.T = size(apres_data.vif,2)/apres_data.fs; # period of sampling (cls test 26 aug 2014)
            # - this makes the amplitude of the fft centred at the right range, but phase wrong

            apres_data.f1 = apres_data.f0 + apres_data.T*apres_data.K/(2*pi); # stop frequency
            #apres_data.f1 = apres_data.f0 + apres_data.dt*(apres_data.SamplesPerChirp-1)*apres_data.K/(2*pi); # stop frequency

            #apres_data.B = apres_data.f1-apres_data.f0; # bandwidth (hz)
            #apres_data.B = apres_data.T*(apres_data.K/(2*pi)); # bandwidth (hz)
            apres_data.B = (size(apres_data.vif,2)/apres_data.fs)*(apres_data.K/(2*pi)); # bandwidth (hz)

            apres_data.fc = mean([apres_data.f0 apres_data.f1]); # Centre frequency
            #apres_data.fc = apres_data.f0 + apres_data.B/2; # Centre frequency
            apres_data.ci = 3e8/sqrt(apres_data.er); # velocity in material
            apres_data.lambdac = apres_data.ci/apres_data.fc; # Centre wavelength

        end
    else
        apres_data.er = 3.18;
        apres_data.dt = 1/apres_data.fs;
        apres_data.ci = 3e8/sqrt(apres_data.er);
        apres_data.lambdac = apres_data.ci/apres_data.fc;
        # Load each chirp into a row

        apres_data.vif = zeros(apres_data.ChirpsInBurst,apres_data.SamplesPerChirp); # preallocate array
        chirpInterval = 1.6384/(24*3600); # days
        apres_data.Endind = apres_data.Startind + apres_data.SamplesPerChirp - 1;
        for chirp = 1:apres_data.ChirpsInBurst
            apres_data.vif(chirp,:) = apres_data.v(apres_data.Startind(chirp):apres_data.Endind(chirp));
            apres_data.chirpNum(chirp,1) = chirp; # chirp number in burst
            apres_data.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); # attenuator setting for chirp
            apres_data.chirpTime(chirp,1) = apres_data.TimeStamp + chirpInterval*(chirp-1); # time of chirp
        end
    end

    # Create time and frequency stamp for samples
    apres_data.t = apres_data.dt*(0:size(apres_data.vif,2)-1); # sampling times (rel to first)
    apres_data.f = apres_data.f0 + apres_data.t.*apres_data.K/(2*pi);

    # Calibrate
    #ca13 = [1 6]; # 2013
    #ca14 = [1 2]; # 2014
    #ca = [1 4];
    #apres_data = fmcw_cal(apres_data,ca13);

    """

# --------------------------------------------------------------------------------------------

def load_burst_rmb(fn,burst,fs,file_format,max_header_len=1500,burst_pointer=0):
    """
    #
    # Read FMCW data file from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)

    # Corrected so that Sampling Frequency has correct use (ie, not used in
    # this case)

    """

    data_code = 0

    try:
        fid = open(fn_apres,'rb')
    except:
        # Unknown file
        data_code = -1;
    fid.seek(0,2)
    file_len = fid.tell()
    burst_count = 1

    while burst_count <= burst and burst_pointer <= file_len - max_header_len:
        fid.seek(burst_pointer)
        header = str(fid.read(max_header_len))

        try:
            strings = ['N_ADC_SAMPLES=','NSubBursts=','Average=','nAttenuators=','Attenuator1=',
                      'AFGain=','TxAnt=','RxAnt=']
            output = np.empty((len(strings))).astype(str)

            for i,string in enumerate(strings):
                if string in header:
                    search_start = header.find(string)
                    search_end = header[search_start:].find('\\')
                    output[i] = header[search_start+len(string):search_end+search_start]
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
                apres_data.chips_in_burst = 1
            else:
                apres_data.chips_in_burst = apres_data.n_subbursts*len(apres_data.tx_ant)*\
                                            len(apres_data.rx_ant)*apres_data.n_attenuators

            search_string = '*** End Header ***'
            search_ind = header.find(search_string)
            burst_pointer += search_ind + len(search_string)

        except:
            data_code = -2
            apres_data.burst = burst_count

        words_per_burst = apres_data.chips_in_burst*apres_data.snum

        if burst_count < burst and burst_pointer <= file_length - max_header_len:
            if apres_data.average != 0:
                burst_pointer += apres_data.chirps_in_burst*apres_data.snum*4
            else:
                burst_pointer += apres_data.chirps_in_burst*apres_data.snum*2

        burst_count += 1

    # -----------------------------------------------------------------------------------------------

    strings = ['Time stamp=','Temp1=','Temp2=','BatteryVoltage=']
    output = np.empty((len(strings))).astype(str)

    for i,string in enumerate(strings):
        if string in header:
            search_start = header.find(string)
            search_end = header[search_start:].find('\\')
            output[i] = header[search_start+len(string):search_end+search_start]

    if 'Time stamp' not in header:
        data_code = -4
    else:
        apres_data.time_stamp = datetime.datetime.strptime(output[0],'%Y-%m-%d %H:%M:%S')

    apres_data.temperature1 = float(output[1])
    apres_data.temperature2 = float(output[2])
    apres_data.battery_voltage = float(output[3])

    # -----------------------------------------------------------------------------------------------

    if burst_count == burst+1:
        if apres_data.average == 2:
            apres_data.data = fid.read(words_per_burst)
        elif apres_data.average == 1:
            fid.seek(burst_pointer+1)
            apres_data.data = fid.read(words_per_burst)
        else:
            #apres_data.data = np.fromfile(fid,dtype='uint16',count=words_per_burst)
            apres_data.data = fid.read(16)
        if fid.tell()-(burst_pointer-1) < words_per_burst:
            data_code = 2

        apres_data.data[apres_data.data<0] = apres_data[apres_data.data<0] + 2**16.
        apres_data.data *= 2.5/2**16.

        if apres_data.average == 2:
            apres_data.data /= (apres_data.n_subbursts*apres_data.n_attenuators)

        apres_data.start_ind = np.transpose(np.arange(0,w_per_chirp_cycle*apres_data.chirps_in_burst,w_per_chirp_cycle))
        apres_data.end_ind = apres_data.start_ind + w_per_chirp_cycle-1
        apres_data.burst = burst

    else:
        apres_data.burst = burst_count -1
        data_code = -4

    fid.close()

    # Clean temperature record (wrong data type?)
    bti1 = find(apres_data.Temperature_1>300);
    if ~isempty(bti1):
        apres_data.Temperature_1(bti1) = apres_data.Temperature_1(bti1)-512;
    bti2 = find(apres_data.Temperature_2>300);
    apres_data.Temperature_2(bti2) = apres_data.Temperature_2(bti2)-512;


