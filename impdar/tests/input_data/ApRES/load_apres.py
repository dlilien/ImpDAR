#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Load ApRES data

uch of this code is either directly referencing or written from older scripts in SeisUnix:
https:ithub.com/JohnWStockwellJr/SeisUnix/wiki
Options are:
    Kirchhoff (diffraction summation)
    Stolt (frequency wavenumber, constant velocity)
    Gazdag (phase shift, either constant or depth-varying velocity)
    SeisUnix (reference su routines directly)

Author:
Benjamin Hills
bhills@uw.edu
University of Washington
Earth and Space Sciences

ar 12 2019

"""

import numpy as np
from ..RadarData import RadarData
from ..RadarFlags import RadarFlags


def load_apres(fn_apres,burst=1,fs=40000):
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
    [path,~,ext] = fileparts(filename);
    if strcmp(ext,'.mat')
        load(filename); # note when you load a mat file you get whatever burst was stored in this - not the one you selected
        FileFormat = 'mat';
    else
        FileFormat = fmcw_file_format(filename);
        vdat = LoadBurstRMB5(filename, burst, SamplingFrequency);# Data from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)
    end
    vdat.FileFormat = FileFormat;

    # Extract just good chirp data from voltage record and rearrange into
    # matrix with one chirp per row
    # note: you can't just use reshape as we are also cropping the 20K samples
    # of sync tone etc which occur after each 40K of chirp.
    AttSet = vdat.Attenuator_1 + 1i*vdat.Attenuator_2; # unique code for attenuator setting


    ## Add metadata to structure


    # Sampling parameters
    vdat.filename = filename;
    if ~ischar(FileFormat)
        vdat.SamplesPerChirp = vdat.Nsamples;
        vdat.fs = 4e4; # sampling frequency
        vdat.f0 = 2e8; # start frequency
        #vdat.fc = 3e8; # start frequency
        vdat.K = 2*pi*2e8; # chirp gradient in rad/s/s (200MHz/s)
        #vdat.f0 = vdat.f0 + (vdat.K/(4*pi))/vdat.fs; # start frequency
        vdat.processing = {};

        if FileFormat == 5
            H = fmcw_ParametersRMB2(vdat.filename);
        elseif FileFormat == 4
            H = fmcw_ParametersRMB1b(vdat.filename);
        end
        if FileFormat == 5 || FileFormat == 4
            vdat.K = H.K;
            vdat.f0 = H.startFreq;
            vdat.fs = H.fs;
            vdat.f1 = H.startFreq + H.chirpLength * H.K/2/pi;
            vdat.SamplesPerChirp = round(H.chirpLength * H.fs);
            vdat.T = H.chirpLength;
            vdat.B = H.chirpLength * H.K/2/pi;
            vdat.fc = H.startFreq + vdat.B/2;
            vdat.dt = 1/H.fs;
            vdat.er = 3.18;
            vdat.ci = 3e8/sqrt(vdat.er);
            vdat.lambdac = vdat.ci/vdat.fc;
            vdat.Nsamples = H.nchirpSamples;
            # Load each chirp into a row
            vdat.Endind = vdat.Startind + vdat.SamplesPerChirp - 1;

            vdat.vif = zeros(vdat.ChirpsInBurst,vdat.SamplesPerChirp); # preallocate array
            chirpInterval = 1.6384/(24*3600); # days
            for chirp = 1:vdat.ChirpsInBurst
                vdat.vif(chirp,:) = vdat.v(vdat.Startind(chirp):vdat.Endind(chirp));
                vdat.chirpNum(chirp,1) = chirp; # chirp number in burst
                vdat.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); # attenuator setting for chirp
                vdat.chirpTime(chirp,1) = vdat.TimeStamp + chirpInterval*(chirp-1); # time of chirp
            end
        else
            vdat.er = 3.18;
            # Load each chirp into a row

            vdat.Endind = vdat.Startind + vdat.SamplesPerChirp - 1;
            vdat.vif = zeros(vdat.ChirpsInBurst,vdat.SamplesPerChirp); # preallocate array
            chirpInterval = 1.6384/(24*3600); # days
            for chirp = 1:vdat.ChirpsInBurst
                vdat.vif(chirp,:) = vdat.v(vdat.Startind(chirp):vdat.Endind(chirp));
                vdat.chirpNum(chirp,1) = chirp; # chirp number in burst
                vdat.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); # attenuator setting for chirp
                vdat.chirpTime(chirp,1) = vdat.TimeStamp + chirpInterval*(chirp-1); # time of chirp
            end
            vdat.ChirpsInBurst = size(vdat.vif,1);
            vdat.SamplesPerChirp = size(vdat.vif,2);
            vdat.dt = 1/vdat.fs; # sample interval (s)
            vdat.T = (size(vdat.vif,2)-1)/vdat.fs; # period between first and last sample
            #vdat.T = size(vdat.vif,2)/vdat.fs; # period of sampling (cls test 26 aug 2014)
            # - this makes the amplitude of the fft centred at the right range, but phase wrong

            vdat.f1 = vdat.f0 + vdat.T*vdat.K/(2*pi); # stop frequency
            #vdat.f1 = vdat.f0 + vdat.dt*(vdat.SamplesPerChirp-1)*vdat.K/(2*pi); # stop frequency

            #vdat.B = vdat.f1-vdat.f0; # bandwidth (hz)
            #vdat.B = vdat.T*(vdat.K/(2*pi)); # bandwidth (hz)
            vdat.B = (size(vdat.vif,2)/vdat.fs)*(vdat.K/(2*pi)); # bandwidth (hz)

            vdat.fc = mean([vdat.f0 vdat.f1]); # Centre frequency
            #vdat.fc = vdat.f0 + vdat.B/2; # Centre frequency
            vdat.ci = 3e8/sqrt(vdat.er); # velocity in material
            vdat.lambdac = vdat.ci/vdat.fc; # Centre wavelength

        end
    else
        vdat.er = 3.18;
        vdat.dt = 1/vdat.fs;
        vdat.ci = 3e8/sqrt(vdat.er);
        vdat.lambdac = vdat.ci/vdat.fc;
        # Load each chirp into a row

        vdat.vif = zeros(vdat.ChirpsInBurst,vdat.SamplesPerChirp); # preallocate array
        chirpInterval = 1.6384/(24*3600); # days
        vdat.Endind = vdat.Startind + vdat.SamplesPerChirp - 1;
        for chirp = 1:vdat.ChirpsInBurst
            vdat.vif(chirp,:) = vdat.v(vdat.Startind(chirp):vdat.Endind(chirp));
            vdat.chirpNum(chirp,1) = chirp; # chirp number in burst
            vdat.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); # attenuator setting for chirp
            vdat.chirpTime(chirp,1) = vdat.TimeStamp + chirpInterval*(chirp-1); # time of chirp
        end
    end



    # Create time and frequency stamp for samples
    vdat.t = vdat.dt*(0:size(vdat.vif,2)-1); # sampling times (rel to first)
    vdat.f = vdat.f0 + vdat.t.*vdat.K/(2*pi);

    # Calibrate
    #ca13 = [1 6]; # 2013
    #ca14 = [1 2]; # 2014
    #ca = [1 4];
    #vdat = fmcw_cal(vdat,ca13);

# --------------------------------------------------------------------------------------------

def load_burst_rmb(fn,burst,fs,version=5):
    """
    # vdat = LoadBurstRMB5(Filename, Burst, SamplesPerChirp)
    #
    # Read FMCW data file from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)

    # Corrected so that Sampling Frequency has correct use (ie, not used in
    # this case)

    """
    MaxHeaderLen = 1500;
    burstpointer = 0;
    vdat.Code = 0;
    fid = fopen(Filename,'r');
    if fid >= 0
        fseek(fid,0,'eof');
        filelength = ftell(fid);
        BurstCount = 1;
        while BurstCount <= Burst && burstpointer <= filelength - MaxHeaderLen
            fseek(fid,burstpointer,'bof');
            A = fread(fid,MaxHeaderLen,'*char');
            A = A';       #'
            SearchString = 'N_ADC_SAMPLES=';
            searchind = strfind(A,SearchString);
            if ~isempty(searchind)
                try
                    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                    vdat.Nsamples = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#d');
                    WperChirpCycle = vdat.Nsamples;
                    SearchString = 'NSubBursts=';
                    searchind = strfind(A,SearchString);
                    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                    vdat.SubBurstsInBurst = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#d');

                    SearchString = 'Average=';
                    searchind = strfind(A, SearchString);
                    if isempty(searchind)
                        vdat.Average = 0; #cls 9/jan/14 -average not included in mooring deploy
                    else
                        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                        vdat.Average = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#d');
                    end

                    SearchString = 'nAttenuators=';
                    searchind = strfind(A, SearchString);
                    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                    vdat.NAttenuators = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#d',1);

                    SearchString = 'Attenuator1=';
                    searchind = strfind(A, SearchString);
                    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                    vdat.Attenuator_1 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#f',vdat.NAttenuators);

                    SearchString = 'AFGain=';
                    searchind = strfind(A, SearchString);
                    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                    vdat.Attenuator_2 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#f',vdat.NAttenuators);

                    SearchString = 'TxAnt=';
                    searchind = strfind(A, SearchString);
                    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                    vdat.TxAnt = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#d',8);

                    SearchString = 'RxAnt=';
                    searchind = strfind(A, SearchString);
                    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                    vdat.RxAnt = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#d',8);

                    ind = find(vdat.TxAnt~=1);
                    vdat.TxAnt(ind) = [];
                    ind = find(vdat.RxAnt~=1);
                    vdat.RxAnt(ind) = [];

                    if vdat.Average
                        vdat.ChirpsInBurst = 1;
                    else
                        vdat.ChirpsInBurst = vdat.SubBurstsInBurst * length(vdat.TxAnt) * ...
                           length(vdat.RxAnt) * vdat.NAttenuators;
                    end

                    SearchString = '*** End Header ***';
                    searchind = strfind(A, SearchString);

                    burstpointer = burstpointer + searchind(1) + length(SearchString);
                catch
                    vdat.Code = -2;
                    vdat.Burst = BurstCount;
                    keyboard
                    return
                end
            end
            WordsPerBurst = vdat.ChirpsInBurst * WperChirpCycle;
            if BurstCount < Burst && burstpointer <= filelength - MaxHeaderLen
                if vdat.Average
                    burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle*4;
                else
                    burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle*2;
                end
            end
            BurstCount = BurstCount + 1;
        end

        # Extract remaining information from header
        SearchString = 'Time stamp=';
        searchind = strfind(A, SearchString);
        if isempty(searchind)
            vdat.Code = -4;
            return
        end
        try
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            td = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),...
                '#d-#d-#d #d:#d:#d');
            vdat.TimeStamp = datenum(td(1),td(2),td(3),td(4),td(5),td(6));
        catch err
            vdat.Code = 1;
        end

        SearchString = 'Temp1=';
        searchind = strfind(A, SearchString);
        try
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.Temperature_1 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#f');
        catch err
            vdat.Code = 1;
        end

        SearchString = 'Temp2=';
        searchind = strfind(A, SearchString);
        try
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.Temperature_2 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#f');
        catch err
            vdat.Code = 1;
        end

        SearchString = 'BatteryVoltage=';
        searchind = strfind(A, SearchString);
        try
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.BatteryVoltage = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'#f');
        catch err
            vdat.Code = 1;
        end


        fseek(fid,burstpointer-1,'bof');
        if BurstCount == Burst+1
            if vdat.Average == 2
                [vdat.v count] = fread(fid,WordsPerBurst,'*uint32','ieee-le');
            elseif vdat.Average == 1
                fseek(fid,burstpointer+1,'bof');
                [vdat.v count] = fread(fid,WordsPerBurst,'*real*4','ieee-le');
            else
                [vdat.v count] = fread(fid,WordsPerBurst,'*uint16','ieee-le');
            end
            if count < WordsPerBurst
                vdat.Code = 2;
            end
            vdat.v(vdat.v<0) = vdat.v(vdat.v<0) + 2^16;
            vdat.v = single(vdat.v);
            vdat.v = vdat.v * 2.5 / 2^16;
            if vdat.Average == 2
                vdat.v = vdat.v / (vdat.SubBurstsInBurst * vdat.NAttenuators);
            end
            vdat.Startind = (1:WperChirpCycle:WperChirpCycle*vdat.ChirpsInBurst)';     #'
            vdat.Endind = vdat.Startind + WperChirpCycle - 1;
            vdat.Burst = Burst;
        else
            # Too few bursts in file
            vdat.Burst = BurstCount - 1;
            vdat.Code = -4;
            #keyboard
        end
        fclose(fid);
    else
        # Unknown file
        vdat.Code = -1;
    end

    # Clean temperature record (wrong data type?)
    bti1 = find(vdat.Temperature_1>300);
    if ~isempty(bti1)
        vdat.Temperature_1(bti1) = vdat.Temperature_1(bti1)-512;
    end
    bti2 = find(vdat.Temperature_2>300);
    vdat.Temperature_2(bti2) = vdat.Temperature_2(bti2)-512;

# --------------------------------------------------------------------------------------------

def file_format(fn):
    """
    # Determine fmcw file format from burst header using keyword presence

    # Craig Stewart
    # 2013-10-20
    #
    # Updated by Keith Nicholls, 2014-10-22: RMB2
    """

    [fid, msg] = fopen(filename,'rt');
    if fid == -1
        error(msg)
    end
    MaxHeaderLen = 2000;
    A = fread(fid,MaxHeaderLen,'*char');
    fclose(fid);
    A = A';      #'
    if ~isempty(strfind(A, 'SW_Issue=')); # Data from RMB2 after Oct 2014
        fmt = 5;
    elseif ~isempty(strfind(A, 'SubBursts in burst:')); # Data from after Oct 2013
        fmt = 4;
    elseif ~isempty(strfind(A, '*** Burst Header ***')); # Data from Jan 2013
        fmt = 3;
    elseif ~isempty(strfind(A, 'RADAR TIME')); # Data from Prototype FMCW radar (nov 2012)
        fmt = 2;
    else
        #fmt = 0; # unknown file format
        error('Unknown file format - check file')
    end
