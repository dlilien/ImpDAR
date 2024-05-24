#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Benjamin Hills <bhills@uw.com>
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

import os

import numpy as np
from scipy.io import loadmat

import datetime
import re
from . import ApresData
from ..ImpdarError import ImpdarError

try:
    from netCDF4 import Dataset
    nc_load = True
except:
    nc_load = False

def load_apres(fns_apres, burst=1, fs=40000, *args, **kwargs):
    """Load and concatenate all apres data from several files

    Parameters
    ----------
    fns_apres: list of file names for ApresData
        each loads object to concatenate
    burst: int

    fs: int
        sampling frequency (per second)

    Returns
    -------
    out
        A single, concatenated output.
    """

    apres_data = []
    for fn in fns_apres:
        try:
            apres_data.append(load_apres_single_file(
                fn, burst=burst, fs=fs, *args, **kwargs))
        except:
            Warning('Cannot load file: '+fn)

    from copy import deepcopy
    out = deepcopy(apres_data[0])
    ext = os.path.splitext(fns_apres[0])[1]

    if len(apres_data)>1 or ext in ['.DAT','.dat']:
        for dat in apres_data[1:]:
            if out.snum != dat.snum:
                raise ValueError('Need the same number of vertical samples in each file')
            if out.cnum != dat.cnum:
                raise ValueError('Need the same number of chirps in each file')
            if not np.all(out.travel_time == dat.travel_time):
                raise ValueError('Need matching travel time vectors')
            if not np.all(out.frequencies == dat.frequencies):
                raise ValueError('Need matching frequency vectors')

        out.data = np.vstack([[dat.data] for dat in apres_data])
        out.chirp_num = np.vstack([[dat.chirp_num] for dat in apres_data])
        out.chirp_att = np.vstack([[dat.chirp_att] for dat in apres_data])
        out.chirp_time = np.vstack([[dat.chirp_time] for dat in apres_data])
        out.decday = np.hstack([dat.decday for dat in apres_data])
        out.time_stamp = np.hstack([dat.time_stamp for dat in apres_data])
        out.lat = np.hstack([dat.lat for dat in apres_data])
        out.long = np.hstack([dat.long for dat in apres_data])
        out.temperature1 = np.hstack([dat.temperature1 for dat in apres_data])
        out.temperature2 = np.hstack([dat.temperature2 for dat in apres_data])
        out.battery_voltage = np.hstack(
            [dat.battery_voltage for dat in apres_data])
        out.bnum = np.shape(out.data)[0]

    else:
        out = deepcopy(apres_data[0])

    out.fn = os.path.splitext(fns_apres[0])[0]

    return out


def load_apres_single_file(fn_apres, burst=1, fs=40000, *args, **kwargs):
    """
    Load ApRES data
    This function calls the load_burst function below

    Parameters
    ---------
    fn_apres: string
        file name
    burst: int
        number of bursts to load
    fs: int
        sampling frequency

    Returns
    ---------
    ApresData: class
        data object


    ### Original Matlab Notes ###

    Craig Stewart
    2013 April 24
    2013 September 30 - corrected error in vif scaling
    2014/5/20 time stamp moved here from fmcw_derive_parameters (so that this
    is not overwritted later)
    2014/5/21 changed how radar chirp is defined (now using chirp gradient as
    fundamental parameter)
    2014/5/22 fixed bug in chirptime
    2014/8/21 moved make odd length to external (called from fmcw_range)
    2014/10/22 KWN - edited to allow for new headers in RMB2 files
    """

    # We can be raw, impdar mat, BAS mat, impdar h5, or BAS nc
    ext = os.path.splitext(fn_apres)[1]
    if ext == '.mat':
        dat = loadmat(fn_apres)

        # The old bas script uses vdat as the main struct
        if 'vdat' in dat:
            impdar_format = False
        else:
            impdar_format = True
        dat = None

        if impdar_format:
            apres_data = ApresData(fn_apres)
        else:
            apres_data = load_BAS_mat(fn_apres)
        return apres_data

    elif ext == '.h5':
        return ApresData(fn_apres)

    elif ext == '.nc':
        return load_BAS_nc(fn_apres)

    elif ext not in ['.dat','.DAT']:
        raise ValueError('Expecting a certain filetype; either .mat, .h5, .dat, .DAT, .nc')

    else:
        # Load data and reshape array
        apres_data = ApresData(None)
        apres_data.header.update_parameters(fn_apres)
        start_ind, end_ind = load_burst(apres_data, burst, fs)

        # Extract just good chirp data from voltage record and rearrange into
        # matrix with one chirp per row
        # note: you can't just use reshape as we are also cropping the 20K samples
        # of sync tone etc which occur after each 40K of chirp.
        AttSet = apres_data.header.attenuator1 + 1.0j * apres_data.header.attenuator2  # unique code for attenuat

        # Sampling parameters
        if apres_data.header.file_format is None:
            raise TypeError("File format is 'None', cannot load")
        else:
            apres_data.header.f1 = apres_data.header.f0 + \
                apres_data.header.chirp_length * apres_data.header.chirp_grad/2./np.pi
            apres_data.header.bandwidth = apres_data.header.chirp_length * \
                apres_data.header.chirp_grad/2/np.pi
            apres_data.header.fc = apres_data.header.f0 + apres_data.header.bandwidth/2.
            apres_data.dt = 1./apres_data.header.fs
            apres_data.header.er = 3.18
            apres_data.header.ci = 3e8/np.sqrt(apres_data.header.er)
            apres_data.header.lambdac = apres_data.header.ci/apres_data.header.fc

            # Load each chirp into a row
            # preallocate array
            data_load = np.zeros((apres_data.cnum, apres_data.snum))
            apres_data.chirp_num = np.arange(apres_data.cnum)
            apres_data.chirp_att = np.zeros((apres_data.cnum)).astype(np.cdouble)
            apres_data.chirp_time = np.zeros((apres_data.cnum))
            apres_data.header.chirp_interval = 1.6384/(24.*3600.)  # in days; apparently this is always the interval?
            for chirp in range(apres_data.cnum):
                data_load[chirp, :] = apres_data.data[start_ind[chirp]: end_ind[chirp]]
                # attenuator setting for chirp
                apres_data.chirp_att[chirp] = AttSet[chirp//apres_data.cnum]
                apres_data.chirp_time[chirp] = apres_data.decday + apres_data.header.chirp_interval*chirp
            apres_data.data = data_load

    # Create time and frequency stamp for samples
    # sampling times (rel to first)
    apres_data.travel_time = apres_data.dt * np.arange(apres_data.snum)
    apres_data.frequencies = apres_data.header.f0 + apres_data.travel_time * apres_data.header.chirp_grad/(2.*np.pi)
    apres_data.travel_time *= 1.0e6

    apres_data.data_dtype = apres_data.data.dtype

    apres_data.check_attrs()

    return apres_data


def load_burst(self, burst=1, fs=40000, max_header_len=2000, burst_pointer=0):
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

    Returns
    ---------
    start_ind: int
        data index for given burst
    end_ind: int
        data index for given burst


    ### Original Matlab Script Notes ###
    Read FMCW data file from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)

    Corrected so that Sampling Frequency has correct use (ie, not used in
    this case)
    """

    if self.header.fn is None:
        raise TypeError('Read in the header before loading data.')
    if self.header.file_format != 5:
        raise TypeError('Loading functions have only been written for rmb5 data.\
                        Look back to the original Matlab scripts if you need to implement earlier formats.')

    try:
        fid = open(self.header.fn, 'rb')
    except FileNotFoundError:
        # Unknown file
        self.flags.file_read_code = 'Unable to read file' + self.header.fn
        raise ImpdarError('Cannot open file', self.header.fn)

    # Get the total length of the file
    fid.seek(0, 2)
    file_len = fid.tell()
    burst_count = 1

    # --- Read bursts in a loop --- #
    while burst_count <= burst and burst_pointer <= file_len - max_header_len:
        # Go to burst pointer and read the header for the burst
        fid.seek(burst_pointer)
        self.header.read_header(self.header.fn, max_header_len)

        try:
            # Read header values
            strings = ['N_ADC_SAMPLES=', 'NSubBursts=', 'Average=', 'nAttenuators=', 'Attenuator1=',
                       'AFGain=', 'TxAnt=', 'RxAnt=']
            output = np.empty((len(strings))).astype(str)
            for i, string in enumerate(strings):
                if string in self.header.header_string:
                    search_start = self.header.header_string.find(string)
                    search_end = self.header.header_string[search_start:].find(
                        '\\')
                    output[i] = self.header.header_string[search_start + len(string):
                                                          search_end+search_start]

            # Write header values to data object
            self.snum = int(output[0])
            self.header.average = int(output[2])
            self.header.n_subbursts = int(output[1])
            self.header.n_attenuators = int(output[3])
            self.header.attenuator1 = np.array(output[4].split(',')).astype(int)[
                :self.header.n_attenuators]
            self.header.attenuator2 = np.array(output[5].split(',')).astype(int)[
                :self.header.n_attenuators]
            self.header.tx_ant = np.array(output[6].split(',')).astype(int)
            self.header.rx_ant = np.array(output[7].split(',')).astype(int)

            self.header.tx_ant = self.header.tx_ant[self.header.tx_ant == 1]
            self.header.rx_ant = self.header.rx_ant[self.header.rx_ant == 1]

            if self.header.average != 0:
                self.cnum = 1
            else:
                self.cnum = self.header.n_subbursts*len(self.header.tx_ant) *\
                    len(self.header.rx_ant)*self.header.n_attenuators

            # End of burst
            search_string = '*** End Header ***'
            search_ind = self.header.header_string.find(search_string)
            burst_pointer += search_ind + len(search_string)

        except ValueError:
            # If the burst read is unsuccessful exit with an updated read code
            self.flags.file_read_code = 'Corrupt header in burst' + \
                str(burst_count) + 'for file' + self.header.fn
            self.bnum = burst_count
            raise ImpdarError('Burst Read Failed.')

        # Move the burst pointer
        if burst_count < burst and burst_pointer <= file_len - max_header_len:
            if self.header.average != 0:
                burst_pointer += self.cnum*self.snum*4
            else:
                burst_pointer += self.cnum*self.snum*2

        burst_count += 1

    # --- Get remaining information from burst header --- #

    # Look for a few different strings and save output
    strings = ['Time stamp=', 'Latitude=', 'Longitude=',
               'Temp1=', 'Temp2=', 'BatteryVoltage=']
    output = []
    for i, string in enumerate(strings):
        if string in self.header.header_string:
            search_start = [m.start() for m in re.finditer(
                string, self.header.header_string)]
            search_end = [self.header.header_string[ind:].find(
                '\\') for ind in search_start]
            out = [self.header.header_string[search_start[i] + len(string):
                                             search_end[i]+search_start[i]]
                   for i in range(len(search_start))]
            output.append(out)

    if 'Time stamp' not in self.header.header_string:
        self.flags.file_read_code = 'Burst' + \
            str(self.bnum) + 'not found in file' + self.header.fn
    else:
        self.time_stamp = np.array([datetime.datetime.strptime(str_time, '%Y-%m-%d %H:%M:%S') for str_time in output[0]])
        timezero = datetime.datetime(1, 1, 1, 0, 0, 0)
        day_offset = self.time_stamp - timezero
        self.decday = np.array([offset.days + offset.seconds/86400. for offset in day_offset]) + 366. # Matlab compatable

    self.lat = np.array(output[1]).astype(float)
    self.long = np.array(output[2]).astype(float)
    self.temperature1 = np.array(output[3]).astype(float)
    self.temperature2 = np.array(output[4]).astype(float)
    self.battery_voltage = np.array(output[5]).astype(float)

    # --- Read in the actual data --- #

    # Go to the end of the header
    end_byte = b'*** End Header ***'
    data_ind = fid.read(max_header_len).rfind(end_byte) + len(end_byte)
    fid.seek(data_ind)

    # Only if all the bursts were read
    if burst_count != burst+1:
        # too few bursts in file
        self.flags.file_read_code = 'Burst' + \
            str(self.bnum) + 'not found in file' + self.header.fn
        self.bnum = burst_count - 1
        raise ImpdarError('Burst {:d} not found in file {:s}'.format(self.bnum, self.header.fn))
    else:
        if self.header.average == 2:
            self.data = np.fromfile(
                fid, dtype='uint32', count=self.cnum*self.snum)
        elif self.header.average == 1:
            fid.seek(burst_pointer+1)
            self.data = np.fromfile(
                fid, dtype='float4', count=self.cnum*self.snum)
        else:
            self.data = np.fromfile(
                fid, dtype='uint16', count=self.cnum*self.snum)

        if fid.tell()-(burst_pointer-1) < self.cnum*self.snum:
            self.flags.file_read_code = 'Corrupt header in burst' + \
                str(burst_count) + 'for file' + self.header.fn

        self.data[self.data < 0] = self.data[self.data < 0] + 2**16.
        self.data = self.data.astype(float) * 2.5/2**16.

        if self.header.average == 2:
            self.data /= (self.header.n_subbursts*self.header.n_attenuators)

        start_ind = np.transpose(np.arange(0, self.snum*self.cnum, self.snum))
        end_ind = start_ind + self.snum
        self.bnum = burst

    fid.close()

    # Clean temperature record (wrong data type?)
    self.temperature1[self.temperature1 > 300] -= 512
    self.temperature2[self.temperature2 > 300] -= 512

    self.flags.file_read_code = 'Successful Read'

    return start_ind, end_ind


def load_BAS_mat(fn, chirp_interval=1.6384/(24.*3600.)):
    """Load apres data from a mat file saved by software from the British Antarctic Survey.

    Parameters
    ----------
    fn: str
        Matlab file name for ApresData
    chirp_interval: float
        set by default

    Returns
    -------
    apres_data: class
        ImpDAR data object
    """

    mat = loadmat(fn)

    apres_data = ApresData(None)

    # Load the header first. Some of these are pre-set, some really loaded
    apres_data.header.f0 = mat['vdat'][0]['f0'][0][0][0]
    apres_data.header.fs = mat['vdat'][0]['fs'][0][0][0]
    apres_data.header.f1 = mat['vdat'][0]['f1'][0][0][0]
    apres_data.header.fc = mat['vdat'][0]['fc'][0][0][0]
    apres_data.header.attenuator1 = mat['vdat'][0]['Attenuator_1'][0][0][0]
    apres_data.header.attenuator2 = mat['vdat'][0]['Attenuator_2'][0][0][0]
    apres_data.header.chirp_length = mat['vdat'][0]['T'][0][0][0]
    apres_data.header.chirp_grad = mat['vdat'][0]['K'][0][0][0]
    apres_data.header.bandwidth = mat['vdat'][0]['B'][0][0][0]
    apres_data.header.lambdac = mat['vdat'][0]['lambdac'][0][0][0]
    apres_data.header.er = mat['vdat'][0]['er'][0][0][0]
    apres_data.header.ci = mat['vdat'][0]['ci'][0][0][0]

    apres_data.snum = mat['vdat'][0]['Nsamples'][0][0][0]
    apres_data.cnum = mat['vdat'][0]['chirpNum'][0][0][0]
    apres_data.bnum = mat['vdat'][0]['Burst'][0][0][0]
    apres_data.header.n_subbursts = mat['vdat'][0]['SubBurstsInBurst'][0][0][0]
    apres_data.header.average = mat['vdat'][0]['Average'][0][0][0]

    apres_data.travel_time = mat['vdat'][0]['t'][0].T
    apres_data.frequencies = mat['vdat'][0]['f'][0].T
    apres_data.dt = 1.0 / apres_data.header.fs

    apres_data.chirp_num = np.arange(apres_data.cnum) + 1
    apres_data.chirp_att = mat['vdat'][0]['chirpAtt'][0]
    apres_data.decday = mat['vdat'][0]['TimeStamp'][0][0][0]
    apres_data.header.chirp_interval = chirp_interval
    apres_data.chirp_time = apres_data.decday + apres_data.header.chirp_interval * np.arange(0.0, apres_data.cnum, 1.0)

    apres_data.data = mat['vdat'][0]['vif'][0]
    if len(apres_data.data.shape) == 2:
        apres_data.data = np.reshape(apres_data.data, (1, apres_data.data.shape[0], apres_data.data.shape[1]))

    apres_data.check_attrs()

    return apres_data


def load_BAS_nc(fn, fs=40000, chirp_interval=1.6384/(24.*3600.), *args, **kwargs):
    """Load apres data from a netcdf file saved by software from the British Antarctic Survey.
    https://github.com/antarctica/bas-apres

    Parameters
    ----------
    fn: str
        netCDF file name for ApresData

    Returns
    -------
    apres_data: class
        ImpDAR data object
    """

    if not nc_load:
        raise ImportError('Need the netCDF4 library to load nc files.')

    apres_data = ApresData(None)

    apres_data.bnum = int(0)

    with Dataset(fn, 'r') as fin:
        # Convert each netCDF group to an ApRESBurst object
        # each key is its own burst
        apres_data.bnum += 1
        if len(fin.groups) > 0:
            key = list(fin.groups.keys())[0]
            attrs = vars(fin.groups[key]).copy()
            apres_data.data = np.array([fin.groups[key].variables['data'][:]])
        else:
            attrs = vars(fin).copy()
            apres_data.data = np.array([fin.variables['data'][:]])

        apres_data.header.fs = fs
        apres_data.header.fn = fn
        apres_data.header.file_format = 'BAS_nc'
        apres_data.header.noDwellHigh = int(attrs['NoDwell'])
        #apres_data.header.noDwellLow = int(attrs['NoDwell'])
        apres_data.header.f0 = float(attrs['StartFreq'])
        apres_data.header.f_stop = float(attrs['StopFreq'])
        apres_data.header.ramp_up_step = float(attrs['FreqStepUp'])
        apres_data.header.ramp_down_step = float(attrs['FreqStepDn'])
        apres_data.header.tstep_up = float(attrs['TStepUp'])
        apres_data.header.tstep_down = float(attrs['TStepDn'])
        apres_data.header.ramp_dir = 'up'
        apres_data.header.nsteps_DDS = round(abs((apres_data.header.f_stop - apres_data.header.f0)/apres_data.header.ramp_up_step)) # abs as ramp could be down
        apres_data.header.chirp_length = int(apres_data.header.nsteps_DDS * apres_data.header.tstep_up)
        apres_data.header.nchirp_samples = round(apres_data.header.chirp_length * apres_data.header.fs)
        # If number of ADC samples collected is less than required to collect
        # entire chirp, set chirp length to length of series actually collected
        apres_data.header.snum = float(attrs['N_ADC_SAMPLES'])
        if apres_data.header.nchirp_samples > apres_data.header.snum:
            apres_data.header.chirp_length = apres_data.header.snum / apres_data.header.fs
        apres_data.header.chirp_grad = 2.*np.pi*(apres_data.header.ramp_up_step/apres_data.header.tstep_up) # chirp gradient (rad/s/s)
        if apres_data.header.f_stop > 400e6:
            apres_data.header.ramp_dir = 'down'
        else:
            apres_data.header.ramp_dir = 'up'

        apres_data.header.f1 = apres_data.header.f0 + \
            apres_data.header.chirp_length * apres_data.header.chirp_grad/2./np.pi
        apres_data.header.bandwidth = apres_data.header.chirp_length * \
            apres_data.header.chirp_grad/2/np.pi
        apres_data.header.fc = apres_data.header.f0 + apres_data.header.bandwidth/2.
        apres_data.header.er = 3.18
        apres_data.header.ci = 3e8/np.sqrt(apres_data.header.er)
        apres_data.header.lambdac = apres_data.header.ci/apres_data.header.fc
        apres_data.header.n_attenuators = int(attrs['nAttenuators'])
        str_attenuators = attrs['Attenuator1']
        apres_data.header.attenuator1 = np.array(str_attenuators.split(',')).astype(int)[
                                                    :apres_data.header.n_attenuators]
        AFGain = attrs['AFGain']
        apres_data.header.attenuator2 = np.array(AFGain.split(',')).astype(int)[
                                                    :apres_data.header.n_attenuators]
        apres_data.header.tx_ant = attrs['TxAnt']
        apres_data.header.rx_ant = attrs['RxAnt']

        # Should be the same for every burst
        apres_data.header.average = float(attrs['Average'])
        apres_data.header.chirp_interval = chirp_interval
        apres_data.dt = 1.0 / apres_data.header.fs
        apres_data.snum = int(attrs['N_ADC_SAMPLES'])
        apres_data.cnum = int(attrs['NSubBursts'])
        apres_data.header.n_subbursts = int(attrs['NSubBursts'])

        # Should be appended to every burst
        apres_data.temperature1 = np.array([float(attrs['Temp1'])])
        apres_data.temperature2 = np.array([float(attrs['Temp2'])])
        apres_data.battery_voltage = np.array([float(attrs['BatteryVoltage'])])

        time_stamp = attrs['Time stamp']
        apres_data.time_stamp = np.array([datetime.datetime.strptime(time_stamp, '%Y-%m-%d %H:%M:%S')])
        timezero = datetime.datetime(1, 1, 1, 0, 0, 0)
        day_offset = apres_data.time_stamp - timezero
        apres_data.decday = np.array([offset.days + offset.seconds/86400. for offset in day_offset]) + 366. # Matlab compatable
        apres_data.chirp_time = apres_data.decday + apres_data.header.chirp_interval * np.arange(0.0, apres_data.cnum, 1.0)

        AttSet = apres_data.header.attenuator1 + 1j * \
                apres_data.header.attenuator2  # unique code for attenuator setting
        apres_data.chirp_att = np.zeros((apres_data.cnum)).astype(np.cdouble)
        for chirp in range(apres_data.cnum):
            apres_data.chirp_att[chirp] = AttSet[chirp//apres_data.cnum]
        apres_data.chirp_num = np.array([np.arange(apres_data.cnum) + 1])

        if len(fin.groups) > 1:
            for key in fin.groups[1:]:
                apres_data.data = np.append(apres_data.data,
                                            [fin.groups[key].variables['data'][:]],
                                            axis=0)

                # Append to these data arrays for subsequent bursts
                apres_data.temperature1 = np.append(apres_data.temperature1, float(attrs['Temp1']))
                apres_data.temperature2 = np.append(apres_data.temperature2, float(attrs['Temp2']))
                apres_data.battery_voltage = np.append(apres_data.battery_voltage, float(attrs['BatteryVoltage']))

                time_stamp = attrs['Time stamp']
                ts = datetime.datetime.strptime(time_stamp, '%Y-%m-%d %H:%M:%S')
                apres_data.time_stamp = np.append(apres_data.time_stamp,ts)
                timezero = datetime.datetime(1, 1, 1, 0, 0, 0)
                day_offset = ts - timezero
                decday = day_offset.days + day_offset.seconds/86400. + 366.
                apres_data.decday = np.append(apres_data.decday, decday)
                apres_data.chirp_time = np.append(apres_data.chirp_time,
                                                  decday + apres_data.header.chirp_interval * np.arange(0.0, apres_data.cnum, 1.0),
                                                  axis=0)

                apres_data.chirp_num = np.append(apres_data.chirp_num,
                                                 [np.arange(apres_data.cnum) + 1],
                                                 axis=0)

    apres_data.travel_time = apres_data.dt * np.arange(apres_data.snum)
    apres_data.frequencies = apres_data.header.f0 + apres_data.travel_time * apres_data.header.chirp_grad/(2.*np.pi)
    apres_data.travel_time *= 1.0e6

    apres_data.data_dtype = apres_data.data.dtype

    apres_data.check_attrs()

    return apres_data
