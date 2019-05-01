#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import io
import sys
import struct
import datetime
import numpy as np

from .RadarData import RadarData


class Olaf(RadarData):

    def __init__(self, fns, Channel_Num=1):
        # We want to be able to use this step concatenate a series of files numbered by the controller
        if type(fns) == str:
            fns = [fns]

        sinfo = []
        s = []
        for i, fn in enumerate(fns):
            # We are going to follow the general format that was used by storead_script_v36
            with io.open(fn, 'rb') as fid:
                lines = fid.read()

            # Header information
            sinfo.append(SInfo(lines))

            # Data is stored per trace. Make a container
            s_i = [ChannelData(lines, sinfo[i]) for j in range(sinfo[i].nChannels)]

            # Read trace-by-trace, channel-by-channel
            for nTrc in range(sinfo[i].tnum):
                for nChan, si in enumerate(s_i):
                    si.read_trace(lines, sinfo[i], nTrc)

            s.append(s_i[Channel_Num - 1])

        # I don't know if we actually want to do this, but the filenaming scheme is wacky and this
        # will make any logical collection look good
        # s.sort(key=lambda x: x.Time[0])

        # Now merge the data into the normal format
        self.dt = 1. / sinfo[0].SampFreq
        self.PreTriggerDepth = sinfo[0].PreTriggerDepth
        self.FileNames = sinfo[0].fn
        self.ant_sep = sinfo[0].AntennaSeparation
        self.freq = sinfo[0].NominalFrequency
        self.travel_time = s[0].TravelTime * 1.0e6
        self.trig_level = s[0].TriggerLevel

        self.fnames = [si.fn for si in sinfo]

        # Data and things we derive from it
        self.data = np.hstack([s_i.Data for s_i in s])
        self.snum = self.data.shape[0]
        self.tnum = self.data.shape[1]
        self.trace_num = np.arange(self.tnum) + 1

        # Other variables that need concatenating
        self.decday = np.hstack([s_i.Time for s_i in s])
        self.elev = np.hstack([s_i.Altitude for s_i in s])
        self.lat = np.hstack([s_i.lat for s_i in s])
        self.long = np.hstack([s_i.long for s_i in s])
        self.trace_int = np.hstack([s_i.TraceInterval for s_i in s])
        self.pressure = np.hstack([s_i.Pressure for s_i in s])


class SInfo:

    def __init__(self, lines):
        self.Version = struct.unpack('<H', lines[0:2])[0] / 100
        self.fn = b''.join(struct.unpack('<64c', lines[2:66])).rstrip(b'\x00')
        try:
            self.fn = self.fn.decode('utf-8')
        except UnicodeDecodeError:
            pass
        self.serialtime = struct.unpack('<d', lines[66:74])[0] + datetime.date.toordinal(datetime.date(1970, 1, 1)) + 366.
        self.timezone = struct.unpack('<H', lines[74:76])[0] / 1440
        self.nChannels = struct.unpack('<B', lines[76:77])[0]
        self.RecordMode = struct.unpack('<B', lines[77:78])[0]
        if self.RecordMode == 0:
            self.RecordModeString = 'Odometer'
        elif self.RecordMode == 1:
            self.RecordModeString = 'Stacks'
        elif self.RecordMode == 2:
            self.RecordModeString = 'Time'
        else:
            print('Unknown Record Mode')
        
        self.RecordInterval = struct.unpack('<H', lines[78:80])[0]
        # Number of stacks per trace
        self.NumberOfStacks = struct.unpack('<H', lines[80:82])[0]
        # Sampling Frequency in MHz
        self.SampFreq = struct.unpack('<H', lines[82:84])[0] * 1.0e6
        # Pretrigger depth
        self.PreTriggerDepth = struct.unpack('<H', lines[84:86])[0]
        # Postrigger depth
        self.PostTriggerDepth = struct.unpack('<H', lines[86:88])[0]

        # Trigger source (1 = Chan A, 2 = Chan B, -1 = External)
        self.TriggerSource = struct.unpack('<B', lines[88:89])[0]
        if self.TriggerSource == 1:
            self.triggersourceString = 'Channel A'
        elif self.TriggerSource == 2:
            self.triggersourceString = 'Channel B'
        elif self.TriggerSource == -1:
            self.triggersourceString = 'External'
        else:
            print('Unknown in Trigger Source')

        # Trigger slope (0 = positive, 1 = negative)
        self.TriggerSlope = struct.unpack('<B', lines[89:90])[0]
        if self.TriggerSlope == 0:
            self.TriggerSlopeString = 'Negative'
        elif self.TriggerSlope == 1:
            self.TriggerSlopeStrong = 'Positive'
        else:
            print('Unknown Trigger Slope')

        # External Trigger range (Full range in in mV)
        self.ExtTriggerRange = struct.unpack('<H', lines[90:92])[0]
        # External trigger coupling (0 = DC, 1 = AC)
        self.ExtTriggerCoupling = struct.unpack('<B', lines[92:93])[0]
        if self.ExtTriggerCoupling == 0:
            self.ExtTriggerCouplingString = 'DC'
        elif self.ExtTriggerCoupling == 1:
            self.ExtTriggerCouplingString = 'AC'
        else:
            print('Unknown External Trigger Coupling')

        # track the index we are reading from since can now vary based upon radar version
        self.offset = 93

        # Odometer calibration constant(meters per trigger)
        # (Only available for pre 3.21 version)
        if self.Version < 3.21:
            self.OdometerCalibration = struct.unpack('<H', lines[self.offset:self.offset + 2])[0]
            self.offset += 2

        # Nominal Frequency (MHz)
        if self.Version < 3.8:
            self.NominalFrequency = struct.unpack('<h', lines[self.offset:self.offset + 2])[0]
            self.offset += 2
        else:
            self.NominalFrequency = struct.unpack('<f', lines[self.offset:self.offset + 4])[0]
            self.offset += 4
        # Antenna separation (m)
        self.AntennaSeparation = struct.unpack('<f', lines[self.offset:self.offset + 4])[0]
        self.offset += 4

        # Read and toss extra blank space
        if self.Version < 3.6:
            self.offset += 27

        # ==================== Channel Headers ==================== %
        for nn in range(self.nChannels):
            # Channel number
            nChan = struct.unpack('<B', lines[self.offset:self.offset + 1])[0]
            self.offset += 1
                        
            if nChan != nn + 1:
                raise ValueError('Corrupt Channel header, ({:d} != {:d})'.format(nChan, nn))

            # Construct channel name
            setattr(self, 'Channel{:d}'.format(nChan), Channel(lines, self.offset, self.Version))
            self.offset = getattr(self, 'Channel{:d}'.format(nChan)).offset

        # Specifics of data for each channel
        if self.Version < 3.6:
            if self.Version < 3.2:
                trace_record_len = 21550
            else:
                trace_record_len = 21552
        else:
            trace_record_len = 21548
        self.snum = self.PreTriggerDepth + self.PostTriggerDepth
        self.tnum = (len(lines) - self.offset) // self.nChannels // trace_record_len


class Channel:

    def __init__(self, lines, offset, version):
        # Full voltage range in mV
        self.VoltRange = struct.unpack('<H', lines[offset: offset + 2])[0]
        offset += 2

        # Channel Impedance (0 = 50 Ohm, 1 = 1 MOhm)
        self.Impedance = struct.unpack('<B', lines[offset:offset + 1])[0]
        offset += 1
        if self.Impedance == 0:
            self.ImpedanceString = '1 MOhm'
        elif self.Impedance == 1:
            self.ImpedanceString = '50 Ohm'
        else:
            print('Unknown Impedance for Channel {:d}'.format(nChan))

        # Channel coupling (0 = DC, 1 = AC)
        self.Coupling = struct.unpack('<B', lines[offset:offset + 1])[0]
        offset += 1
        if self.Coupling == 0:
            self.CouplingString = 'DC'
        elif self.Coupling == 1:
            self.CouplingString = 'AC'
        else:
            print('Unknown Coupling for Channel {:d}'.format(nChan))
        
        # Read and toss extra blank space
        # (Only needed for pre 3.6 version)
        if version < 3.6:
            offset += 27
        self.offset = offset


class ChannelData:

    def __init__(self, lines, sinfo):
        # Travel time
        self.TravelTime = np.arange(-sinfo.PreTriggerDepth, (sinfo.PostTriggerDepth)) * 1. / sinfo.SampFreq

        # I am going to try to preallocate because I think it should be possible
        # This is not actually true since we could have comments
        self.nTrace = np.zeros((sinfo.tnum, ))
        self.Time = np.zeros((sinfo.tnum, ))
        self.TraceInterval = np.zeros((sinfo.tnum, ))
        self.TriggerLevel = np.zeros((sinfo.tnum, ))
        self.lat = np.zeros((sinfo.tnum, ))
        self.long = np.zeros((sinfo.tnum, ))
        self.Altitude = np.zeros((sinfo.tnum, ))
        self.GPSResolution = np.zeros((sinfo.tnum, ))
        self.Data = np.zeros((sinfo.snum, sinfo.tnum))

        # These will often be empty, but leave here so we don't have missing attributes
        self.Odometer = np.zeros((sinfo.tnum, ))
        self.Pressure = np.zeros((sinfo.tnum, ))

    def read_trace(self, lines, sinfo, nTrc):
        nHeaderType = struct.unpack('<B', lines[sinfo.offset: sinfo.offset + 1])[0]
        nChannel = struct.unpack('<B', lines[sinfo.offset + 1: sinfo.offset + 2])[0]
        # print(nHeaderType, nChannel)

        # I'm rolling with a separate offset counter so i can work on guessing at tnum
        offset = 2

        # Trace number in file set
        self.nTrace[nTrc] = struct.unpack('<i', lines[sinfo.offset + offset:sinfo.offset + offset + 4])[0]
        offset += 4

        # Decimal day from 1 Jan 1970. 
        nTime = struct.unpack('<d', lines[sinfo.offset + offset:sinfo.offset + offset + 8])[0]
        offset += 8
        # We add an offset to 1 Jan 1970 to get MATLAB date numbers
        self.Time[nTrc] = nTime + datetime.date.toordinal(datetime.date(1970, 1, 1)) + 366.

        # Stacks/trace unless record mode is stacks, when it is
        # time/trace
        self.TraceInterval[nTrc] = struct.unpack('<f', lines[sinfo.offset + offset:sinfo.offset + offset + 4])[0]
        offset += 4

        # Trigger level in percentage of input range in mV
        self.TriggerLevel[nTrc] = struct.unpack('<H', lines[sinfo.offset + offset:sinfo.offset + offset + 2])[0]
        offset += 2

        if sinfo.Version < 3.21:
            # Odometer readings (0 if no Odometer is used)
            self.Odometer[nTrc] = struct.unpack('<f', lines[sinfo.offset + offset:sinfo.offset + offset + 4])[0]
            offset += 4

            # Pressure gauge (0 if no pressure gauge is used)
            self.Pressure[nTrc] = struct.unpack('<f', lines[sinfo.offset + offset:sinfo.offset + offset + 4])[0]
            offset += 4

        # GPS Latitude
        self.lat[nTrc] = struct.unpack('<d', lines[sinfo.offset + offset:sinfo.offset + offset + 8])[0]
        offset += 8

        # GPS longitude
        self.long[nTrc] = struct.unpack('<d', lines[sinfo.offset + offset:sinfo.offset + offset + 8])[0]
        offset += 8

        # GPS altitude
        self.Altitude[nTrc] = struct.unpack('<f', lines[sinfo.offset + offset:sinfo.offset + offset + 4])[0]
        offset += 4

        # GPS accuracy
        self.GPSResolution[nTrc] = struct.unpack('<f', lines[sinfo.offset + offset:sinfo.offset + offset + 4])[0]
        offset += 4
        
        # Read and toss last blank bytes
        # (Only needed for pre 3.6 version)
        if sinfo.Version < 3.6:
            if sinfo.Version < 3.2:
                offset += 12
            else:
                offset += 14

        # If it is actual radar data, not a comment or marker
        if nHeaderType == 0:
            # Read the trace data
            newdata = struct.unpack('<{:d}h'.format(sinfo.snum), lines[sinfo.offset + offset:sinfo.offset + offset + 2 * sinfo.snum])
            offset += 2 * sinfo.snum
            # Store data
            # Total number of data points 
            self.Data[:, nTrc] = newdata
        elif nHeaderType == 1:
            # Toss marker information
            offset += 38

        sinfo.offset += offset


def load_olaf(fn, channel=1):
    return Olaf(fn, channel)
