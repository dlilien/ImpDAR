#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Read the data from a St. Olaf/Gecko file
"""
import io
import struct
import datetime
import numpy as np

from ..RadarData import RadarData
from .loading_utils import common_start


class SInfo:
    """Information about a single profile line"""

    def __init__(self, lines):
        """Get information about the collection

        Parameters
        ----------
        lines: bytes
            The binary data to read
        """
        self.version = struct.unpack('<H', lines[0:2])[0] / 100
        self.fn_in = b''.join(
            struct.unpack('<64c', lines[2:66])).rstrip(b'\x00')
        try:
            self.fn_in = self.fn_in.decode('utf-8')
        except UnicodeDecodeError:
            pass

        # get time in useful way
        self.serialtime = struct.unpack('<d', lines[66:74])[0]
        self.serialtime += datetime.date.toordinal(datetime.date(
            1970, 1, 1)) + 366.

        self.timezone = struct.unpack('<H', lines[74:76])[0] / 1440
        self.n_channels = struct.unpack('<B', lines[76:77])[0]
        self.record_mode = struct.unpack('<B', lines[77:78])[0]
        self.get_record_mode()

        self.record_interval = struct.unpack('<H', lines[78:80])[0]
        # Number of stacks per trace
        self.number_of_stacks = struct.unpack('<H', lines[80:82])[0]
        # Sampling Frequency in MHz
        self.samp_freq = struct.unpack('<H', lines[82:84])[0] * 1.0e6
        # Pretrigger depth
        self.pre_trigger_depth = struct.unpack('<H', lines[84:86])[0]
        # Postrigger depth
        self.post_trigger_depth = struct.unpack('<H', lines[86:88])[0]

        # Trigger source (1 = Chan A, 2 = Chan B, -1 = External)
        self.trigger_source = struct.unpack('<B', lines[88:89])[0]
        self.get_trigger_source_string()

        # Trigger slope (0 = positive, 1 = negative)
        self.trigger_slope = struct.unpack('<B', lines[89:90])[0]
        self.get_trigger_slope_string()

        # External Trigger range (Full range in in mV)
        self.ext_trigger_range = struct.unpack('<H', lines[90:92])[0]
        # External trigger coupling (0 = DC, 1 = AC)
        self.ext_trigger_coupling = struct.unpack('<B', lines[92:93])[0]
        if self.ext_trigger_coupling == 0:
            self.ext_trigger_coupling_string = 'DC'
        elif self.ext_trigger_coupling == 1:
            self.ext_trigger_coupling_string = 'AC'
        else:
            print('Unknown External Trigger Coupling')

        # track the index we are reading from since can now vary on version
        self.offset = 93

        # Odometer calibration constant(meters per trigger)
        # (Only available for pre 3.21 version)
        if self.version < 3.21:
            self.odometer_calibration = struct.unpack('<H', lines[
                self.offset:self.offset + 2])[0]
            self.offset += 2

        # Nominal Frequency (MHz)
        if self.version < 3.8:
            self.nominal_frequency = struct.unpack('<h', lines[
                self.offset:self.offset + 2])[0]
            self.offset += 2
        else:
            self.nominal_frequency = struct.unpack('<f', lines[
                self.offset:self.offset + 4])[0]
            self.offset += 4
        # Antenna separation (m)
        self.antenna_separation = struct.unpack('<f', lines[
            self.offset:self.offset + 4])[0]
        self.offset += 4

        # Read and toss extra blank space
        if self.version < 3.6:
            self.offset += 27

        # ==================== Channel Headers ==================== %
        for i in range(self.n_channels):
            # Channel number
            n_chan = struct.unpack('<B', lines[self.offset:self.offset + 1])[0]
            self.offset += 1

            if n_chan != i + 1:
                raise ValueError(
                    'Corrupt Channel header, ({:d} != {:d})'.format(n_chan, i))

            # Construct channel name
            setattr(self,
                    'Channel{:d}'.format(n_chan),
                    Channel(lines, self.offset, self.version, n_chan))
            self.offset = getattr(self, 'Channel{:d}'.format(n_chan)).offset

        # Specifics of data for each channel
        if self.version < 3.6:
            if self.version < 3.2:
                trace_record_len = 21550
            else:
                trace_record_len = 21552
        else:
            trace_record_len = 21548
        # define number of samples and traces to preallocate arrays
        self.snum = self.pre_trigger_depth + self.post_trigger_depth
        self.tnum = (len(lines) - self.offset
                     ) // self.n_channels // (2*self.snum)

    def get_trigger_source_string(self):
        """Turn an integer source info into a string"""
        if self.trigger_source == 1:
            self.trigger_source_string = 'Channel A'
        elif self.trigger_source == 2:
            self.trigger_source_string = 'Channel B'
        elif self.trigger_source == -1:
            self.trigger_source_string = 'External'
        else:
            print('Unknown in Trigger Source')

    def get_trigger_slope_string(self):
        """Turn the integer slope info into a string"""
        if self.trigger_slope == 0:
            self.trigger_slope_string = 'Negative'
        elif self.trigger_slope == 1:
            self.trigger_slope_string = 'Positive'
        else:
            print('Unknown Trigger Slope')

    def get_record_mode(self):
        """Integer record mode to string"""
        if self.record_mode == 0:
            self.record_mode_string = 'Odometer'
        elif self.record_mode == 1:
            self.record_mode_string = 'Stacks'
        elif self.record_mode == 2:
            self.record_mode_string = 'Time'
        else:
            print('Unknown Record Mode')


class Channel:
    """Information about a channel."""
    def __init__(self, lines, offset, version, n_chan):
        # Full voltage range in mV
        self.volt_range = struct.unpack('<H', lines[offset: offset + 2])[0]
        offset += 2

        # Channel Impedance (0 = 50 Ohm, 1 = 1 MOhm)
        self.impedance = struct.unpack('<B', lines[offset:offset + 1])[0]
        offset += 1
        if self.impedance == 0:
            self.impedance_string = '1 MOhm'
        elif self.impedance == 1:
            self.impedance_string = '50 Ohm'
        else:
            print('Unknown Impedance for Channel {:d}'.format(n_chan))

        # Channel coupling (0 = DC, 1 = AC)
        self.coupling = struct.unpack('<B', lines[offset:offset + 1])[0]
        offset += 1
        if self.coupling == 0:
            self.coupling_string = 'DC'
        elif self.coupling == 1:
            self.coupling_string = 'AC'
        else:
            print('Unknown Coupling for Channel {:d}'.format(n_chan))

        # Read and toss extra blank space
        # (Only needed for pre 3.6 version)
        if version < 3.6:
            offset += 27
        self.offset = offset


class ChannelData:
    """Full data for radar channel."""

    def __init__(self, lines, sinfo):
        """Read in binary data.

        Parameters
        ----------
        lines: bytes
            the binary data
        sinfo: SInfo
            information needed to make sense of the binary data
        """
        # Travel time
        self.travel_time = np.arange(-sinfo.pre_trigger_depth,
                                     (sinfo.post_trigger_depth)
                                     ) * 1. / sinfo.samp_freq

        # I am going to preallocate because I think it should be possible
        # This is not actually true since we could have comments
        self.n_trace = np.zeros((sinfo.tnum, ))
        self.time = np.zeros((sinfo.tnum, ))
        self.trace_interval = np.zeros((sinfo.tnum, ))
        self.trigger_level = np.zeros((sinfo.tnum, ))
        self.lat = np.zeros((sinfo.tnum, ))
        self.long = np.zeros((sinfo.tnum, ))
        self.altitude = np.zeros((sinfo.tnum, ))
        self.gps_resolution = np.zeros((sinfo.tnum, ))
        self.data = np.zeros((sinfo.snum, sinfo.tnum))

        # These will often be empty, but leave here so we don't have missing attributes
        self.odometer = np.zeros((sinfo.tnum, ))
        self.pressure = np.zeros((sinfo.tnum, ))

    def read_trace(self, lines, sinfo, n_trc):
        """Ingest binary trace information

        Parameters
        ----------
        lines: bytes
            The binary data
        sinfo: SInfo
            The overall collection info
        n_trc: int
            The trace under consideration
        """
        n_header_type = struct.unpack('<B', lines[sinfo.offset: sinfo.offset + 1])[0]

        # I'm rolling with a separate offset counter so i can work on guessing at tnum
        offset = 2

        # Trace number in file set
        self.n_trace[n_trc] = struct.unpack('<i', lines[sinfo.offset + offset:
                                                        sinfo.offset + offset + 4])[0]
        offset += 4

        # Decimal day from 1 Jan 1970.
        n_time = struct.unpack('<d', lines[sinfo.offset + offset:sinfo.offset + offset + 8])[0]
        offset += 8
        # We add an offset to 1 Jan 1970 to get MATLAB date numbers
        self.time[n_trc] = n_time + datetime.date.toordinal(datetime.date(1970, 1, 1)) + 366.

        # Stacks/trace unless record mode is stacks, when it is
        # time/trace
        self.trace_interval[n_trc] = struct.unpack('<f', lines[sinfo.offset + offset:
                                                               sinfo.offset + offset + 4])[0]
        offset += 4

        # Trigger level in percentage of input range in mV
        self.trigger_level[n_trc] = struct.unpack('<H', lines[sinfo.offset + offset:
                                                              sinfo.offset + offset + 2])[0]
        offset += 2

        if sinfo.version < 3.21:
            # Odometer readings (0 if no Odometer is used)
            self.odometer[n_trc] = struct.unpack('<f', lines[sinfo.offset + offset:
                                                             sinfo.offset + offset + 4])[0]
            offset += 4

            # Pressure gauge (0 if no pressure gauge is used)
            self.pressure[n_trc] = struct.unpack('<f', lines[sinfo.offset + offset:
                                                             sinfo.offset + offset + 4])[0]
            offset += 4

        # GPS Latitude
        self.lat[n_trc] = struct.unpack('<d', lines[sinfo.offset + offset:
                                                    sinfo.offset + offset + 8])[0]
        offset += 8

        # GPS longitude
        self.long[n_trc] = struct.unpack('<d', lines[sinfo.offset + offset:
                                                     sinfo.offset + offset + 8])[0]
        offset += 8

        # GPS altitude
        self.altitude[n_trc] = struct.unpack('<f', lines[sinfo.offset + offset:
                                                         sinfo.offset + offset + 4])[0]
        offset += 4

        # GPS accuracy
        self.gps_resolution[n_trc] = struct.unpack('<f', lines[sinfo.offset + offset:
                                                               sinfo.offset + offset + 4])[0]
        offset += 4

        # Read and toss last blank bytes
        # (Only needed for pre 3.6 version)
        if sinfo.version < 3.6:
            if sinfo.version < 3.2:
                offset += 12
            else:
                offset += 14

        # If it is actual radar data, not a comment or marker
        if n_header_type == 0:
            # Read the trace data
            newdata = struct.unpack('<{:d}h'.format(sinfo.snum),
                                    lines[sinfo.offset + offset:
                                          sinfo.offset + offset + 2 * sinfo.snum])
            offset += 2 * sinfo.snum
            # Store data
            # Total number of data points
            self.data[:, n_trc] = newdata
        elif n_header_type == 1:
            # Toss marker information
            offset += 38

        sinfo.offset += offset


def load_olaf(fns_olaf, channel=1):
    """Read data from a gecko recording"""
    olaf_data = RadarData(None)
    # We want to be able to use this step concatenate a series of files numbered by the controller
    if isinstance(fns_olaf, str):
        fns_olaf = [fns_olaf]
        olaf_data.fn = fns_olaf[0]
    else:
        olaf_data.fn = common_start(fns_olaf).rstrip('[')

    sinfo = []
    stacks = []
    for i, fn_i in enumerate(fns_olaf):
        # We are going to follow the general format that was used by storead_script_v36
        with io.open(fn_i, 'rb') as fid:
            lines = fid.read()

        # Header information
        sinfo.append(SInfo(lines))

        # Data is stored per trace. Make a container
        s_i = [ChannelData(lines, sinfo[i]) for j in range(sinfo[i].n_channels)]

        # Read trace-by-trace, channel-by-channel
        for n_trc in range(sinfo[i].tnum):
            try:
                for s_j in s_i:
                    s_j.read_trace(lines, sinfo[i], n_trc)
            except:
                continue

        stacks.append(s_i[channel - 1])

    # I don't know if we actually want to do this, but the filenaming scheme is wacky and this
    # will make any logical collection look good
    sort_idx = np.argsort(np.array([(lambda x: x.serialtime)(s) for s in sinfo]))
    sinfo = [sinfo[i] for i in sort_idx]
    stacks = [stacks[i] for i in sort_idx]

    # Data and things we derive from it
    olaf_data.chan = channel
    olaf_data.data = np.hstack([s_i.data for s_i in stacks])
    olaf_data.snum = olaf_data.data.shape[0]
    olaf_data.tnum = olaf_data.data.shape[1]
    olaf_data.trace_num = np.arange(olaf_data.tnum) + 1

    # Now merge the data into the normal format
    olaf_data.dt = 1. / sinfo[0].samp_freq
    olaf_data.fns_in = sinfo[0].fn_in
    olaf_data.ant_sep = sinfo[0].antenna_separation
    olaf_data.freq = sinfo[0].nominal_frequency
    olaf_data.travel_time = stacks[0].travel_time * 1.0e6
    olaf_data.trig_level = stacks[0].trigger_level
    olaf_data.trig = sinfo[0].pre_trigger_depth * np.ones(olaf_data.tnum)
    olaf_data.fnames = [si.fn_in for si in sinfo]

    # Other variables that need concatenating
    olaf_data.decday = np.hstack([s_i.time for s_i in stacks])
    olaf_data.elev = np.hstack([s_i.altitude for s_i in stacks])
    olaf_data.lat = np.hstack([s_i.lat for s_i in stacks])
    olaf_data.long = np.hstack([s_i.long for s_i in stacks])
    olaf_data.trace_int = np.hstack([s_i.trace_interval for s_i in stacks])
    olaf_data.pressure = np.hstack([s_i.pressure for s_i in stacks])
    try:
        olaf_data.get_projected_coords()
    except ImportError:
        pass
    olaf_data.check_attrs()
    return olaf_data
