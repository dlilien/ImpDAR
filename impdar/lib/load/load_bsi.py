#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""Load BSI IceRadar h5 files and convert to the .mat ImpDAR file."""

import os
import numpy as np
from scipy.interpolate import interp1d

from ..RadarData import RadarData
import re
from .. import gpslib
import datetime


try:
    import h5py
    H5 = True
except ImportError:
    H5 = False


def _xmlGetVal(xml, name):
    """Look up a value in an XML fragment. Mod from Nat Wilson's irlib."""
    m = re.search(r'<Name>{0}</Name>[\r]?\n<Val>'.format(
        name.replace(' ', r'\s')), xml, flags=re.IGNORECASE)
    if m is not None:
        tail = xml[m.span()[1]:]
        end = tail.find('</Val')
        return tail[:end]
    else:
        return None


def _dm2dec(dms):
    """Convert the degree - decimal minute GGA to a decimal."""
    return ((dms - dms % 100) / 100 + (dms % 100) / 60)


def _dt_from_comment(dset):
    """Return the day when the first location in the dset was collected."""
    low_level_group = h5py.h5g.open(dset['location_0'].id, b'.')
    group_comment = low_level_group.get_comment(b'.').decode('utf-8')
    # this string is variable length, so we need complicated parsing
    group_comment = group_comment[group_comment.find(']') + 1:]
    group_comment = group_comment[group_comment.find(']') + 1:]
    group_comment = group_comment[:group_comment.find(' ')]
    dmy = list(map(int, group_comment.split('/')))
    return datetime.datetime(dmy[2], dmy[0], dmy[1], 0, 0, 0)


def load_bsi(fn_h5, nans=None, *args, **kwargs):
    """Load a BSI IceRadar file, which is just an h5 file, into ImpDAR."""
    if not H5:
        raise ImportError('You need H5 to load bsi')

    # This is different from most input filetypes in that we can have multiple
    # lines stored in one file. We want to preserve these since concatenation
    # may be illogical
    h5_data_list = []

    # open the h5 file
    with h5py.File(fn_h5, 'r') as f_in:
        dset_names = [key for key in f_in.keys()]
        for dset_name in dset_names:
            print('Loading {:s} from {:s}'.format(dset_name, fn_h5))
            # Just in case there is something else that can be here
            if 'line_' not in dset_name:
                continue

            dset = f_in[dset_name]
            h5_data = RadarData(None)
            # We need this for logical file naming later on
            h5_data.fn = os.path.splitext(fn_h5)[0] + dset_name + '.h5'
            h5_data.tnum = len(list(dset.keys()))
            h5_data.snum = len(
                dset['location_0']['datacapture_0']['echogram_0'])
            lat = np.zeros((h5_data.tnum,))
            lon = np.zeros((h5_data.tnum,))
            h5_data.elev = np.zeros((h5_data.tnum,))
            time = np.zeros((h5_data.tnum,))
            h5_data.data = np.zeros((h5_data.snum, h5_data.tnum))

            if type(dset['location_0']['datacapture_0']['echogram_0'].attrs['Digitizer-MetaData_xml']) == str:
                digitizer_data = dset['location_0']['datacapture_0'][
                    'echogram_0'].attrs['Digitizer-MetaData_xml']
            else:
                digitizer_data = dset['location_0']['datacapture_0'][
                    'echogram_0'].attrs['Digitizer-MetaData_xml'].decode('utf-8')
            for location_num in range(h5_data.tnum):
                # apparently settings can change mid-line
                nsamps = dset['location_{:d}'.format(location_num)]['datacapture_0']['echogram_0'].shape[0]
                if nsamps > h5_data.snum:
                    h5_data.data = np.vstack((h5_data.data, np.zeros((nsamps - h5_data.snum, h5_data.tnum))))
                    h5_data.snum = nsamps
                h5_data.data[:, location_num] = dset[
                    'location_{:d}'.format(location_num)][
                        'datacapture_0']['echogram_0']
                if type(dset['location_{:d}'.format(location_num)][
                    'datacapture_0']['echogram_0'].attrs[
                        'GPS Cluster- MetaData_xml']) == str:
                    gps_data = dset['location_{:d}'.format(location_num)][
                        'datacapture_0']['echogram_0'].attrs[
                            'GPS Cluster- MetaData_xml']
                else:
                    gps_data = dset['location_{:d}'.format(location_num)][
                        'datacapture_0']['echogram_0'].attrs[
                            'GPS Cluster- MetaData_xml'].decode('utf-8')
                if (float(_xmlGetVal(gps_data, 'GPS Fix valid')) > 0) and (
                        float(_xmlGetVal(gps_data, 'GPS Message ok')) > 0):
                    # sometimes, there are bad entries that are unmarked
                    try:
                        lat[location_num] = float(_xmlGetVal(gps_data, 'Lat_N'))
                        lon[location_num] = float(
                            _xmlGetVal(gps_data, 'Long_ W'))
                        time[location_num] = float(
                            _xmlGetVal(gps_data, 'GPS_timestamp_UTC'))
                        h5_data.elev[location_num] = float(
                            _xmlGetVal(gps_data, 'Alt_asl_m'))
                    except:
                        lat[location_num] = np.nan
                        lon[location_num] = np.nan
                        time[location_num] = np.nan
                        h5_data.elev[location_num] = np.nan
                else:
                    lat[location_num] = np.nan
                    lon[location_num] = np.nan
                    time[location_num] = np.nan
                    h5_data.elev[location_num] = np.nan

            h5_data.dt = 1.0 / float(
                _xmlGetVal(digitizer_data, ' Sample Rate'))
            h5_data.travel_time = np.arange(h5_data.snum) * h5_data.dt * 1.0e6

            # Other information that ImpDAR currently cannot use
            # _xmlGetVal(digitizer_data, 'vertical range')
            h5_data.trig_level = float(
                _xmlGetVal(digitizer_data, 'trigger level'))
            time_offset = float(_xmlGetVal(digitizer_data, 'relativeInitialX'))
            h5_data.travel_time = h5_data.travel_time + time_offset * 1.0e6


            mask = ~np.isnan(time)
            if nans == 'interp':
                if np.any(~mask):
                    print('Interpolating traces with bad GPS in {:s}'.format(dset_name))
                # get rid of bad GPS locations
                h5_data.trace_num = np.arange(h5_data.tnum).astype(int) + 1
                time = interp1d(h5_data.trace_num[mask],
                                time[mask],
                                fill_value='extrapolate')(h5_data.trace_num)
                h5_data.lat = interp1d(h5_data.trace_num[mask],
                                       _dm2dec(lat[mask]),
                                       fill_value='extrapolate')(h5_data.trace_num)
                h5_data.long = interp1d(h5_data.trace_num[mask],
                                        -_dm2dec(lon[mask]),
                                        fill_value='extrapolate')(h5_data.trace_num)
                h5_data.elev = interp1d(h5_data.trace_num[mask],
                                        h5_data.elev[mask],
                                        fill_value='extrapolate')(h5_data.trace_num)
            elif nans == 'delete':
                if np.any(~mask):
                    print('Deleting traces with bad GPS in {:s}'.format(dset_name))
                h5_data.lat = _dm2dec(lat[mask])
                h5_data.long = -_dm2dec(lon[mask])
                h5_data.elev = h5_data.elev[mask]
                h5_data.data = h5_data.data[:, mask]
                time = time[mask]

                # Deal with this here in case tnum changed due to bad traces
                h5_data.tnum = h5_data.data.shape[1]
                h5_data.trace_num = np.arange(h5_data.tnum).astype(int) + 1
            else:
                h5_data.lat = _dm2dec(lat)
                h5_data.long = -_dm2dec(lon)
                h5_data.trace_num = np.arange(h5_data.tnum).astype(int) + 1

            h5_data.trig = np.floor(np.ones((h5_data.tnum, )) * np.abs(
                time_offset) / h5_data.dt)

            # need to access a comment to get the day
            day_collection = _dt_from_comment(dset)
            day_offset = (day_collection - datetime.datetime(1, 1, 1, 0, 0, 0)
                          ).days
            h5_data.decday = gpslib.hhmmss2dec(time) + day_offset

            try:
                h5_data.get_projected_coords()
            except ImportError:
                temp_x = h5_data.long * 110. * np.cos(h5_data.lat)
                temp_y = 110. * h5_data.lat
                h5_data.dist = np.hstack((
                    [0], np.cumsum(np.sqrt(
                        np.diff(temp_x) ** 2.0 + np.diff(temp_y) ** 2.0))))
                print('Geographic coordinate conversion skipped: no ogr')

            h5_data.trace_int = np.hstack(
                (np.array(np.nanmean(np.diff(h5_data.dist))),
                 np.diff(h5_data.dist)))
            h5_data.pressure = np.zeros_like(h5_data.lat)

            h5_data.chan = 0
            h5_data.check_attrs()
            h5_data_list.append(h5_data)
    return h5_data_list
