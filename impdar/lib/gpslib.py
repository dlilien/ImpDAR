#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""
Some classes and functions to handle different types of GPS data

The workhorse of this library, nmea_info, is not designed to be created directly. See `gssilib.get_dzg_data`
"""
import numpy as np
import osr

from scipy.interpolate import interp1d


def get_utm_conversion(lat, lon):   
    def utm_getZone(longitude):
        return (int(1 + (longitude + 180.0) / 6.0))

    def utm_isNorthern(latitude):
        if (latitude < 0.0):
            return False
        else:
            return True

    utm_zone = utm_getZone(lon)
    is_northern = utm_isNorthern(lat)
   
    utm_cs = osr.SpatialReference()
    utm_cs.SetWellKnownGeogCS('WGS84')
    utm_cs.SetUTM(utm_zone, is_northern)

    wgs84_cs = utm_cs.CloneGeogCS()
    wgs84_cs.ExportToPrettyWkt()
   
    transform_WGS84_To_UTM = osr.CoordinateTransformation(wgs84_cs, utm_cs)
    return transform_WGS84_To_UTM.TransformPoints


class nmea_info:
    """Container for general information about lat, lon, etc.
    
    Attributes
    ----------
    lat: `np.ndarray`
        wgs84 latitude of points
    lon: `np.ndarray`
        wgs84 longitude of points
    x: `np.ndarray`
        Projected x coordinates of points
    y: `np.ndarray`
        Projected y coordinates of points
    z: `np.ndarray`
        Projected z coordinates of points
    """
    all_data = None
    lat = None
    lon = None
    qual = None
    sats = None
    x = None
    y = None
    z = None
    geo_offset = None
    times = None
    scans = None

    def rev(self):
        self.all_data = np.flipud(self.all_data)
        if self.lat is not None:
            self.lat = np.flip(self.lat)
        if self.lon is not None:
            self.lon = np.flip(self.lon)

    def get_all(self):
        self.glat()
        self.glon()
        self.gqual()
        self.gsats()
        self.gz()
        self.ggeo_offset()
        self.gtimes()
        self.get_utm()
        self.get_dist()

    def glat(self):
        self.lat = self.all_data[:, 2] * ((self.all_data[:, 1] - self.all_data[:, 1] % 100) / 100 + (self.all_data[:, 1] % 100) / 60)
        if self.y is None:
            self.y = self.lat
        return self.lat

    def glon(self):
        self.lon = self.all_data[:, 4] * ((self.all_data[:, 3] - self.all_data[:, 3] % 100) / 100 + (self.all_data[:, 3] % 100) / 60)
        if self.x is None:
            self.x = self.lon
        return self.lon

    def gqual(self):
        self.qual = self.all_data[:, 5]
        return self.qual

    def gsats(self):
        self.sats = self.all_data[:, 6]
        return self.sats

    def gz(self):
        self.z = self.all_data[:, 8]
        return self.z

    def ggeo_offset(self):
        self.geo_offset = self.all_data[:, 8]
        return self.geo_offset

    def gtimes(self):
        self.times = self.all_data[:, 0]
        return self.times

    def get_dist(self):
        if self.x is None:
            self.glat()
        if self.y is None:
            self.glon()
        self.dist = np.zeros((len(self.y), 1))
        for i in range(1, len(self.dist)):
            self.dist[i] = self.dist[i - 1] + np.sqrt((self.x[i] - self.x[i - 1]) ** 2.0 + (self.y[i] - self.y[i - 1]) ** 2.0)

    def get_utm(self):
        transform = get_utm_conversion(np.mean(self.lat), np.mean(self.lon))
        pts = np.array(transform(np.vstack((self.lon, self.lat)).transpose()))
        self.x, self.y = pts[:, 0], pts[:, 1]

    @property
    def dectime(self):
        s = self.times % 100
        m = (self.times % 10000 - s) / 100
        h = (self.times - m * 100 - s) / 10000
        return (h + m / 60.0 + s / 3600.0) / 24.0

def nmea_to_ll(list_of_sentences):
    """Take in a list of raw sentences in GGA format and return a lon, lat list"""
    def _gga_sentence_to_ll(sentence):
        vals = sentence.split(',')
        lat = float(vals[2])
        lon = float(vals[4])
        return (lon, lat)

    if list_of_sentences[0].split(',')[0] == '$GPGGA':
        return np.array([_gga_sentence_to_ll(sentence) for sentence in list_of_sentences]) / 100.0
    else:
        raise ValueError('I can only do gga sentences right now')


def nmea_all_info(list_of_sentences):
    """Return an object with the nmea info from a given list of sentences"""
    def _gga_sentence_split(sentence):
        all = sentence.split(',')
        numbers = list(map(float, all[1:3] + [1] + [all[4]] + [1] + all[6:10] + [all[11]]))
        if all[3] == 'S':
            numbers[2] = -1
        if all[5] == 'W':
            numbers[4] = -1
        return numbers

    if list_of_sentences[0].split(',')[0] == '$GPGGA':
        data = nmea_info()
        data.all_data = np.array([_gga_sentence_split(sentence) for sentence in list_of_sentences])
        return data
    else:
        print(list_of_sentences[0].split(',')[0])
        raise ValueError('I can only do gga sentences right now')


class RadarGPS(nmea_info):
    
    def __init__(self, gga, scans, trace_num):
        self.nmea_info = nmea_all_info(gga)
        self.nmea_info.scans = scans
        self.nmea_info.get_all()

        kgps_indx = np.hstack((np.array([0]), np.where(np.diff(self.nmea_info.times) != 0)[0]))
        self.lat = interp1d(self.nmea_info.scans[kgps_indx], self.nmea_info.lat[kgps_indx], kind='linear', fill_value='extrapolate')(trace_num)
        self.lon = interp1d(self.nmea_info.scans[kgps_indx], self.nmea_info.lon[kgps_indx], kind='linear', fill_value='extrapolate')(trace_num)
        self.z = interp1d(self.nmea_info.scans[kgps_indx], self.nmea_info.z[kgps_indx], kind='linear', fill_value='extrapolate')(trace_num)
        self.times = interp1d(self.nmea_info.scans[kgps_indx], self.nmea_info.times[kgps_indx], kind='linear', fill_value='extrapolate')(trace_num)
        self.get_utm()
        self.get_dist()
