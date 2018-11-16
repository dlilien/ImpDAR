#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2016 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3 license.

"""
Some classes and functions to handle different types of GPS data

The workhorse of this library, nmea_info, is not designed to be created directly. Use RadarGPS class, which has an __init__ method, instead.

Additional methods in this library are used to read the filetypes from StoDeep. These can then be used to redo the GPS info on another object
"""
import numpy as np
try:
    import osr
    conversions_enabled = True
except ImportError:
    conversions_enabled = False

from scipy.interpolate import interp1d

if conversions_enabled:
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
else:
    def get_utm_conversion(lat, lon):
        raise ImportError('Cannot convert coordinates: osr not importable')


class nmea_info:
    """Container for general information about lat, lon, etc.
    
    Attributes
    ----------
    lat: `ndarray <https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html>`_
        wgs84 latitude of points
    lon: `ndarray <https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html>`_
        wgs84 longitude of points
    x: `ndarray <https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html>`_
        Projected x coordinates of points
    y: `ndarray <https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html>`_
        Projected y coordinates of points
    z: `ndarray <https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html>`_
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
        if conversions_enabled:
            self.get_utm()
        self.get_dist()

    def glat(self):
        if self.lat is None:
            self.lat = self.all_data[:, 2] * ((self.all_data[:, 1] - self.all_data[:, 1] % 100) / 100 + (self.all_data[:, 1] % 100) / 60)
        if self.y is None:
            self.y = self.lat
        return self.lat

    def glon(self):
        if self.lon is None:
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
            self.glon()
        if self.y is None:
            self.glat()
        self.dist = np.zeros((len(self.y), ))
        self.dist[1:] = np.cumsum(np.sqrt((self.x[1:] - self.x[:-1]) ** 2.0 + (self.y[1:] - self.y[:-1]) ** 2.0)) / 1000.0

    def get_utm(self):
        transform = get_utm_conversion(np.nanmean(self.lat), np.nanmean(self.lon))
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
        numbers = list(map(lambda x: float(x) if x != '' else 0, all[1:3] + [1] + [all[4]] + [1] + all[6:10] + [all[11]]))
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

        kgps_indx = np.hstack((np.array([0]), 1 + np.where(np.logical_and(np.diff(self.nmea_info.times) != 0, np.logical_and(~np.isnan(self.nmea_info.times[1:]), np.diff(self.nmea_info.scans) != 0)))[0]))
        self.lat = interp1d(self.nmea_info.scans[kgps_indx], self.nmea_info.lat[kgps_indx], kind='linear', fill_value='extrapolate')(trace_num)
        self.lon = interp1d(self.nmea_info.scans[kgps_indx], self.nmea_info.lon[kgps_indx], kind='linear', fill_value='extrapolate')(trace_num)
        self.z = interp1d(self.nmea_info.scans[kgps_indx], self.nmea_info.z[kgps_indx], kind='linear', fill_value='extrapolate')(trace_num)
        self.times = interp1d(self.nmea_info.scans[kgps_indx], self.nmea_info.times[kgps_indx], kind='linear', fill_value='extrapolate')(trace_num)
        if conversions_enabled:
            self.get_utm()
        self.get_dist()


def kinematic_gps_control(dats, lat, lon, elev, decday, offset=0.0):
    """Use new, better GPS data for lat, lon, and elevation

    The interpolation in this function is done using the time since the radar has accurate timing from its GPS. The old version of this function in StoDeep required redundant variables (x_coord, y_coord, dist). I've dropped that dependency.

    Parameters
    ----------
    dats: list of impdar.RadarData or impdar.RadarData
        The data to act upon
    lat: :class:`numpy.ndarray`
        Latitude from kinematic
    lon: :class:`numpy.ndarray`
        Longitude from kinematic
    elev: :class:`numpy.ndarray`
        Elevation from kinematic
    decday: :class:`numpy.ndarray`
        Decimal day. You need to reference this to match up with what the radar uses using offset
    offset: float, optional
        Translate the GPS times by this amount for alignment with the radar
    """
    int_lat = interp1d(decday, lat)
    int_long = interp1d(decday, lon)
    int_elev = interp1d(decday, elev)
    if type(dats) not in [list, tuple]:
        dats = [dats]
    for dat in dats:
        dat.lat = int_lat(dat.decday)
        dat.long = int_long(dat.decday)
        dat.elev = int_elev(dat.decday)
        gpsdat = nmea_info()
        gpsdat.lat = dat.lat
        gpsdat.lon = dat.long
        gpsdat.get_utm()
        gpsdat.get_dist()
        dat.x_coord, dat.y_coord = gpsdat.x, gpsdat.y
        dat.dist = gpsdat.dist


def kinematic_gps_mat(dats, mat_fn, offset=0.0):
    """Use a matlab file with gps info to redo radar GPS

    Parameters
    ----------
    dats: impdar.RadarData or list of impdar.RadarData
        The radar data to redo
    mat_fn: str
        The matlab file, containing lat, long, elev, and decday, to use
    offset: float, optional
        Change decday by this much to match the radar's gps
    """
    from scipy.io import loadmat
    mat = loadmat(mat_fn)
    for val in ['lat', 'long', 'elev', 'decday']:
        if val not in mat:
            raise ValueError('{:s} needs to be contained in matlab input file'.format(val))
    kinematic_gps_control(dats, mat['lat'].flatten(), mat['long'].flatten(), mat['elev'].flatten(), mat['decday'].flatten(), offset=offset)


def kinematic_gps_csv(dats, csv_fn, offset=0, names='decday,long,lat,elev', **genfromtxt_flags):
    """Use a csv gps file to redo the GPS on radar data.

    The csv is read using numpy.genfromtxt, which supports a number of options. One, 'names', is set explicitly by the argument 'names' to this function: you can change the value of that string to True to read column names from the first post-header line in the file. You can also manually change the column names by giving a different comma-separated string. The names must contain 'decday', 'long', 'lat', and 'elev'.

    Parameters
    ----------
    dats: impdar.RadarData or list of impdar.RadarData
        The radar data to redo
    csv_fn: str
        The filename to act upon
    offset: float, optional
        Change decday by this much to match the radar's gps
    names: str, bool, list, or None
        names argument to numpy.genfromtxt used to read the csv.

    Any additional kwargs are passed to numpy.genfromtxt
    """
    data = np.genfromtxt(csv_fn, names=names, **genfromtxt_flags)
    for val in ['lat', 'long', 'elev', 'decday']:
        if val not in data.names:
            raise ValueError('{:s} needs to be contained in csv'.format(val))
    kinematic_gps_control(dats, data['lat'].flatten(), data['lon'].flatten(), data['decday'].flatten(), data['decday'].flatten(), offset=offset)


def interp(dats, spacing, fn=None, fn_type=None, offset=0.0, min_movement=1.0e-2, genfromtxt_kwargs={}, **kwargs):
    """Do kinematic GPS control then interpolate the data to constant spacing

    Parameters
    ----------
    spacing: float
        Target distance spacing in meters
    fn: str, optional
        If this is None, no control. Otherwise, this should be a mat or csv file with lat, long, elev, and decday fields
    fn_type: str, optional
        csv or mat? Ignored if fn is None. Will guess based on extension if this is None
    offset: float, optional
        move the decday by this much to match gps
    min_movement: float, optional
        use this separation to try to cull stationary entries in the gps
    genfromtxt_kwargs: dict, optional
        kwargs to pass to genfromtxt when reading a csv. Ignored otherwise.
    """
    if fn is not None:
        if fn_type == 'mat' or (fn_type is None and fn[-4:] == '.mat'):
            kinematic_gps_mat(dats, fn, offset=offset)
        if fn_type == 'csv' or (fn_type is None and fn[-4:] in ['.csv', '.txt']):
            kinematic_gps_csv(dats, fn, offset=offset, **genfromtxt_kwargs)
        else:
            raise ValueError('fn_type must be mat or csv')
    for dat in dats:
        dat.constant_space(spacing, min_movement=min_movement)
