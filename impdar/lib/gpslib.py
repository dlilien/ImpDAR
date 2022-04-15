#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3 license.

"""
Some classes and functions to handle different types of GPS data.

The workhorse of this library, nmea_info, is not designed to be created
directly. Use RadarGPS class, which has an __init__ method, instead.

Additional methods in this library are used to read the filetypes from StoDeep.
These can then be used to redo the GPS info on another object
"""
import numpy as np
try:
    import osr
    conversions_enabled = True
except ImportError:
    try:
        from osgeo import osr
        conversions_enabled = True
    except ImportError:
        conversions_enabled = False

from scipy.interpolate import interp1d

if conversions_enabled:
    def get_utm_conversion(lat, lon):
        """Retrun the gdal transform to convert coords."""
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

        # On newer versions of osr we need this, but on old versions it will fail
        try:
            utm_cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass

        wgs84_cs = utm_cs.CloneGeogCS()
        wgs84_cs.ExportToPrettyWkt()
        try:
            wgs84_cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass

        transform_WGS84_To_UTM = osr.CoordinateTransformation(wgs84_cs, utm_cs)
        return transform_WGS84_To_UTM.TransformPoints, utm_cs.ExportToPrettyWkt()

    def get_conversion(t_srs):
        out_cs = osr.SpatialReference()
        out_cs.SetFromUserInput(t_srs)
        try:
            out_cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass

        wgs84_cs = out_cs.CloneGeogCS()
        wgs84_cs.ExportToPrettyWkt()
        try:
            wgs84_cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass

        transform_WGS84_To_srs = osr.CoordinateTransformation(wgs84_cs, out_cs)
        return transform_WGS84_To_srs.TransformPoints, out_cs.ExportToPrettyWkt()

    def get_rev_conversion(t_srs):
        out_cs = osr.SpatialReference()
        out_cs.SetFromUserInput(t_srs)
        try:
            out_cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass

        wgs84_cs = out_cs.CloneGeogCS()
        wgs84_cs.ExportToPrettyWkt()
        try:
            wgs84_cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass

        transform_srs_to_WGS84 = osr.CoordinateTransformation(out_cs, wgs84_cs)
        return transform_srs_to_WGS84.TransformPoints, out_cs.ExportToPrettyWkt()

else:
    def get_utm_conversion(lat, lon):
        """Just raise an exception since we cannot really convert."""
        raise ImportError('Cannot convert coordinates: osr not importable')

    def get_conversion(t_srs):
        """Just raise an exception since we cannot really convert."""
        raise ImportError('Cannot convert coordinates: osr not importable')

    def get_rev_conversion(t_srs):
        """Just raise an exception since we cannot really convert."""
        raise ImportError('Cannot convert coordinates: osr not importable')


def hhmmss2dec(times):
    """Deal with one of the many weird time formats. 6-char string to day."""
    s = times % 100
    m = (times % 10000 - s) / 100
    h = (times - m * 100 - s) / 10000
    return (h + m / 60.0 + s / 3600.0) / 24.0


class nmea_info:
    """Container for general information about lat, lon, etc.

    Attributes
    ----------
    lat: np.ndarray
        wgs84 latitude of points
    lon: np.ndarray_
        wgs84 longitude of points
    x: np.ndarray
        Projected x coordinates of points
    y: np.ndarray
        Projected y coordinates of points
    z: np.ndarray
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

    def get_all(self):
        """Populate all the values from the input data."""
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
        """Populate lat(itude)."""
        if self.lat is None:
            self.lat = self.all_data[:, 2] * (
                (self.all_data[:, 1] - self.all_data[:, 1] % 100) / 100 + (
                    self.all_data[:, 1] % 100) / 60)
        if self.y is None:
            self.y = self.lat * 110000.0  # Temporary guess using earths radius
        return self.lat

    def glon(self):
        """Populate lon(gitude)."""
        if self.lon is None:
            self.lon = self.all_data[:, 4] * (
                (self.all_data[:, 3] - self.all_data[:, 3] % 100) / 100 + (
                    self.all_data[:, 3] % 100) / 60)
        if self.x is None:
            # Temporary guess using radius of the earth
            if self.lat is None:
                self.glat()
            self.x = self.lon * 110000.0 * \
                np.abs(np.cos(self.lat * np.pi / 180.0))
        return self.lon

    def gqual(self):
        """Populate qual(ity)."""
        self.qual = self.all_data[:, 5]
        return self.qual

    def gsats(self):
        """Populate sats (number of satellites)."""
        self.sats = self.all_data[:, 6]
        return self.sats

    def gz(self):
        """Populate z (elevation)."""
        self.z = self.all_data[:, 8]
        return self.z

    def ggeo_offset(self):
        """Populate geo_offset (Distance between ellipsoid and geoid)."""
        self.geo_offset = self.all_data[:, 8]
        return self.geo_offset

    def gtimes(self):
        """Populate times."""
        self.times = self.all_data[:, 0]
        return self.times

    def get_dist(self):
        """Calculate distance."""
        if self.y is None:
            self.glat()
        if self.x is None:
            self.glon()
        if conversions_enabled:
            self.get_utm()

        self.dist = np.zeros((len(self.y), ))
        self.dist[1:] = np.cumsum(np.sqrt(np.diff(self.x) ** 2.0 + np.diff(self.y) ** 2.0)) / 1000.0

    def get_utm(self):
        """Transform lat and lon to utm coords in a nice way."""
        transform, _ = get_utm_conversion(np.nanmean(self.lat),
                                          np.nanmean(self.lon))
        pts = np.array(transform(np.vstack((self.lon, self.lat)).transpose()))
        self.x, self.y = pts[:, 0], pts[:, 1]

    @property
    def dectime(self):
        """Convert the nasty 6-char time to something usable."""
        return hhmmss2dec(self.times)


def nmea_all_info(list_of_sentences):
    """
    Return an object with the nmea info from a given list of sentences.

    Parameters
    ----------
    list_of_sentences : list of strs
        NMEA output.

    Raises
    ------
    ValueError
        If the NMEA output does not contain GGA strings.

    Returns
    -------
    np.ndarray
        An array of the useful information in the NMEA sentences.
    """
    def _gga_sentence_split(sentence):
        all = sentence.split(',')
        if len(all) > 5:
            # We can have corrupted lines--just ignore these and continue
            try:
                numbers = list(map(lambda x: float(x) if x != '' else np.nan, all[1:3] + [1] + [all[4]] + [1] + all[6:10] + [all[11]]))
                if all[3] == 'S':
                    numbers[2] = -1
                if all[5] == 'W':
                    numbers[4] = -1
            except (ValueError, IndexError):
                numbers = [np.nan] * 10
        elif len(all) > 2:
            try:
                numbers = list(map(lambda x: float(x) if x != '' else np.nan, all[1:3] + [1]))
                if all[3] == 'S':
                    numbers[2] = -1
            except (ValueError, IndexError):
                numbers = [np.nan] * 10
        else:
            numbers = [np.nan] * 10
        return numbers

    if np.all([sentence.split(',')[0] == '$GPGGA' for sentence in list_of_sentences]):
        data = nmea_info()
        data.all_data = np.array([_gga_sentence_split(sentence)
                                  for sentence in list_of_sentences])
        return data
    else:
        raise ValueError('I can only do gga sentences right now')


class RadarGPS(nmea_info):
    """
    A container to make nmea info useful.

    This should handle frequency mismatch between radar and gps.

    Parameters
    ----------
    gga : list of strs
        The GPS data
    scans : TYPE
        DESCRIPTION.
    trace_num : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    def __init__(self, gga, scans, trace_num):

        self.nmea_info = nmea_all_info(gga)
        self.nmea_info.scans = scans
        self.nmea_info.get_all()

        kgps_mask = np.logical_and(~np.isnan(self.nmea_info.times[1:]),
                                   np.diff(self.nmea_info.scans) != 0)
        kgps_mask = np.logical_and(np.diff(self.nmea_info.times) != 0,
                                   kgps_mask)
        kgps_where = np.where(kgps_mask)[0]
        kgps_indx = np.hstack((np.array([0]), 1 + kgps_where))
        self.lat = interp1d(self.nmea_info.scans[kgps_indx],
                            self.nmea_info.lat[kgps_indx],
                            kind='linear',
                            fill_value='extrapolate')(trace_num)
        self.lon = interp1d(self.nmea_info.scans[kgps_indx],
                            self.nmea_info.lon[kgps_indx],
                            kind='linear',
                            fill_value='extrapolate')(trace_num)
        self.z = interp1d(self.nmea_info.scans[kgps_indx],
                          self.nmea_info.z[kgps_indx], kind='linear',
                          fill_value='extrapolate')(trace_num)
        self.times = interp1d(self.nmea_info.scans[kgps_indx],
                              self.nmea_info.times[kgps_indx],
                              kind='linear',
                              fill_value='extrapolate')(trace_num)
        if conversions_enabled:
            self.get_utm()
        self.get_dist()


def kinematic_gps_control(dats, lat, lon, elev, decday, offset=0.0,
                          extrapolate=False, guess_offset=True, old_gps_gaps=False):
    """Use new, better GPS data for lat, lon, and elevation.

    The interpolation in this function is done using the time since the radar
    has accurate timing from its GPS.
    The old version of this function in StoDeep required redundant variables
    (x_coord, y_coord, dist).
    I've dropped that dependency.

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
        Decimal day. You need to reference this to match up with what the
        radar uses using offset
    offset: float, optional
        Translate the GPS times by this amount for alignment with the radar
    extrapoloate: bool, optional
        If true, extrapolate data to fill values rather than using NaNs.
        Desirable for small offsets with the GPS, but dangerous since you
        can totally screw up the geolocation and not get an error.
        USE WITH CAUTION.
    guess_offset: bool, optional
        If true, ImpDAR will attempt to find the offset between the GPS and
        Radar times using the cross-correlation between
        the x coordinates in the two datasets. If the guess at the offset is
        nonzero, we look at 1000 offsets within 10% of
        the offset. Else we look at +/- 0.001 days

    """
    if extrapolate:
        fill_value = 'extrapolate'
    else:
        fill_value = np.NaN

    if type(dats) not in [list, tuple]:
        dats = [dats]

    for in_dat in [lat, lon, elev]:
        if len(decday) != len(in_dat):
            raise IndexError('lat, lon, elev, and decday must be the same len')
    offsets = [offset for i in dats]
    if guess_offset:
        print('CC search')
        for j, dat in enumerate(dats):
            decday_interp = dat.decday.copy()
            if old_gps_gaps:
                for i,dday in enumerate(decday_interp):
                    if np.min(abs(dday-decday))>1./(24*3600.):
                        decday_interp[i] = np.nan
                dat.lat[dat.lat == 0.] = np.nan
                dat.long[dat.long == 0.] = np.nan
            for i in range(5):
                if (min(lon % 360) - max(dat.long % 360)) > 0. or (min(dat.long % 360) - max(lon % 360)) > 0.:
                    raise ValueError('No overlap in longitudes')

                if offsets[j] != 0.0:
                    search_vals = np.linspace(-0.1 * abs(offsets[j]),
                                              0.1 * abs(offsets[j]),1001)
                else:
                    search_vals = np.linspace(-0.1, 0.1, 5001)

                cc_coeffs = np.zeros_like(search_vals)
                for i_search,inc_offset in enumerate(search_vals):
                    lat_interp = interp1d(decday + inc_offset + offsets[j], lat, kind='linear', fill_value=fill_value)(decday_interp)
                    long_interp = interp1d(decday + inc_offset + offsets[j], lon%360, kind='linear', fill_value=fill_value)(decday_interp)
                    idxs_lat = ~np.isnan(lat_interp) & ~np.isnan(dat.lat)
                    idxs_long = ~np.isnan(long_interp) & ~np.isnan(dat.long)
                    cc_coeffs[i_search] = np.corrcoef(lat_interp[idxs_lat], dat.lat[idxs_lat])[0, 1] + \
                                 np.corrcoef(long_interp[idxs_long], dat.long[idxs_long] % 360)[0, 1]
                offsets[j] += search_vals[np.argmax(cc_coeffs)]
                print('Maximum correlation at offset: {:f}'.format(offsets[j]))

    for j, dat in enumerate(dats):
        decday_interp = dat.decday.copy()
        lat_interpolator = interp1d(decday + offsets[j],
                           lat,
                           kind='linear',
                           fill_value=fill_value)
        long_interpolator = interp1d(decday + offsets[j],
                            lon % 360, kind='linear',
                            fill_value=fill_value)
        elev_interpolator = interp1d(decday + offsets[j],
                            elev,
                            kind='linear',
                            fill_value=fill_value)
        if old_gps_gaps:
            lat_interp = lat_interpolator(decday_interp)
            long_interp = long_interpolator(decday_interp)
            elev_interp = elev_interpolator(decday_interp)
            lat_interp[np.isnan(decday_interp)] = dat.lat[np.isnan(decday_interp)]
            long_interp[np.isnan(decday_interp)] = dat.long[np.isnan(decday_interp)]
            elev_interp[np.isnan(decday_interp)] = dat.elev[np.isnan(decday_interp)]
            dat.lat = lat_interp
            dat.long = long_interp%360
            dat.elev = elev_interp
        else:
            dat.lat = lat_interpolator(decday_interp)
            dat.long = long_interpolator(decday_interp)
            dat.elev = elev_interpolator(decday_interp)
        if conversions_enabled:
            dat.get_projected_coords()


def kinematic_gps_mat(dats, mat_fn, offset=0.0, extrapolate=False,
                      guess_offset=False, old_gps_gaps=False):
    """Use a matlab file with gps info to redo radar GPS.

    Parameters
    ----------
    dats: impdar.RadarData or list of impdar.RadarData
        The radar data to redo
    mat_fn: str
        The matlab file, containing lat, long, elev, and decday, to use
    offset: float, optional
        Change decday by this much to match the radar's gps
    extrapoloate: bool, optional
        If true, extrapolate data to fill values rather than using NaNs.
        Desirable for small offsets with the GPS,
        but dangerous since you can totally screw up the geolocation and
        not get an error.
        USE WITH CAUTION.
    """
    from scipy.io import loadmat
    mat = loadmat(mat_fn)
    for val in ['lat', 'long', 'elev', 'decday']:
        if val not in mat:
            raise ValueError('{:s} needs to be contained in matlab \
                               input file'.format(val))
    kinematic_gps_control(dats,
                          mat['lat'].flatten(),
                          mat['long'].flatten(),
                          mat['elev'].flatten(),
                          mat['decday'].flatten(),
                          offset=offset, extrapolate=extrapolate,
                          guess_offset=guess_offset, old_gps_gaps=old_gps_gaps)


def kinematic_gps_csv(dats, csv_fn, offset=0, names='decday,long,lat,elev',
                      extrapolate=False, guess_offset=False, old_gps_gaps=False,
                      **genfromtxt_flags):
    """Use a csv gps file to redo the GPS on radar data.

    The csv is read using numpy.genfromtxt, which supports a number of options.
    One, 'names', is set explicitly by the argument 'names' to this function:
        you can change the value of that string to True to read column names
        from the first post-header line in the file.
        You can also manually change the column names by giving a different
        comma-separated string.
        The names must contain 'decday', 'long', 'lat', and 'elev'.

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
    extrapoloate: bool, optional
        If true, extrapolate data to fill values rather than using NaNs.
        Desirable for small offsets with the GPS, but dangerous since you can
        totally screw up
        the geolocation and not get an error.
        USE WITH CAUTION.


    Any additional kwargs are passed to numpy.genfromtxt
    """
    data = np.genfromtxt(csv_fn, names=names, **genfromtxt_flags)
    kinematic_gps_control(dats,
                          data['lat'].flatten(),
                          data['long'].flatten(),
                          data['elev'].flatten(),
                          data['decday'].flatten(),
                          offset=offset,
                          extrapolate=extrapolate,
                          guess_offset=guess_offset,
                          old_gps_gaps=old_gps_gaps)


def interp(dats, spacing=None, fn=None, fn_type=None, offset=0.0,
           min_movement=1.0e-2, genfromtxt_kwargs={}, extrapolate=False,
           guess_offset=False, **kwargs):
    """Do kinematic GPS control then interpolate the data to constant spacing.

    Parameters
    ----------
    spacing: float, optional
        Target distance spacing in meters
    fn: str, optional
        If this is None, no control.
        Otherwise, this should be a mat or csv file with lat, long, elev,
        and decday fields
    fn_type: str, optional
        csv or mat? Ignored if fn is None.
        Will guess based on extension if this is None
    offset: float, optional
        move the decday by this much to match gps
    min_movement: float, optional
        use this separation to try to cull stationary entries in the gps
    genfromtxt_kwargs: dict, optional
        kwargs to pass to genfromtxt when reading a csv. Ignored otherwise.
    extrapoloate: bool, optional
        If true, extrapolate data to fill values rather than using NaNs.
        Desirable for small offsets with the GPS, but dangerous since you can
        totally screw up the geolocation and not get an error.
        USE WITH CAUTION.
    """
    if fn is not None:
        if fn_type == 'mat' or ((fn_type is None) and (fn[-4:] == '.mat')):
            kinematic_gps_mat(dats,
                              fn,
                              offset=offset,
                              extrapolate=extrapolate,
                              guess_offset=guess_offset)
        elif fn_type == 'csv' or (fn_type is None and fn[-4:] in ['.csv',
                                                                  '.txt']):
            kinematic_gps_csv(dats,
                              fn,
                              offset=offset,
                              extrapolate=extrapolate,
                              guess_offset=guess_offset,
                              **genfromtxt_kwargs)
        else:
            raise ValueError('Cannot identify fn filetype, must be mat or csv')
    if spacing is not None:
        for dat in dats:
            if dat.dist is None:
                kinematic_gps_control(dat,dat.lat,dat.long,dat.elev,dat.decday,extrapolate=extrapolate,guess_offset=False)
            dat.constant_space(spacing, min_movement=min_movement)
