#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""Methods for saving radar data in different formats."""
from ..gpslib import get_conversion
import numpy as np
from scipy.io import savemat
from ..RadarFlags import RadarFlags

# Try to enable saving to shapefiles
try:
    from osgeo import osr, ogr
    CONVERSIONS_ENABLED = True
except ImportError:
    CONVERSIONS_ENABLED = False

# Try to enable saving to segy
try:
    import segyio
    SEGY = True
except ImportError:
    SEGY = False


def save(self, fn):
    """Save the radar data.

    Parameters
    ----------
    fn: str
        Filename. Should have a .mat extension
    """
    mat = {}

    for attr in self.attrs_guaranteed:
        if getattr(self, attr) is not None:
            mat[attr] = getattr(self, attr)
        else:
            # this guards against error in matlab format
            mat[attr] = 0
    for attr in self.attrs_optional + self.stodeep_attrs:
        if hasattr(self, attr) and getattr(self, attr) is not None:
            mat[attr] = getattr(self, attr)
    if hasattr(self, 'picks') and self.picks is not None:
        mat['picks'] = self.picks.to_struct()
    if self.flags is not None:
        mat['flags'] = self.flags.to_matlab()
    else:
        # We want the structure available to prevent read errors from corrupt files
        mat['flags'] = RadarFlags().to_matlab()

    # Make sure not to expand the size of the data due to type conversion
    if hasattr(self, 'data_dtype') and (
            self.data_dtype is not None) and (self.data_dtype != mat['data'].dtype):
        # Be careful of obliterating NaNs
        # We will use singles instead of ints for this guess
        if (self.data_dtype in [int, np.int8, np.int16]) and np.any(np.isnan(mat['data'])):
            print('Warning: new file is float16 rather than ',
                  self.data_dtype, ' since we now have NaNs')
            mat['data'] = mat['data'].astype(np.float16)
        elif (self.data_dtype in [np.int32]) and np.any(np.isnan(mat['data'])):
            print('Warning: new file is float32 rather than ',
                  self.data_dtype, ' since we now have NaNs')
            mat['data'] = mat['data'].astype(np.float32)
        elif (self.data_dtype in [np.int64]) and np.any(np.isnan(mat['data'])):
            print('Warning: new file is float64 rather than ',
                  self.data_dtype, ' since we now have NaNs')
            mat['data'] = mat['data'].astype(np.float64)
        else:
            mat['data'] = mat['data'].astype(self.data_dtype)
    savemat(fn, mat)


def save_as_segy(self, fn):
    """Save as a (non standard-compliant) SEGY file that can be used with, e.g., SeisUNIX.

    Parameters
    ----------
    fn: str
        The filename to save as.

    Raises
    ------
    ImportError
        If segyio cannot be imported.
    """
    if not SEGY:
        raise ImportError('segyio failed to import, cannot save as segy')

    segyio.tools.from_array2D(fn,
                              np.ascontiguousarray(self.data.transpose(), np.float32),
                              dt=self.dt * 1.0e12)


def output_shp(self, fn, t_srs=None, target_out=None):
    """Output a shapefile of the traces.

    If there are any picks, we want to output these.
    If not, we will only output the tracenumber.
    This function requires osr/gdal for shapefile creation.
    I suggest exporting a csv if you don't want to deal with gdal.

    Parameters
    ----------
    fn: str
        The filename of the output
    t_srs: int, optional
        EPSG number of the target spatial reference system. Default 4326 (wgs84)
    target_out: str, optional
        Used to overwrite the default output format of picks.
        By default, try to write depth and if there is no nmo_depth use TWTT.
        You might want to use this to get the output in TWTT or sample number
        (options are depth, elev, twtt, snum)

    Raises
    ------
    ImportError
        If osgeo cannot be imported
    """
    if not CONVERSIONS_ENABLED:
        raise ImportError('osgeo could not be imported')

    if t_srs is not None:
        # We overwrite the t_srs with the WKT version
        cT, t_srs = get_conversion(t_srs=t_srs)
        pts = np.array(cT(np.vstack((self.long, self.lat)).transpose()))
    else:
        pts = np.vstack((self.long, self.lat)).transpose()
        t_srs = 'EPSG:4326'

    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.CreateDataSource(fn)
    out_srs = osr.SpatialReference()
    out_srs.SetFromUserInput(t_srs)
    layer = data_source.CreateLayer('traces', out_srs, ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn('TraceNum', ogr.OFTInteger))

    if self.picks is not None and self.picks.samp2 is not None:
        out_name, target_out_array = self._get_pick_targ_info(target_out)
        for picknum in self.picks.picknums:
            layer.CreateField(ogr.FieldDefn('L{:d}_{:s}'.format(picknum, out_name), ogr.OFTReal))

    # Process the text file and add the attributes and features to the shapefile
    for trace in range(self.tnum):
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetField('TraceNum', trace + 1)
        if self.picks is not None and self.picks.samp2 is not None:
            for i, picknum in enumerate(self.picks.picknums):
                if out_name != 'elev':
                    if not np.isnan(self.picks.samp2[i, trace]):
                        feature.SetField('L{:d}_{:s}'.format(picknum, out_name),
                                         target_out_array[int(self.picks.samp2[i, trace])])
                    else:
                        feature.SetField('L{:d}_{:s}'.format(picknum, out_name), np.nan)

                else:
                    if not np.isnan(self.picks.samp2[i, trace]):
                        feature.SetField('L{:d}_{:s}'.format(picknum, out_name),
                                         self.elev[trace] - target_out_array[
                                             int(self.picks.samp2[i, trace])])
                    else:
                        feature.SetField('L{:d}_{:s}'.format(picknum, out_name), np.nan)

        x, y = pts[trace, 0], pts[trace, 1]
        wkt = 'POINT({:f} {:f})'.format(x, y)
        point = ogr.CreateGeometryFromWkt(wkt)
        feature.SetGeometry(point)
        layer.CreateFeature(feature)
        feature = None
    data_source = None


def output_csv(self, fn, target_out=None, delimiter=','):
    """Output a csv of the traces.

    If there are any picks, we want to output these.
    If not, we will only output the tracenumber.

    Parameters
    ----------
    fn: str
        The filename of the output
    target_out: str, optional
        Used to overwrite the default output format of picks.
        By default, try to write depth and if there is no nmo_depth use TWTT.
        You might want to use this to get the output in TWTT or sample number
        (options are depth, elev, twtt, snum)
    delimiter: str, optional
        Delimiter for csv. Default ','.
    """
    header = delimiter.join(['lat', 'lon', 'tnum'])
    outs = np.vstack((self.lat, self.long, np.arange(self.tnum) + 1))

    if self.picks is not None and self.picks.samp2 is not None:
        out_name, target_out_array = self._get_pick_targ_info(target_out)
        for picknum in self.picks.picknums:
            header += (delimiter + 'Layer_{:d}_{:s}'.format(picknum, out_name))

        out_ind_picks = self.picks.samp2.astype(int)
        out_ind_picks_viableind = out_ind_picks.copy()
        out_ind_picks_viableind[out_ind_picks_viableind < 0] = 0
        out_arr_picks = target_out_array[out_ind_picks_viableind]
        out_arr_picks[out_ind_picks < 0] = np.nan
        outs = np.vstack((outs, out_arr_picks))
    np.savetxt(fn, outs.transpose(), header=header, delimiter=delimiter)


def _get_pick_targ_info(self, target_out):
    """Get the rate type of pick information for returning.

    Use this code to eliminate some overlap in exporting csvs and shps.

    Parameters
    ----------
    target_out: str
        The z coordinate to output. Options are depth, twtt, elev, snum.

    Raises
    ------
    ValueError
        If the target output is bad.
    """
    if target_out is None:
        if self.nmo_depth is not None:
            out_name = 'depth'
            target_out_array = self.nmo_depth
        else:
            out_name = 'twtt'
            target_out_array = self.travel_time
    else:
        out_name = target_out
        if target_out == 'depth':
            if (not hasattr(self, 'nmo_depth')) or self.nmo_depth is None:
                raise AttributeError('Cannot do depth output with no nmo_depth')
            target_out_array = self.nmo_depth
        elif target_out == 'elev':
            if (not hasattr(self, 'elev')) or self.elev is None:
                raise AttributeError('Cannot do depth output with no nmo_depth')
            target_out_array = self.nmo_depth
        elif target_out == 'twtt':
            target_out_array = self.travel_time
        elif target_out == 'snum':
            target_out_array = np.arange(self.snum)
        else:
            raise ValueError('target_out must be snum, twtt, depth, or elev')
    return out_name, target_out_array
