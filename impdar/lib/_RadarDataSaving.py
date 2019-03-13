#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
A class that RadarData inherits so that it has convenient saving methods
"""

from scipy.io import savemat
from .RadarFlags import RadarFlags

# Try to enable saving to shapefiles
try:
    from osgeo import osr, ogr
    conversions_enabled = True
except ImportError:
    conversions_enabled = False

# Try to enable saving to segy
try:
    import segyio
    segy = True
except ImportError:
    segy = False


class RadarDataSaving:

    def save(self, fn):
        """Save the radar data

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
        for attr in self.attrs_optional:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
        if hasattr(self, 'picks') and self.picks is not None:
            mat['picks'] = self.picks.to_struct()
        if self.flags is not None:
            mat['flags'] = self.flags.to_matlab()
        else:
            # We want the structure available to prevent read errors from corrupt files
            mat['flags'] = RadarFlags().to_matlab()
        savemat(fn, mat)

    def save_as_segy(self, fn):
        if not segy:
            raise ImportError('segyio failed to import, cannot save as segy')

        spec = segyio.spec()
        spec.sorting = 2
        spec.format = 1
        spec.samples = range(self.snum)
        spec.tracecount = self.tnum
        spec.ilines = range(self.tnum)
        
        # We assume that this is radar data on a line, so there is no cross-line
        spec.xlines = [0]
        with segyio.create(fn, spec) as f:
            for il in spec.ilines:
                f.header[il] = {segyio.su.offset: 1, segyio.su.iline: il, segyio.su.xline: 0}
                f.trace[il] = self.data[:, il]
                f.bin.update(tsort=segyio.TraceSortingFormat.INLINE_SORTING)

    def output_shp(self, fn, t_srs=4326):
        if not conversions_enabled:
            raise ImportError('osgeo was not imported')
        out_srs = osr.SpatialReference()
        out_srs.ImportFromEPSG(t_srs)
        in_srs = osr.SpatialReference()
        in_srs.ImportFromEPSG(4326)
        cT = osr.CoordinateTransformation(in_srs, out_srs)

        driver = ogr.GetDriverByName('ESRI Shapefile')
        data_source = driver.CreateDataSource(fn)
        layer = data_source.CreateLayer('traces', out_srs, ogr.wkbPoint)
        layer.CreateField(ogr.FieldDefn('TraceNum', ogr.OFTInteger))
        if hasattr(self, 'layers') and self.layers is not None:
            for i in range(len(self.layers)):
                layer.CreateField(ogr.FieldDefn('Layer {:d}'.format(i), ogr.OFTReal))

        # Process the text file and add the attributes and features to the shapefile
        for trace in range(self.tnum):
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetField('TraceNum', trace)
            if hasattr(self, 'layers') and self.layers is not None:
                for i in range(len(self.layers)):
                    feature.SetField('Layer {:d}'.format(i), self.layers[i][trace])
            x, y, _ = cT.TransformPoint(self.long[trace], self.lat[trace])        
            wkt = 'POINT({:f} {:f})'.format(x, y)
            point = ogr.CreateGeometryFromWkt(wkt)
            feature.SetGeometry(point)
            layer.CreateFeature(feature)
            feature = None
        data_source = None
