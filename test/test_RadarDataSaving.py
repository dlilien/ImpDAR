#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""

"""
import warnings
import os
import unittest
import numpy as np
from impdar.lib.RadarData import RadarData
from impdar.lib.NoInitRadarData import NoInitRadarData
from impdar.lib.RadarData._RadarDataSaving import CONVERSIONS_ENABLED
from impdar.lib.RadarFlags import RadarFlags
from impdar.lib.Picks import Picks

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestRadarDataSaving(unittest.TestCase):

    def test_WriteNoFLags(self):
        rd = NoInitRadarData()
        rd.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        os.remove(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    def testWriteWithFlags(self):
        rd = NoInitRadarData()
        rd.flags = RadarFlags()
        rd.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    def testWriteWithPicksBlank(self):
        rd = NoInitRadarData()
        rd.picks = Picks(rd)
        rd.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        self.assertTrue(data.picks is not None)
        self.assertTrue(data.picks.lasttrace is not None)
        self.assertTrue(data.picks.lasttrace.tnum is None)
        self.assertTrue(data.picks.samp1 is None)
        self.assertTrue(data.picks.samp2 is None)
        self.assertTrue(data.picks.samp3 is None)

    def testWriteWithPicksFull(self):
        rd = NoInitRadarData()
        rd.picks = Picks(rd)
        rd.picks.add_pick()
        rd.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        data = RadarData(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        self.assertTrue(data.picks is not None)
        self.assertTrue(data.picks.lasttrace is not None)
        self.assertTrue(data.picks.samp1 is not None)
        self.assertTrue(data.picks.samp2 is not None)
        self.assertTrue(data.picks.samp3 is not None)

    def test_WriteRead(self):
        # We are going to create a really bad file (most info missing) and see if we recover it or get an error
        rd = NoInitRadarData()
        rd.save(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))
        RadarData(os.path.join(THIS_DIR, 'input_data', 'test_out.mat'))

    def tearDown(self):
        for fn in ['test_out.mat', 'test.shp', 'test.shx', 'test.prj', 'test.dbf']:
            if os.path.exists(os.path.join(THIS_DIR, 'input_data', fn)):
                os.remove(os.path.join(THIS_DIR, 'input_data', fn))


class TestRadarDataExports(unittest.TestCase):

    def test__get_pick_targ_infoAutoselect(self):
        # Make sure that we are selecting the proper output format
        rd = NoInitRadarData()

        # With no depth, we should output travel time
        rd.nmo_depth = None
        out_name, tout = rd._get_pick_targ_info(None)
        self.assertEqual(out_name, 'twtt')
        self.assertTrue(np.all(tout == rd.travel_time))

        # With depth, return depth
        rd.nmo_depth = np.arange(len(rd.travel_time)) * 1.1
        out_name, tout = rd._get_pick_targ_info(None)
        self.assertEqual(out_name, 'depth')
        self.assertTrue(np.all(tout == rd.nmo_depth))

    def test__get_pick_targ_infoBadSelections(self):
        # Make sure that we are selecting the proper output format
        rd = NoInitRadarData()

        # Try depth when there is no depth
        rd.nmo_depth = None
        with self.assertRaises(AttributeError):
            out_name, tout = rd._get_pick_targ_info('depth')

        # Elevation with no depth or elevation
        with self.assertRaises(AttributeError):
            out_name, tout = rd._get_pick_targ_info('elev')

        # Elevation with depth but not elevation
        rd.nmo_depth = np.arange(len(rd.travel_time)) * 1.1
        with self.assertRaises(AttributeError):
            out_name, tout = rd._get_pick_targ_info('elev')

        # Now try to pass a bad value for the selection
        with self.assertRaises(ValueError):
            out_name, tout = rd._get_pick_targ_info('dummy')
        with self.assertRaises(ValueError):
            out_name, tout = rd._get_pick_targ_info(['dummy', 'snum'])

    def test__get_pick_targ_infoGoodSelections(self):
        # Make sure that we are selecting the proper output format
        rd = NoInitRadarData()
        rd.nmo_depth = np.arange(len(rd.travel_time)) * 1.1
        rd.elev = np.arange(rd.tnum) * 1001

        out_name, tout = rd._get_pick_targ_info('twtt')
        self.assertEqual(out_name, 'twtt')
        self.assertTrue(np.all(tout == rd.travel_time))

        out_name, tout = rd._get_pick_targ_info('depth')
        self.assertEqual(out_name, 'depth')
        self.assertTrue(np.all(tout == rd.nmo_depth))

        out_name, tout = rd._get_pick_targ_info('snum')
        self.assertEqual(out_name, 'snum')
        self.assertTrue(np.all(tout == np.arange(rd.snum)))

        out_name, tout = rd._get_pick_targ_info('elev')
        self.assertEqual(out_name, 'elev')

    @unittest.skipIf(not CONVERSIONS_ENABLED, 'No GDAL on this version')
    def test_output_shp_nolayers(self):
        rd = NoInitRadarData()
        with warnings.catch_warnings(record=True) as w:
            rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test.shp'))
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertIn("deprecated", str(w[-1].message))

    @unittest.skipIf(not CONVERSIONS_ENABLED, 'No GDAL on this version')
    def test_output_shp_picks(self):

        # Make sure that we are selecting the proper output format
        rd = NoInitRadarData()
        rd.nmo_depth = np.arange(len(rd.travel_time)) * 1.1
        rd.elev = np.arange(rd.tnum) * 1001
        rd.picks = Picks(rd)
        rd.picks.add_pick()

        # First, export with NaNs, both with normal field (depth) and elev
        rd.picks.samp2[:] = np.nan
        with warnings.catch_warnings(record=True) as w:
            rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test0.shp'))
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertIn("deprecated", str(w[-1].message))
        with warnings.catch_warnings(record=True) as w:
            rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test1.shp'), target_out='elev')
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertIn("deprecated", str(w[-1].message))

        # Fill in NaNs
        rd.picks.samp2[:] = 1
        with warnings.catch_warnings(record=True) as w:
            rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test2.shp'))
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertIn("deprecated", str(w[-1].message))
        with warnings.catch_warnings(record=True) as w:
            rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test3.shp'), target_out='elev')
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertIn("deprecated", str(w[-1].message))

        # Check geometry
        with warnings.catch_warnings(record=True) as w:
            rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test4.shp'), t_srs='EPSG:3413')
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertIn("deprecated", str(w[-1].message))

    @unittest.skipIf(CONVERSIONS_ENABLED, 'Version has GDAL, just checking we fail without')
    def test_output_shp_nolayers_nogdal(self):
        rd = NoInitRadarData()
        with self.assertRaises(ImportError):
            rd.output_shp(os.path.join(THIS_DIR, 'input_data', 'test.shp'))

    def test_output_csv(self):
        # Make sure that we are selecting the proper output format
        rd = NoInitRadarData()
        rd.nmo_depth = np.arange(len(rd.travel_time)) * 1.1
        rd.elev = np.arange(rd.tnum) * 1001
        rd.picks = Picks(rd)
        rd.picks.add_pick()

        # First, export with NaNs
        rd.picks.samp2[:] = np.nan
        rd.output_csv(os.path.join(THIS_DIR, 'input_data', 'test.csv'))
        with open(os.path.join(THIS_DIR, 'input_data', 'test.csv')) as fin:
            lines = fin.readlines()
            # we should have four entries: lat, lon, trace, and the one pick in header and data
            self.assertEqual(len(lines[0].split(',')), 4)
            self.assertEqual(len(lines[1].split(',')), 4)

            # we should have a row per trace, plus a header
            self.assertEqual(len(lines), rd.tnum + 1)

            # The final header should be in terms of depth
            self.assertTrue(lines[0].index('depth') > 0)

        # Fill in NaNs
        rd.picks.samp2[:] = 1
        rd.output_csv(os.path.join(THIS_DIR, 'input_data', 'test.csv'))
        with open(os.path.join(THIS_DIR, 'input_data', 'test.csv')) as fin:
            lines = fin.readlines()
            # we should have four entries: lat, lon, trace, and the one pick in header and data
            self.assertEqual(len(lines[0].split(',')), 4)
            self.assertEqual(len(lines[1].split(',')), 4)

            # we should have a row per trace, plus a header
            self.assertEqual(len(lines), rd.tnum + 1)

            # The final header should be in terms of depth
            self.assertTrue(lines[0].index('depth') > 0)

        # Check output target for elevation, which is the only weird one
        rd.output_csv(os.path.join(THIS_DIR, 'input_data', 'test.csv'), target_out='elev')
        with open(os.path.join(THIS_DIR, 'input_data', 'test.csv')) as fin:
            lines = fin.readlines()
            # we should have four entries: lat, lon, trace, and the one pick in header and data
            self.assertEqual(len(lines[0].split(',')), 4)
            self.assertEqual(len(lines[1].split(',')), 4)

            # we should have a row per trace, plus a header
            self.assertEqual(len(lines), rd.tnum + 1)

            # The final header should be in terms of elev
            self.assertTrue(lines[0].index('elev') > 0)

    def test_output_csv_nolayers(self):
        rd = NoInitRadarData()
        rd.output_csv(os.path.join(THIS_DIR, 'input_data', 'test.csv'))
        with open(os.path.join(THIS_DIR, 'input_data', 'test.csv')) as fin:
            lines = fin.readlines()
            # we should only have three entries: lat, lon, trace in header and data
            self.assertEqual(len(lines[0].split(',')), 3)
            self.assertEqual(len(lines[1].split(',')), 3)

            # we should have a row per trace, plus a header
            self.assertEqual(len(lines), rd.tnum + 1)

    def tearDown(self):
        for i in range(6):
            for fn in ['test_out.mat', 'test{:d}.shp'.format(i), 'test{:d}.shx'.format(i), 'test{:d}.prj'.format(i), 'test{:d}.dbf'.format(i), 'test.csv']:
                if os.path.exists(os.path.join(THIS_DIR, 'input_data', fn)):
                    os.remove(os.path.join(THIS_DIR, 'input_data', fn))


if __name__ == '__main__':
    unittest.main()
