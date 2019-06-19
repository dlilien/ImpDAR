#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Test the migration routines

Author:
Benjamin Hills
benjaminhhills@gmail.com
University of Washington
Earth and Space Sciences

Mar 12 2019

"""

import unittest
import os
from impdar.lib import migration_routines, load_gprMax
from impdar.lib.RadarData import RadarData
import numpy as np

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(THIS_DIR, 'Migration_tests')
out_prefix = 'rectangle'
# in_file = out_prefix+'_gprMax_Bscan.h5'


class NoInitRadarData(RadarData):
    # This only exists so we can do tests on writing without reading

    def __init__(self):
        self.data = np.zeros((10, 20))
        # need to set this to avoid divide by zero later
        self.dt = 1
        self.dist = np.arange(20)
        self.tnum = 20
        self.trace_num = np.arange(self.tnum) + 1.
        self.snum = 10
        self.travel_time = np.arange(10)
        self.trace_int = 1


class TestMigration(unittest.TestCase):

    def test_check_data_shape(self):
        data = NoInitRadarData()

        # should pass, i.e. nothing happens
        migration_routines._check_data_shape(data)

        # make it fail
        data.data = np.ones((1, 1))
        with self.assertRaises(ValueError):
            migration_routines._check_data_shape(data)

    def test_getVelocityProfile(self):
        data = NoInitRadarData()
        self.assertEqual(1.68e8, migration_routines.getVelocityProfile(data, 1.68e8))

        # need reasonable input here for 2d. Needs a different travel time.
        data.travel_time = data.travel_time / 10.
        migration_routines.getVelocityProfile(data, np.genfromtxt(os.path.join(THIS_DIR, 'input_data', 'velocity_layers.txt')))

        # this should still work since we are close
        data.travel_time = data.travel_time / 10.
        twod = np.genfromtxt(os.path.join(THIS_DIR, 'input_data', 'velocity_layers.txt'))
        twod = twod * 0.0045 + 1.0e-7 * twod[1]
        migration_routines.getVelocityProfile(data, twod)

        # need reasonable input here for 3d
        data = NoInitRadarData()
        migration_routines.getVelocityProfile(data, np.genfromtxt(os.path.join(THIS_DIR, 'input_data', 'velocity_lateral.txt')))

        # Bad distance with good 3d grid
        data.dist = None
        with self.assertRaises(ValueError):
            migration_routines.getVelocityProfile(data, np.genfromtxt(os.path.join(THIS_DIR, 'input_data', 'velocity_lateral.txt')))
        data = NoInitRadarData()

        # this should fail on bad z
        twod_vel = 1.68e8 * np.ones((10, 2))
        twod_vel[:, 1] = 0.
        with self.assertRaises(ValueError):
            migration_routines.getVelocityProfile(data, twod_vel)
        
        # Use some bad x values
        with self.assertRaises(ValueError):
            migration_routines.getVelocityProfile(data, 1.68e8 * np.ones((10, 3)))
        # bad z values
        threed_vel = 1.68e8 * np.ones((10, 3))
        threed_vel[:, -1] = np.arange(10) * 1000.
        with self.assertRaises(ValueError):
            migration_routines.getVelocityProfile(data, threed_vel)

        # Make sure we reject bad input shapes
        with self.assertRaises(ValueError):
            migration_routines.getVelocityProfile(data, 1.68e8 * np.ones((8,)))
        with self.assertRaises(ValueError):
            migration_routines.getVelocityProfile(data, 1.68e8 * np.ones((8, 1)))
        with self.assertRaises(ValueError):
            migration_routines.getVelocityProfile(data, 1.68e8 * np.ones((1, 2)))
        with self.assertRaises(ValueError):
            migration_routines.getVelocityProfile(data, 1.68e8 * np.ones((8, 4)))
        
    def test_Stolt(self):
        data = NoInitRadarData()
        data = migration_routines.migrationStolt(data)

    def test_Kirchhoff(self):
        data = NoInitRadarData()
        data = migration_routines.migrationKirchhoff(data)

    def test_PhaseShiftConstant(self):
        data = NoInitRadarData()
        data = migration_routines.migrationPhaseShift(data)

    def test_PhaseShiftVariable(self):
        data = NoInitRadarData()
        data.travel_time = data.travel_time / 10.
        data = migration_routines.migrationPhaseShift(data, vel_fn=os.path.join(THIS_DIR, 'input_data', 'velocity_layers.txt'))

        data = NoInitRadarData()
        with self.assertRaises(TypeError):
            data = migration_routines.migrationPhaseShift(data, vel_fn=os.path.join(THIS_DIR, 'input_data', 'notafile.txt'))

    def test_PhaseShiftLateral(self):
        data = NoInitRadarData()
        data = migration_routines.migrationPhaseShift(data, vel_fn=os.path.join(THIS_DIR, 'input_data', 'velocity_lateral.txt'))


if __name__ == '__main__':
    unittest.main()
