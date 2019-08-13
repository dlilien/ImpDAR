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

import os
import unittest
import subprocess as sp
import numpy as np

from impdar.lib import migrationlib
from impdar.lib.migrationlib import mig_python
from impdar.lib.load import load_segy
from impdar.lib.NoInitRadarData import NoInitRadarData

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(THIS_DIR, 'Migration_tests')
OUT_PREFIX = 'rectangle'
# in_file = out_prefix+'_gprMax_Bscan.h5'


class TestMigration(unittest.TestCase):

    def test_check_data_shape(self):
        data = NoInitRadarData(big=True)

        # should pass, i.e. nothing happens
        mig_python._check_data_shape(data)

        # make it fail
        data.data = np.ones((1, 1))
        with self.assertRaises(ValueError):
            mig_python._check_data_shape(data)

    def test_getVelocityProfile(self):
        data = NoInitRadarData(big=True)
        self.assertEqual(1.68e8, mig_python.getVelocityProfile(data, 1.68e8))

        # need reasonable input here for 2d. Needs a different travel time.
        data.travel_time = data.travel_time / 10.
        mig_python.getVelocityProfile(data, np.genfromtxt(os.path.join(THIS_DIR, 'input_data', 'velocity_layers.txt')))

        # this should still work since we are close
        data.travel_time = data.travel_time / 10.
        twod = np.genfromtxt(os.path.join(THIS_DIR, 'input_data', 'velocity_layers.txt'))
        twod = twod * 0.0045 + 1.0e-7 * twod[1]
        mig_python.getVelocityProfile(data, twod)

        # need reasonable input here for 3d
        data = NoInitRadarData(big=True)
        mig_python.getVelocityProfile(data, np.genfromtxt(os.path.join(THIS_DIR, 'input_data', 'velocity_lateral.txt')))

        # Bad distance with good 3d grid
        data.dist = None
        with self.assertRaises(ValueError):
            mig_python.getVelocityProfile(data, np.genfromtxt(os.path.join(THIS_DIR, 'input_data', 'velocity_lateral.txt')))
        data = NoInitRadarData(big=True)

        # this should fail on bad z
        twod_vel = 1.68e8 * np.ones((10, 2))
        twod_vel[:, 1] = 0.
        with self.assertRaises(ValueError):
            mig_python.getVelocityProfile(data, twod_vel)

        # Use some bad x values
        with self.assertRaises(ValueError):
            mig_python.getVelocityProfile(data, 1.68e8 * np.ones((10, 3)))
        # bad z values
        threed_vel = 1.68e8 * np.ones((10, 3))
        threed_vel[:, -1] = np.arange(10) * 1000.
        with self.assertRaises(ValueError):
            mig_python.getVelocityProfile(data, threed_vel)

        # Make sure we reject bad input shapes
        with self.assertRaises(ValueError):
            mig_python.getVelocityProfile(data, 1.68e8 * np.ones((8,)))
        with self.assertRaises(ValueError):
            mig_python.getVelocityProfile(data, 1.68e8 * np.ones((8, 1)))
        with self.assertRaises(ValueError):
            mig_python.getVelocityProfile(data, 1.68e8 * np.ones((1, 2)))
        with self.assertRaises(ValueError):
            mig_python.getVelocityProfile(data, 1.68e8 * np.ones((8, 4)))

    def test_Stolt(self):
        data = NoInitRadarData(big=True)
        data = migrationlib.migrationStolt(data)

    def test_Kirchhoff(self):
        data = NoInitRadarData(big=True)
        data = migrationlib.migrationKirchhoff(data)

    def test_TimeWavenumber(self):
        data = NoInitRadarData(big=True)
        data = migrationlib.migrationTimeWavenumber(data)

    def test_PhaseShiftConstant(self):
        data = NoInitRadarData(big=True)
        data = migrationlib.migrationPhaseShift(data)

    def test_PhaseShiftVariable(self):
        data = NoInitRadarData(big=True)
        data.travel_time = data.travel_time / 10.
        data = migrationlib.migrationPhaseShift(data, vel_fn=os.path.join(THIS_DIR, 'input_data', 'velocity_layers.txt'))

        data = NoInitRadarData(big=True)
        with self.assertRaises(TypeError):
            data = migrationlib.migrationPhaseShift(data, vel_fn=os.path.join(THIS_DIR, 'input_data', 'notafile.txt'))

    def test_PhaseShiftLateral(self):
        data = NoInitRadarData(big=True)
        data = migrationlib.migrationPhaseShift(data, vel_fn=os.path.join(THIS_DIR, 'input_data', 'velocity_lateral.txt'))

    @unittest.skipIf(sp.Popen(['which', 'sumigtk']).wait() != 0 or (not load_segy.SEGY), 'SeisUnix or segy not found')
    def test_sumigtk(self):
        data = NoInitRadarData(big=True)
        data.dt = 1.0e-9
        data.travel_time = data.travel_time * 1.0e-9
        data.fn = os.path.join(THIS_DIR, 'input_data', 'rectangle_sumigtk.mat')
        migrationlib.migrationSeisUnix(data)

    @unittest.skipIf(sp.Popen(['which', 'sustolt']).wait() != 0 or (not load_segy.SEGY), 'SeisUnix or segy not found')
    def test_sustolt(self):
        data = NoInitRadarData(big=True)
        data.dt = 1.0e-9
        data.travel_time = data.travel_time * 1.0e-9
        data.fn = os.path.join(THIS_DIR, 'input_data', 'rectangle_sustolt.mat')
        migrationlib.migrationSeisUnix(data)

    @unittest.skipIf(sp.Popen(['which', 'sustolt']).wait() != 0 or load_segy.SEGY, 'Test of edge case')
    def test_sustolt_nosegy(self):
        data = NoInitRadarData(big=True)
        data.dt = 1.0e-9
        data.travel_time = data.travel_time * 1.0e-9
        data.fn = os.path.join(THIS_DIR, 'input_data', 'rectangle_sustolt.mat')
        with self.assertRaises(ImportError):
            migrationlib.migrationSeisUnix(data)

    @unittest.skipIf(sp.Popen(['which', 'sustolt']).wait() == 0, 'Test for no SeisUnix')
    def test_sustolt_seisunix(self):
        data = NoInitRadarData(big=True)
        data.dt = 1.0e-9
        data.travel_time = data.travel_time * 1.0e-9
        data.fn = os.path.join(THIS_DIR, 'input_data', 'rectangle_sustolt.mat')
        with self.assertRaises(FileNotFoundError):
            migrationlib.migrationSeisUnix(data)

    def tearDown(self):
        for suff in ['PhaseShiftLateral', 'PhaseShiftConstant', 'PhaseShiftVariable', 'Kirchoff', 'Stolt', 'sumigtk', 'sustolt']:
            if os.path.exists(os.path.join(THIS_DIR, 'input_data', 'rectangle_' + suff + '.mat')):
                os.remove(os.path.join(THIS_DIR, 'input_data', 'rectangle_' + suff + '.mat'))


if __name__ == '__main__':
    unittest.main()
