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
    def test_Stolt(self):
        data = NoInitRadarData()
        data = migration_routines.migrationStolt(data)
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR, out_prefix + '_Stolt.mat')
        data.save(out_fn)
        os.remove(out_fn)

    def test_Kirchhoff(self):
        data = NoInitRadarData()
        data = migration_routines.migrationKirchhoff(data)
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR, out_prefix + '_Kirchhoff.mat')
        data.save(out_fn)
        os.remove(out_fn)

    def test_PhaseShiftConstant(self):
        data = NoInitRadarData()
        data = migration_routines.migrationPhaseShift(data)
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR, out_prefix + '_PhaseShiftConstant.mat')
        data.save(out_fn)
        os.remove(out_fn)

    """
    def test_PhaseShiftVariable(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', in_file))
        data = migration_routines.migrationPhaseShift(data,vel_fn='./input_data/velocity_layers.txt')
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR,out_prefix+'_PhaseShiftVariable.mat')
        data.save(out_fn)
    def test_PhaseShiftLateral(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', in_file))
        data = migration_routines.migrationPhaseShift(data,vel_fn='./input_data/velocity_lateral.txt')
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR,out_prefix+'_PhaseShiftLateral.mat')
        data.save(out_fn)
    """


if __name__ == '__main__':
    unittest.main()
