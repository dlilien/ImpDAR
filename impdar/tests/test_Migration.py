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

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(THIS_DIR,'Migration_tests')
out_prefix = 'rectangle'
in_file = out_prefix+'_gprMax_Bscan.h5'

class TestMigration(unittest.TestCase):
    """
    def test_Stolt(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', in_file))
        data = migration_routines.migrationStolt(data)
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR,out_prefix+'_Stolt.mat')
        data.save(out_fn)

    def test_Kirchhoff(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', in_file))
        data = migration_routines.migrationKirchhoff(data)
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR,out_prefix+'_Kirchhoff.mat')
        data.save(out_fn)

    def test_PhaseShiftConstant(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', in_file))
        data = migration_routines.migrationPhaseShift(data)
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR,out_prefix+'_PhaseShiftConstant.mat')
        data.save(out_fn)

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
    pass


if __name__ == '__main__':
    unittest.main()
