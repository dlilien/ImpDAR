#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL-3.0 license.

"""
Test the basics of RadarData
"""

import unittest
import os
from impdar.lib import migration_routines, load_gprMax

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(THIS_DIR,'Migration_tests')

class TestMigration(unittest.TestCase):
    def test_Stolt(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', 'rectangle_gprMax_Bscan.h5'))
        data = migration_routines.migrationStolt(data)
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR,'rectangle_Stolt.mat')
        data.save(out_fn)
    """
    def test_Kirchhoff(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', 'rectangle_gprMax_Bscan.h5'))
        data = migration_routines.migrationKirchhoff(data)

    def test_GazdagConstant(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', 'rectangle_gprMax_Bscan.h5'))
        data = migration_routines.migrationGazdag(data)
    def test_GazdagVariable(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', 'rectangle_gprMax_Bscan.h5'))
        data = migration_routines.migrationGazdag(data,vel_fn='./input_data/velocity_layers.txt')
    """
if __name__ == '__main__':
    unittest.main()
