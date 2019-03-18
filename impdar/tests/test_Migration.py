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
out_prefix = 'rectangle'
in_file = out_prefix+'_gprMax_Bscan.h5'

class TestMigration(unittest.TestCase):
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

    def test_GazdagConstant(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', in_file))
        data = migration_routines.migrationGazdag(data)
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR,out_prefix+'_GazdagConstant.mat')
        data.save(out_fn)

    def test_GazdagVariable(self):
        data = load_gprMax.load_gprMax(os.path.join(THIS_DIR, 'input_data', in_file))
        data = migration_routines.migrationGazdag(data,vel_fn='./input_data/velocity_layers.txt')
        if not os.path.isdir(OUT_DIR):
            os.mkdir(OUT_DIR)
        out_fn = os.path.join(OUT_DIR,out_prefix+'_GazdagVariable.mat')
        data.save(out_fn)

if __name__ == '__main__':
    unittest.main()
