#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Test the methods on the Pick object
"""

import os
import unittest
import numpy as np
from impdar.lib.LastTrace import LastTrace


class TestPickMods(unittest.TestCase):

    def test_mod_line(self):
        lt = LastTrace()
        with self.assertRaises(AttributeError):
            lt.mod_line(0, 1, 1)

        lt.snum = [0]
        lt.tnum = [0]
        with self.assertRaises(ValueError):
            lt.mod_line(1, 50, 40)

        lt.mod_line(0, 50, 40)
        self.assertEqual(lt.snum[0], 50)
        self.assertEqual(lt.tnum[0], 40)

    def test_add_pick(self):
        lt = LastTrace()
        lt.add_pick(0, 10)
        self.assertEqual(len(lt.snum), 1)
        self.assertEqual(len(lt.tnum), 1)
        self.assertEqual(lt.snum[0], 0)
        self.assertEqual(lt.tnum[0], 10)

        lt.add_pick(50, 40)
        self.assertEqual(len(lt.snum), 2)
        self.assertEqual(len(lt.tnum), 2)
        self.assertEqual(lt.snum, [0, 50])
        self.assertEqual(lt.tnum, [10, 40])

        with self.assertRaises(TypeError):
            lt.add_pick([12, 15.5], 0)


if __name__ == '__main__':
    unittest.main()
