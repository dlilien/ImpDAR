#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Testing loading_utils
"""

import unittest
from impdar.lib.load import loading_utils


class TestLoadGecko(unittest.TestCase):

    def test_common_start(self):
        start = loading_utils.common_start(['abra', 'abracadabra'])
        self.assertEqual('abra', start)

        start = loading_utils.common_start(['abra', 'abra'])
        self.assertEqual('abra', start)

        start = loading_utils.common_start(['abra', 'abra', 'abracad'])
        self.assertEqual('abra', start)

        start = loading_utils.common_start(['abra'])
        self.assertEqual('abra', start)

        start = loading_utils.common_start(['', 'abra'])
        self.assertEqual('', start)


if __name__ == '__main__':
    unittest.main()
