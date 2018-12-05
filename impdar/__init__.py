#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3 license.

"""
Skeleton import impdar
"""

import platform
import matplotlib
# if platform.system() == 'Darwin':
#     matplotlib.use('macosx')
# else:
#     matplotlib.use('gtk3agg')
matplotlib.use('Qt5Agg')

from .lib import load, process, plot, convert
from .lib.RadarData import RadarData
