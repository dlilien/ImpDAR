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

# We are going to do this here so that it comes first. We need Qt5 for the picker
# but we can break things with no-qt installations that way
matplotlib.use('Qt5Agg')
try:
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

from .lib import load, process, plot, convert
from .lib.RadarData import RadarData
