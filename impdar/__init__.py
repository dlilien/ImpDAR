#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3 license.

"""Skeleton import impdar."""

import matplotlib

# We are going to do this here so that it comes first.
# We need Qt5 for the picker
# but we can break things with no-qt installations that way
try:
    import PyQt5
    matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib.use('AGG')
    import matplotlib.pyplot as plt
