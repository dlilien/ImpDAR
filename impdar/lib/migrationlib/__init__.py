#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Migration is a bit messy in order not to be too slow. This should collect the different options.
"""

from .mig_su import migrationSeisUnix
from .mig_python import migrationStolt, migrationPhaseShift, migrationTimeWavenumber

try:
    from .mig_cython import migrationKirchhoff
except ImportError:
    from .mig_python import migrationKirchhoff
